"""Conservative structural audit and input checks; not a minimization pipeline."""
import csv,json,math,hashlib,re
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

K_RESTRAINT=4184.0 # kJ mol^-1 nm^-2 for 0.5*k*r^2, 10 kcal mol^-1 A^-2
HEADERS=['variant_id','group','rna_sequence','modifications_json','replicate_id','notes']
RNA={'A','C','G','U'}
SUGAR={"O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"}
BASE={'A':{'N9','C8','N7','C5','C6','N6','N1','C2','N3','C4'}, 'G':{'N9','C8','N7','C5','C6','O6','N1','C2','N2','N3','C4'},'C':{'N1','C2','O2','N3','C4','N4','C5','C6'},'U':{'N1','C2','O2','N3','C4','O4','C5','C6'}}

def sha256(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def validate_variants(path,expected_counts={'benchmark':34,'control':27}):
 with open(path,newline='',encoding='utf-8-sig') as f:
  rd=csv.DictReader(f)
  if rd.fieldnames!=HEADERS: raise ValueError('CSV headers must exactly equal '+','.join(HEADERS))
  rows=list(rd)
 seen=set()
 for i,r in enumerate(rows,2):
  if None in r or any(v is None for v in r.values()): raise ValueError(f'Row {i}: incorrect number of fields')
  v=r['variant_id'].strip(); group=r['group']; seq=r['rna_sequence'].strip().upper()
  if not re.fullmatch(r'[A-Za-z0-9_.-]+',v) or v in ('.','..') or v in seen: raise ValueError(f'Row {i}: invalid or duplicate variant_id')
  seen.add(v)
  if group not in ('benchmark','control') or not r['replicate_id'].strip(): raise ValueError(f'Row {i}: invalid group or replicate ID')
  if len(seq)!=31 or set(seq)-RNA: raise ValueError(f'Row {i}: sequence must contain exactly 31 A/C/G/U characters; variable-length structures require a separate validated builder')
  try: mods=json.loads(r['modifications_json'])
  except json.JSONDecodeError as e: raise ValueError(f'Row {i}: malformed modifications JSON: {e}') from e
  if not isinstance(mods,list): raise ValueError(f'Row {i}: modifications must be an array')
  used=set()
  for m in mods:
   if not isinstance(m,dict) or set(m)!={'position','modification','parent_base'}: raise ValueError(f'Row {i}: each modification requires position, modification, parent_base')
   pos=m['position']; parent=m['parent_base']; mod=m['modification']
   if type(pos) is not int or not 1<=pos<=31 or pos in used or parent not in RNA or parent!=seq[pos-1] or not isinstance(mod,str) or not mod.strip(): raise ValueError(f'Row {i}: invalid modification at {pos}')
   if mod=='m6A' and parent!='A': raise ValueError(f'Row {i}: m6A requires A')
   used.add(pos)
  if (seq[16:]!='A'*15 or any(p>=17 for p in used)) and '[allow_A15_changes]' not in r['notes']: raise ValueError(f'Row {i}: A15 changes require [allow_A15_changes] in notes')
  r['rna_sequence']=seq; r['modification_status']='requires_parameterization' if mods else 'unmodified'; r['parsed_modifications']=mods
 if expected_counts is not None and Counter(r['group'] for r in rows)!=Counter(expected_counts): raise ValueError(f'Unexpected group counts: {dict(Counter(r["group"] for r in rows))}; expected {expected_counts}')
 return rows

def audit_pdb(path,uniprot_sequence=None):
 path=Path(path); lines=path.read_text(errors='replace').splitlines(); cryst=next((x for x in lines if x.startswith('CRYST1')),None)
 box=None
 if cryst:
  try: box=[float(cryst[6+i:15+i]) for i in (0,9,18)]
  except ValueError: box=None
 placeholder=bool(box and min(box)<=2.0)
 models=sum(x.startswith('MODEL ') for x in lines)
 atomlines=[x for x in lines if x.startswith(('ATOM  ','HETATM'))]
 keys=[(x[21:22],x[22:27],x[12:16],x[16:17],x[17:20]) for x in atomlines]
 duplicates=[str(k) for k,n in Counter(keys).items() if n>1]
 altloc=sum(bool(x[16:17].strip()) for x in atomlines)
 parser=PDBParser(QUIET=True,PERMISSIVE=False)
 structure=parser.get_structure('input',str(path)); model=next(structure.get_models())
 chains={c.id:list(c.get_residues()) for c in model}; protein=chains.get('A',[]); rna=chains.get('E',[])
 def one(r): return seq1(r.resname,custom_map={'UNK':'X'},undef_code='X')
 protseq=''.join(one(r) for r in protein); rnaseq=''.join(r.resname.strip() if r.resname.strip() in RNA else '?' for r in rna)
 ids=[r.id[1] for r in rna]
 missing=[]
 for j,r in enumerate(rna):
  base=r.resname.strip()
  if base not in BASE: missing.append({'residue':r.id[1],'error':'unsupported residue '+base}); continue
  need=SUGAR|BASE[base]
  if j>0: need=need|{'P','OP1','OP2'}
  m=sorted(need-{a.name for a in r})
  if m: missing.append({'residue':r.id[1],'missing_heavy_atoms':m})
 for r in protein:
  m=sorted({'N','CA','C','O'}-{a.name for a in r})
  if m: missing.append({'chain':'A','residue':r.id[1],'missing_backbone_atoms':m})
 def heavy(rs): return [a for r in rs for a in r if a.element not in ('H','D') and a.element]
 pa,ra=heavy(protein),heavy(rna); contacts={'pairs_lt_2A':None,'minimum_A':None,'closest_pair':None,'criterion':'all protein–RNA heavy-atom pairs < 2.0 Å; no interchain covalent bonds assumed'}
 if pa and ra:
  px=np.array([a.coord for a in pa]); rx=np.array([a.coord for a in ra]); tree=cKDTree(px); close=tree.query_ball_point(rx,2.0)
  contacts['pairs_lt_2A']=sum(float(np.linalg.norm(rx[j]-px[i]))<2.0 for j,ix in enumerate(close) for i in ix)
  dist,idx=tree.query(rx); j=int(np.argmin(dist)); a,b=pa[int(idx[j])],ra[j]
  def atomid(x): return {'chain':x.parent.parent.id,'residue':x.parent.id[1],'resname':x.parent.resname,'atom':x.name}
  contacts['minimum_A']=float(dist[j]);contacts['closest_pair']=[atomid(a),atomid(b)]
 bonds=[]
 for left,right in zip(rna,rna[1:]):
  d=float(np.linalg.norm(left["O3'"].coord-right['P'].coord)) if "O3'" in left and 'P' in right else None
  bonds.append({'from_residue':left.id[1],'to_residue':right.id[1],'O3_P_A':d,'within_1_3_to_2_0_A':d is not None and 1.3<=d<=2.0})
 peptide=[]
 for left,right in zip(protein,protein[1:]):
  d=float(np.linalg.norm(left['C'].coord-right['N'].coord)) if 'C' in left and 'N' in right else None
  if d is None or not 1.1<=d<=1.8: peptide.append({'from':left.id[1],'to':right.id[1],'C_N_A':d})
 errors=[]
 if models>1: errors.append('Multiple MODEL records')
 if duplicates or altloc: errors.append('Duplicate atoms or alternate locations require explicit review')
 if len(protein)!=360 or [r.id[1] for r in protein]!=list(range(1,361)): errors.append('Chain A is not the expected residues 1–360')
 if ids!=list(range(37,68)) or rnaseq!='A'*31: errors.append('Chain E is not expected 31 adenines E37–E67')
 if uniprot_sequence is not None and protseq!=uniprot_sequence[10:370]: errors.append('Chain A sequence does not exactly match P11940 residues 11–370')
 if missing: errors.append('Missing expected heavy atoms or unsupported residues')
 if any(not b['within_1_3_to_2_0_A'] for b in bonds): errors.append('RNA linkage outlier or missing bond atoms')
 if peptide: errors.append('Protein peptide connectivity outlier')
 if contacts['minimum_A'] is not None and contacts['minimum_A']<2: errors.append('Protein–RNA contacts < 2 Å require structural review')
 return {'input_file':path.name,'sha256':sha256(path),'model_records':models,'chain_residue_counts':{k:len(v) for k,v in chains.items()},'protein_sequence':protseq,'uniprot_mapping_11_370_match':None if uniprot_sequence is None else protseq==uniprot_sequence[10:370], 'rna_sequence':rnaseq,'rna_mapping':[{'position':j+1,'chain':'E','residue':r.id[1],'base':r.resname} for j,r in enumerate(rna)],'cryst1':{'line':cryst,'lengths_A':box,'placeholder_detected':placeholder},'duplicate_atom_records':duplicates,'alternate_location_atom_records':altloc,'missing_atoms':missing,'rna_linkages':bonds,'peptide_linkage_outliers':peptide,'contacts':contacts,'preflight_errors':errors,'limitations':['Side-chain completeness, absolute stereochemistry and force-field matching are not fully validated by this audit.','A 5-prime hydroxyl and 3-prime hydroxyl are intended but terminal chemistry needs explicit atom/bond and template verification.','No relaxation or repair was performed.']}

def normalized_length_analysis(utr50,utr15):
 if not isinstance(utr50,str) or not isinstance(utr15,str): raise ValueError('missing sequence')
 s=utr50.upper().strip(); short=utr15.upper().strip()
 if len(s)!=50 or len(short)!=15 or s[-15:]!=short or set(s)-RNA: raise ValueError('UTR50 does not agree with UTR15')
 run=len(s)-len(s.rstrip('A')); end=50-run
 if end<16: raise ValueError('fewer than 16 observed bases through the last non-A; cannot derive construct')
 return {'native_terminal_A_count':run,'normalized_UTR16':s[end-16:end],'candidate_sequence':s[end-16:]+'A'*15,'candidate_length':31+run,'fixed_31_supported':run==0}
