"""Experimental reference-only restrained preparation. No affinity prediction.

Generates *candidate* structures and detailed quality reports. Geometry and
convergence checks are necessary, not sufficient, for biological validation.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, random, shutil, time, urllib.request
import numpy as np
from scipy.spatial import cKDTree
import openmm as mm
from openmm import app, unit
import preflight
from preflight import audit_pdb

RNA_BB={"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"}
FF_NAMES=['amber19/protein.ff19SB.xml','amber14/RNA.OL3.xml','amber19/tip3p.xml']

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic_json(p,obj):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2,allow_nan=False));t.replace(p)
def context_for(system,integrator,platform):
 order=['CUDA','CPU'] if platform=='auto' else [platform]
 errors=[]
 for name in order:
  try:
   p=mm.Platform.getPlatformByName(name);opts={'Threads':os.environ.get('OPENMM_CPU_THREADS','2')} if name=='CPU' else {}
   return mm.Context(system,integrator,p,opts),name,errors
  except Exception as e:errors.append(f'{name}: {type(e).__name__}: {e}')
 raise RuntimeError('No requested OpenMM platform available: '+'; '.join(errors))

def selections(topology):
 atoms=list(topology.atoms()); protein=[];rna=[];ca=[];protein_bb=[];bb=[];stereo=[];links=[]
 residues=list(topology.residues())
 for a in atoms:
  if a.element is None:raise ValueError(f'Atom has no element: {a}')
  if a.residue.chain.id=='A':
   if a.element.symbol not in ('H','D'):protein.append(a.index)
   if a.name=='CA':ca.append(a.index)
   if a.name in ('N','CA','C','O'):protein_bb.append(a.index)
  if a.residue.chain.id=='E':
   if a.element.symbol not in ('H','D'):rna.append(a.index)
   if a.name in RNA_BB:bb.append(a.index)
 rna_res=[r for r in residues if r.chain.id=='E'];protein_res=[r for r in residues if r.chain.id=='A']
 if len(rna_res)!=31 or len(protein_res)!=360 or len(ca)!=360:raise ValueError('Reference chain/residue counts differ from E31/A360')
 for r in rna_res+protein_res:
  ns={a.name:a.index for a in r.atoms()};spec=[]
  if r.chain.id=='E':
   gly='N9' if r.name in ('A','G') else 'N1'
   spec=[("C1'",("O4'","C2'",gly),1),("C2'",("C1'","C3'","O2'"),-1),("C3'",("C2'","C4'","O3'"),-1),("C4'",("C3'","O4'","C5'"),1)]
  elif r.name!='GLY':spec=[('CA',('N','C','CB'),1)]
  if r.chain.id=='A' and r.name=='ILE':spec.append(('CB',('CA','CG1','CG2'),1))
  if r.chain.id=='A' and r.name=='THR':spec.append(('CB',('CA','OG1','CG2'),1))
  for center,group,sgn in spec:
   if not all(k in ns for k in (center,*group)):raise ValueError(f'Cannot check chirality {r} {center}')
   stereo.append({'residue':str(r),'center':center,'indices':[ns[k] for k in (center,*group)],'expected_sign':sgn})
 for l,r in zip(rna_res,rna_res[1:]):
  x={a.name:a.index for a in l.atoms()};y={a.name:a.index for a in r.atoms()};links.append((x["O3'"],y['P']))
 return {'protein':protein,'rna':rna,'protein_ca':ca,'protein_backbone':protein_bb,'rna_backbone':bb,'stereo':stereo,'links':links}

def rms(a):return float(np.sqrt(np.mean(np.sum(a*a,axis=1))))
def diagnostics(x,initial,sel):
 # All positions in nm, converted to Å for geometry.
 xyz=x*10; ref=initial*10;ca=sel['protein_ca'];bb=sel['rna_backbone']
 ref_c=ref[ca].mean(0);mob_c=xyz[ca].mean(0);u,s,vt=np.linalg.svd((xyz[ca]-mob_c).T@(ref[ca]-ref_c));rot=u@np.diag([1,1,np.linalg.det(u@vt)])@vt; aligned=(xyz-mob_c)@rot+ref_c
 rna=sel['rna'];protein=sel['protein'];tree=cKDTree(xyz[protein]);near=tree.query_ball_point(xyz[rna],2.0)
 pair_count=int(sum(np.linalg.norm(xyz[rna[j]]-xyz[protein[i]])<2.0 for j,ns in enumerate(near) for i in ns));minimum=float(tree.query(xyz[rna])[0].min())
 bad=[];volume_min=float('inf')
 for row in sel['stereo']:
  center,*nb=row['indices'];v=float(np.linalg.det(np.array([xyz[i]-xyz[center] for i in nb])));signed=v*row['expected_sign'];volume_min=min(volume_min,signed)
  if signed<0.5:bad.append({'residue':row['residue'],'center':row['center'],'signed_volume_A3':v})
 distances=[float(np.linalg.norm(xyz[i]-xyz[j])) for i,j in sel['links']]
 return {'protein_CA_RMSD_protein_aligned_A':rms(aligned[ca]-ref[ca]),'protein_backbone_RMSD_protein_aligned_A':rms(aligned[sel['protein_backbone']]-ref[sel['protein_backbone']]),'RNA_backbone_RMSD_protein_aligned_A':rms(aligned[bb]-ref[bb]),'protein_CA_max_displacement_unaligned_A':float(np.linalg.norm(xyz[ca]-ref[ca],axis=1).max()),'RNA_backbone_max_displacement_unaligned_A':float(np.linalg.norm(xyz[bb]-ref[bb],axis=1).max()),'interchain_minimum_A':minimum,'interchain_pairs_lt_2A':pair_count,'RNA_O3_P_range_A':[min(distances),max(distances)],'RNA_bad_link_count':sum(not 1.3<=x<=2.0 for x in distances),'checked_chiral_centers':len(sel['stereo']),'smallest_expected_signed_volume_A3':volume_min,'stereochemistry_outliers':bad}

def validate_geometry(d,require_no_interchain_clash):
 issues=[]
 if d['stereochemistry_outliers']:issues.append('inverted or nearly planar stereocenters')
 if d['RNA_bad_link_count']:issues.append('RNA phosphodiester distance outliers')
 if d['protein_CA_RMSD_protein_aligned_A']>1.0 or d['RNA_backbone_RMSD_protein_aligned_A']>1.0:issues.append('protein or RNA backbone RMSD exceeds the conservative 1 Å review threshold')
 if d['protein_CA_max_displacement_unaligned_A']>3.0 or d['RNA_backbone_max_displacement_unaligned_A']>3.0:issues.append('restrained backbone atom displacement exceeds the conservative 3 Å review threshold')
 if require_no_interchain_clash and d['interchain_pairs_lt_2A']:issues.append('interchain heavy-atom pairs remain below 2 Å')
 if not math.isfinite(d['interchain_minimum_A']):issues.append('nonfinite geometry')
 if issues:raise RuntimeError('; '.join(issues))

def convergence_observed(last,tolerance):
 return bool(last is not None and last['objective_gradient_rms_kJ_mol_nm']<=tolerance and last['max_constraint_error']<=1e-5)

class Tracker(mm.MinimizationReporter):
 def __init__(self):super().__init__();self.calls=0;self.last=None
 def report(self,iteration,x,grad,args):
  self.calls+=1;g=np.asarray(grad,dtype=float)
  self.last={'iteration':int(iteration),'objective_gradient_rms_kJ_mol_nm':float(np.sqrt(np.mean(g*g))),'max_constraint_error':float(args['max constraint error'])}
  return False

def stereo_barrier(sel):
 # Numerical anti-inversion aid. It is not part of ff19SB/OL3 and must be
 # reported separately; zero final energy does not make the pathway unbiased.
 expr=('0.5*kc*step(vmin-s*v)*(vmin-s*v)^2;'
       'v=(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))'
       '-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))'
       '+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))')
 force=mm.CustomCompoundBondForce(4,expr);force.addGlobalParameter('kc',1e12);force.addGlobalParameter('vmin',.001);force.addPerBondParameter('s');force.setForceGroup(3)
 for row in sel['stereo']:force.addBond(row['indices'],[row['expected_sign']])
 return force

def restraint_force(name,indices,reference,group):
 f=mm.CustomExternalForce(f'0.5*{name}*((x-x0)^2+(y-y0)^2+(z-z0)^2)');f.addGlobalParameter(name,0.0)
 for p in ('x0','y0','z0'):f.addPerParticleParameter(p)
 for i in indices:f.addParticle(i,reference[i])
 f.setForceGroup(group);return f

def save_pdb(path,top,pos):
 with open(path,'w') as fh:app.PDBFile.writeFile(top,pos,fh,keepIds=True)
 loaded=app.PDBFile(str(path))
 if loaded.topology.getNumAtoms()!=top.getNumAtoms():raise RuntimeError(f'Reloaded atom count mismatch in {path}')

def minimize_with_retries(ctx,tolerance,iterations):
 # Resetting L-BFGS can help when it exits early. It does not prove convergence.
 remaining=iterations;attempts=[];tracker=None
 for _ in range(3):
  if remaining<=0:break
  tracker=Tracker();budget=remaining;mm.LocalEnergyMinimizer.minimize(ctx,tolerance,budget,tracker)
  attempts.append({'iteration_budget':budget,'reporter_calls':tracker.calls,'last_reporter':tracker.last})
  if convergence_observed(tracker.last,tolerance) or tracker.calls==0 or tracker.calls>=budget:break
  remaining-=tracker.calls
 return tracker,attempts

def minimize_stage(name,ctx,top,sel,original,out,iterations,all_k,anchor_k,tolerance=10.0,soft_stage=False,allow_unconverged=False,upstream_sha=None):
 ctx.setParameter('k_all',all_k);ctx.setParameter('k_anchor',anchor_k); t=time.monotonic()
 tracker,attempts=minimize_with_retries(ctx,tolerance,iterations)
 state=ctx.getState(getPositions=True,getEnergy=True,getForces=True,enforcePeriodicBox=False)
 pos=state.getPositions(asNumpy=True);xyz=np.asarray(pos.value_in_unit(unit.nanometer));force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))
 if not np.isfinite(xyz).all() or not np.isfinite(force).all():raise RuntimeError('Nonfinite coordinates or forces in '+name)
 def energy(groups):return float(ctx.getState(getEnergy=True,groups=groups).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
 qa=diagnostics(xyz,original,sel); qa.update({'stage':name,'energy_interpretation':'artificial soft repulsion and bonded terms; not full force field' if soft_stage else 'full force field','seconds':time.monotonic()-t,'max_iterations':iterations,'base_system_energy_kJ_mol':energy({0}),'physical_forcefield_energy_kJ_mol':None if soft_stage else energy({0}),'positional_restraint_energy_kJ_mol':energy({1,2}),'stereochemistry_penalty_kJ_mol':energy({3}),'total_energy_kJ_mol':energy({0,1,2,3}),'maximum_unprojected_force_kJ_mol_nm':float(np.linalg.norm(force,axis=1).max()),'reporter_calls':sum(a['reporter_calls'] for a in attempts),'minimizer_attempts':attempts,'last_reporter':tracker.last,'convergence_status':'not established; reporter gradient includes constraint penalties and early termination is possible','restraint_definition':'E=0.5*k*sum(displacement^2)','k_all_kJ_mol_nm2':all_k,'k_anchor_kJ_mol_nm2':anchor_k})
 observed=convergence_observed(tracker.last,tolerance)
 qa['reporter_threshold_observed']=observed
 # A diagnostic is the commit marker: it is written only after the PDB and
 # exact OpenMM state are saved, reloaded and hashed.
 pdbpath=out/(name+'.pdb');statepath=out/(name+'_state.xml')
 save_pdb(pdbpath,top,pos)
 tmppath=out/(name+'_state.xml.tmp');tmppath.write_text(mm.XmlSerializer.serialize(state));tmppath.replace(statepath)
 qa['artifact_sha256']={'pdb':digest(pdbpath),'state':digest(statepath)};qa['upstream_sha256']=upstream_sha
 atomic_json(out/(name+'_diagnostics.json'),qa)
 validate_geometry(qa,not soft_stage)
 if not soft_stage and abs(qa['stereochemistry_penalty_kJ_mol'])>1e-6:raise RuntimeError('Stereochemistry penalty remains active; do not use this stage as an unbiased force-field minimum')
 if not soft_stage and not observed and not allow_unconverged:raise RuntimeError(f'{name} did not show the requested gradient/constraint thresholds; inspect diagnostics; early optimizer termination or poor geometry may persist even with more iterations')
 print(name, 'iterations reported',tracker.calls,'min distance Å',round(qa['interchain_minimum_A'],3),'RNA bb RMSD Å',round(qa['RNA_backbone_RMSD_protein_aligned_A'],3),'seconds',round(qa['seconds'],1),flush=True)
 return qa,pos

def advance_stage(name,ctx,top,sel,original,out,iterations,all_k,anchor_k,soft_stage=False,allow_unconverged=False,upstream_sha=None):
 diag=out/(name+'_diagnostics.json');state=out/(name+'_state.xml');pdb=out/(name+'.pdb')
 if diag.exists() and state.exists() and pdb.exists():
  try:
   old=json.loads(diag.read_text());h=old['artifact_sha256']
   intact=h['state']==digest(state) and h['pdb']==digest(pdb) and old.get('upstream_sha256')==upstream_sha
  except (ValueError,KeyError,OSError):intact=False
  if intact:
   validate_geometry(old,not soft_stage)
   if not soft_stage and abs(old['stereochemistry_penalty_kJ_mol'])>1e-6:raise RuntimeError(f'{name}: saved stereochemistry penalty is active; review before continuing')
   ctx.setState(mm.XmlSerializer.deserialize(state.read_text()))
   if soft_stage or old['reporter_threshold_observed']:
    print('Reused completed checkpoint:',name,flush=True)
    return old,ctx.getState(getPositions=True).getPositions(asNumpy=True)
   history=out/(name+'_history.jsonl')
   with open(history,'a') as fh:fh.write(json.dumps(old,allow_nan=False)+'\n')
   print('Continuing unconverged stage from saved state:',name,flush=True)
  else:print('Ignoring incomplete checkpoint and restarting this stage:',name,flush=True)
 return minimize_stage(name,ctx,top,sel,original,out,iterations,all_k,anchor_k,soft_stage=soft_stage,allow_unconverged=allow_unconverged,upstream_sha=upstream_sha)

def main(source,out,platform='auto',dry_steps=2000,water_steps=2000,stage_steps=2000,seed=9181,allow_unconverged=False):
 if min(dry_steps,water_steps,stage_steps)<=0:raise ValueError('Iteration budgets must be positive; unlimited minimization is disabled')
 source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 identity={'source_sha256':digest(source),'runner_sha256':digest(__file__),'preflight_sha256':digest(preflight.__file__),'openmm_version':mm.__version__,'seed':seed,'forcefields':FF_NAMES,'experimental_allow_unconverged':allow_unconverged}
 identity_path=out/'run_identity.json';initial_pdb=out/'initial_protonated_solute.pdb';initial_npy=out/'initial_positions_nm.npy'
 ff=app.ForceField(*FF_NAMES)
 if identity_path.exists():
  old=json.loads(identity_path.read_text())
  if old['identity']!=identity:raise ValueError('Cannot resume: input, code, seed or OpenMM version differs from the original run')
  if old['initial_pdb_sha256']!=digest(initial_pdb) or old['initial_positions_sha256']!=digest(initial_npy):raise ValueError('Cannot resume: initial checkpoint checksum mismatch')
  loaded=app.PDBFile(str(initial_pdb));original=np.load(initial_npy,allow_pickle=False)
  if original.shape!=(loaded.topology.getNumAtoms(),3):raise ValueError('Cannot resume: initial atom count mismatch')
  mod=app.Modeller(loaded.topology,original*unit.nanometer);print('Resuming run from verified checkpoints:',out,flush=True)
 else:
  if any(out.iterdir()):raise FileExistsError('Output folder contains files but no valid run identity; choose a new folder')
  with urllib.request.urlopen('https://rest.uniprot.org/uniprotkb/P11940.fasta',timeout=30) as r:fasta=r.read().decode()
  seq=''.join(fasta.splitlines()[1:]);a=audit_pdb(source,seq);atomic_json(out/'initial_audit.json',a)
  allowed_errors=('5-prime terminal phosphate atoms','CRYST1 placeholder box','Protein–RNA contacts < 2')
  other_errors=[error for error in a['preflight_errors'] if not error.startswith(allowed_errors)]
  if other_errors:raise ValueError('Input identity or geometry failures: '+str(other_errors))
  if a['uniprot_mapping_11_370_match'] is not True or a['rna_sequence']!='A'*31 or len(a['rna_mapping'])!=31 or [r['residue'] for r in a['rna_mapping']]!=list(range(37,68)) or a['missing_atoms'] or a['duplicate_atom_records'] or a['alternate_location_atom_records'] or a['peptide_linkage_outliers'] or any(not z['within_1_3_to_2_0_A'] for z in a['rna_linkages']) or a['model_records']>1:
   raise ValueError('Input fails mandatory identity, missing atom, or backbone checks; inspect initial_audit.json')
  lines=source.read_text().splitlines(keepends=True);terminal=[s for s in lines if s.startswith('ATOM  ') and s[21:22]=='E' and s[22:26].strip()=='37' and s[12:16].strip() in {'P','OP1','OP2','OP3'}]
  names={s[12:16].strip() for s in terminal}
  if names not in (set(),{'P','OP1','OP2'}):raise ValueError('Unexpected 5-prime terminal phosphate atoms; manual chemistry review required')
  prepared=out/'input_5OH_no_box.pdb';prepared.write_text(''.join(s for s in lines if not s.startswith('CRYST1') and s not in terminal))
  atomic_json(out/'input_conversion.json',{'source_sha256':digest(source),'intended_5prime':'OH','removed_terminal_atoms':[s.strip() for s in terminal],'removed_CRYST1_records':[s.strip() for s in lines if s.startswith('CRYST1')],'note':'explicit chemical conversion per user request; original is not overwritten'})
  pdb=app.PDBFile(str(prepared));mod=app.Modeller(pdb.topology,pdb.positions);random.seed(seed)
  hp=mm.Platform.getPlatformByName('CPU');mod.addHydrogens(ff,pH=7.4,platform=hp)
  if ff.getUnmatchedResidues(mod.topology):raise ValueError('Unmatched force-field residues: '+str(ff.getUnmatchedResidues(mod.topology)))
  def heavy_coordinates(top,positions):
   pos=np.asarray(positions.value_in_unit(unit.nanometer));return {(a.residue.chain.id,a.residue.id,a.name):pos[a.index] for a in top.atoms() if a.element is not None and a.element.symbol not in ('H','D')}
  original_heavy=heavy_coordinates(pdb.topology,pdb.positions);with_hydrogen=heavy_coordinates(mod.topology,mod.positions)
  if original_heavy.keys()!=with_hydrogen.keys() or any(np.linalg.norm(original_heavy[k]-with_hydrogen[k])>1e-7 for k in original_heavy):raise ValueError('Adding hydrogens changed an existing heavy atom or its coordinates')
  rna_res=[r for r in mod.topology.residues() if r.chain.id=='E'];first={a.name for a in rna_res[0].atoms()};last={a.name for a in rna_res[-1].atoms()}
  if not {"O5'","HO5'"}<=first or 'P' in first or not {"O3'","HO3'"}<=last:raise ValueError('Requested 5-prime OH or assumed 3-prime OH could not be confirmed')
  original=np.asarray(mod.positions.value_in_unit(unit.nanometer));save_pdb(initial_pdb,mod.topology,mod.positions);np.save(initial_npy,original,allow_pickle=False)
  atomic_json(identity_path,{'identity':identity,'initial_pdb_sha256':digest(initial_pdb),'initial_positions_sha256':digest(initial_npy)})
 sel=selections(mod.topology);natoms=len(original)
 validate_geometry(diagnostics(original,original,sel),False)
 if (out/'final_status.json').exists():
  done=json.loads((out/'final_status.json').read_text());candidate=out/'reference_candidate_solute.pdb'
  if candidate.exists() and done.get('candidate_sha256')==digest(candidate):print('Verified final candidate already exists:',candidate,flush=True);return
  raise ValueError('Final status exists but candidate is missing or changed')
 platform_info={'requested':platform,'openmm_version':mm.__version__,'forcefields':FF_NAMES,'seed':seed,'salt':'0.089 M KCl; a proxy for potassium acetate, without magnesium, glycerol or other assay components','water':'TIP3P; rigidWater=True','ions':'K+, Cl-','experimental_allow_unconverged':allow_unconverged,'terminal_chemistry':'5-prime OH; 3-prime OH assumed','pH_for_hydrogen_placement':7.4,'padding_nm_per_axis':1.0,'caveats':['PABPC1 is a truncated model with artificial charged termini.','Force field TIP3P is a modeling choice; ff19SB was developed with OPC.','Soft repulsion, vacuum preconditioning and the anti-inversion penalty are artificial preparation steps; pathway bias remains even when the penalty is zero at the end.','Reference binding geometry is not experimentally validated by minimization.']}
 atomic_json(out/'configuration.json',platform_info)
 # Reduce catastrophic overlaps using an explicitly nonphysical repulsive ramp.
 # Keep canonical bonded forces; omit electrostatics and replace nonbonded LJ
 # only for this preparation stage. No such energy is a force-field energy.
 softsys=ff.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False)
 nbidx=next(i for i in range(softsys.getNumForces()) if isinstance(softsys.getForce(i),mm.NonbondedForce));nb=softsys.getForce(nbidx)
 exclusions=[tuple(int(x) for x in nb.getExceptionParameters(i)[:2]) for i in range(nb.getNumExceptions())]
 soft=mm.CustomNonbondedForce('0.5*krep*step(0.5*(diam1+diam2)-r)*(r-0.5*(diam1+diam2))^2');soft.addGlobalParameter('krep',1000.);soft.addPerParticleParameter('diam');soft.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic);soft.setCutoffDistance(.30)
 for atom in mod.topology.atoms():soft.addParticle([.17 if atom.element.symbol=='H' else .27])
 for i,j in exclusions:soft.addExclusion(i,j)
 softsys.removeForce(nbidx);softsys.addForce(soft);softsys.addForce(restraint_force('k_all',[],original,1));softsys.addForce(restraint_force('k_anchor',sel['protein_ca']+sel['rna_backbone'],original,2));softsys.addForce(stereo_barrier(sel))
 inte=mm.VerletIntegrator(.001);ctx,platform_used,errs=context_for(softsys,inte,platform);ctx.setPositions(mod.positions)
 platform_info.update({'soft_platform':platform_used,'soft_platform_fallback_errors':errs,'soft_krep_schedule':[1000.,10000.,100000.,1000000.], 'soft_steps_per_stage':300,'soft_particle_diameters_nm':{'H':.17,'other':.27},'soft_cutoff_nm':.30,'soft_exclusions':'all standard nonbonded exceptions, including 1-4 terms','backbone_RMSD_review_threshold_A':1.0,'max_restrained_atom_displacement_review_threshold_A':3.0,'chirality_barrier_k_kJ_mol_nm6':1e12,'chirality_barrier_min_volume_nm3':.001});atomic_json(out/'configuration.json',platform_info)
 for i,repulsion in enumerate((1000.,10000.,100000.,1000000.),1):
  ctx.setParameter('krep',repulsion);upstream=digest(initial_npy if i==1 else out/f'soft_precondition_{i-1}_state.xml');d,positions=advance_stage(f'soft_precondition_{i}',ctx,mod.topology,sel,original,out,300,0.,4184.,soft_stage=True,upstream_sha=upstream)
 del ctx,inte
 dry=ff.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False)
 dry.addForce(restraint_force('k_all',[],original,1));dry.addForce(restraint_force('k_anchor',sel['protein_ca']+sel['rna_backbone'],original,2));dry.addForce(stereo_barrier(sel))
 inte=mm.VerletIntegrator(.001);ctx,platform_used,errs=context_for(dry,inte,platform);ctx.setPositions(positions)
 platform_info.update({'dry_platform':platform_used,'dry_platform_fallback_errors':errs});atomic_json(out/'configuration.json',platform_info)
 d,positions=advance_stage('dry_precondition_exploratory',ctx,mod.topology,sel,original,out,dry_steps,0.,4184.,allow_unconverged=allow_unconverged,upstream_sha=digest(out/'soft_precondition_4_state.xml'));del ctx,inte
 mod.positions=positions
 # Cache the exact solvated starting coordinates and System before stage 0.
 # A marker with hashes is written last so incomplete writes cannot be reused.
 solvent_meta=out/'solvent_initial_manifest.json';solvent_pdb=out/'solvent_initial.pdb';solvent_npy=out/'solvent_initial_positions_nm.npy';system_path=out/'solvated_system.xml'
 if solvent_meta.exists():
  metadata=json.loads(solvent_meta.read_text())
  if metadata.get('dry_state_sha256')!=digest(out/'dry_precondition_exploratory_state.xml'):raise ValueError('Cannot resume: solvated checkpoint was built from a different dry state; choose a new output folder')
  for key,path in [('pdb',solvent_pdb),('positions',solvent_npy),('system',system_path)]:
   if not path.exists() or digest(path)!=metadata['sha256'][key]:raise ValueError('Cannot resume: solvated '+key+' checkpoint checksum mismatch')
  saved=app.PDBFile(str(solvent_pdb));sol_x=np.load(solvent_npy,allow_pickle=False)
  if sol_x.shape!=(saved.topology.getNumAtoms(),3):raise ValueError('Cannot resume: solvated atom count mismatch')
  def signatures(top):return [(a.residue.chain.id,a.residue.id,a.residue.name,a.name) for a in top.atoms()]
  if signatures(saved.topology)[:natoms]!=signatures(mod.topology):raise ValueError('Cannot resume: solvated solute atom mapping changed')
  mod=app.Modeller(saved.topology,sol_x*unit.nanometer);system=mm.XmlSerializer.deserialize(system_path.read_text())
  print('Reused verified solvated system:',mod.topology.getNumAtoms(),'atoms',flush=True)
 else:
  span=np.ptp(np.asarray(positions.value_in_unit(unit.nanometer)),axis=0);box=mm.Vec3(*(span+2.0))*unit.nanometer
  random.seed(seed);mod.addSolvent(ff,model='tip3p',boxSize=box,positiveIon='K+',negativeIon='Cl-',ionicStrength=0.089*unit.molar,neutralize=True)
  system=ff.createSystem(mod.topology,nonbondedMethod=app.PME,nonbondedCutoff=.9*unit.nanometer,constraints=None,rigidWater=True,removeCMMotion=False)
  sol_x=np.asarray(mod.positions.value_in_unit(unit.nanometer));system.addForce(restraint_force('k_all',list(range(natoms)),sol_x,1));system.addForce(restraint_force('k_anchor',sel['protein_ca']+sel['rna_backbone'],original,2));system.addForce(stereo_barrier(sel))
  save_pdb(solvent_pdb,mod.topology,mod.positions);np.save(solvent_npy,sol_x,allow_pickle=False);system_path.write_text(mm.XmlSerializer.serialize(system))
  atomic_json(solvent_meta,{'dry_state_sha256':digest(out/'dry_precondition_exploratory_state.xml'),'sha256':{'pdb':digest(solvent_pdb),'positions':digest(solvent_npy),'system':digest(system_path)},'total_atoms':mod.topology.getNumAtoms(),'solute_atoms':natoms})
 inte=mm.VerletIntegrator(.001);ctx,platform_used,errs=context_for(system,inte,platform);ctx.setPositions(mod.positions)
 vectors=system.getDefaultPeriodicBoxVectors();lengths=[float(v[i].value_in_unit(unit.nanometer)) for i,v in enumerate(vectors)]
 platform_info.update({'solvated_platform':platform_used,'solvated_platform_fallback_errors':errs,'solute_atoms':natoms,'total_atoms':mod.topology.getNumAtoms(),'box_lengths_nm':lengths});atomic_json(out/'configuration.json',platform_info)
 stages=[('stage0_water',water_steps,418400.,0.),('stage1_strong',stage_steps,0.,4184.),('stage2_moderate',stage_steps,0.,1255.2)]
 results=[]
 for idx,(name,iters,kall,kanchor) in enumerate(stages):
  upstream=digest(solvent_npy if idx==0 else out/(stages[idx-1][0]+'_state.xml'))
  d,positions=advance_stage(name,ctx,mod.topology,sel,original,out,iters,kall,kanchor,allow_unconverged=allow_unconverged,upstream_sha=upstream);results.append(d)
 mod.positions=positions
 # Make a separate solute-only file; the full solvated system is retained above.
 mod.delete([a for a in mod.topology.atoms() if a.index>=natoms]);save_pdb(out/'reference_candidate_solute.pdb',mod.topology,mod.positions)
 atomic_json(out/'final_status.json',{'status':'candidate_for_review_not_scientifically_validated','candidate_sha256':digest(out/'reference_candidate_solute.pdb'),'all_requested_stages_ran':True,'reporter_thresholds_observed_all_solvated':all(r['reporter_threshold_observed'] for r in results),'convergence_proven':False,'notes':['Gradient and constraint details are in stage diagnostics; iteration caps do not prove convergence.','No claim of biological pose validation or binding affinity.','Within-stage interruption replays that stage from its last saved checkpoint; platform changes may alter numerical results.','All-atom clash and independent stereochemical review remain required.'],'stages':results})
 print('Candidate completed; independent review and convergence assessment remain required.',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('out');p.add_argument('--platform',choices=['auto','CPU','CUDA','Reference'],default='auto');p.add_argument('--dry-steps',type=int,default=2000);p.add_argument('--water-steps',type=int,default=2000);p.add_argument('--stage-steps',type=int,default=2000);p.add_argument('--seed',type=int,default=9181);p.add_argument('--experimental-allow-unconverged',action='store_true');a=p.parse_args();main(a.source,a.out,a.platform,a.dry_steps,a.water_steps,a.stage_steps,a.seed,a.experimental_allow_unconverged)
