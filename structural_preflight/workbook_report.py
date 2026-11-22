"""Analyze actual workbooks; planning report only, no structures or accepted variants CSV."""
import sys,csv,json
from collections import Counter
from pathlib import Path
from openpyxl import load_workbook
from preflight import normalized_length_analysis

def read(path,sheet,idcol,longcol,shortcol,group):
 ws=load_workbook(path,read_only=True,data_only=True)[sheet]; it=ws.values; headers=next(it)
 rows=[]
 for values in it:
  r=dict(zip(headers,values))
  if not r[idcol]: continue
  info=normalized_length_analysis(r[longcol],r[shortcol])
  rows.append({'group':group,'source_id':r[idcol],**info,'status':'planning_only_template_extension_required' if not info['fixed_31_supported'] else 'planning_only_not_structurally_validated'})
 return rows

def write_candidates(rows,path):
 headers=['variant_id','group','rna_sequence','modifications_json','replicate_id','notes']
 seen=set()
 with open(path,'w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=headers);writer.writeheader()
  for r in rows:
   identifier=('BM_' if r['group']=='benchmark' else 'CTRL_')+str(r['source_id'])
   if identifier in seen:raise ValueError('Duplicate candidate identifier: '+identifier)
   seen.add(identifier);n=int(r['native_terminal_A_count']);seq=r['candidate_sequence']
   if len(seq)!=31+n or seq[-(15+n):]!='A'*(15+n):raise ValueError('Inconsistent candidate sequence: '+identifier)
   condition='template extension required' if n else '31-nt construct; structural preparation still required'
   writer.writerow({'variant_id':identifier,'group':r['group'],'rna_sequence':seq,'modifications_json':'[]','replicate_id':'1','notes':f"Planning only; source {r['source_id']}; 16 nt to last non-A + {n} native A + A15; {condition}"})

def main(b,c,out,candidate_path=None):
 rows=read(b,'Selected panel','RNA_ID','UTR50','UTR15','benchmark')+read(c,'Unmodified controls','RNA ID','50-nt native UTR','15-nt native UTR','control')
 if Counter(r['group'] for r in rows)!=Counter({'benchmark':34,'control':27}): raise ValueError('Workbook counts differ from expected 34 and 27')
 with open(out,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 if candidate_path:write_candidates(rows,candidate_path)
 print(json.dumps({'counts':dict(Counter(r['group'] for r in rows)),'length_counts':dict(sorted(Counter(r['candidate_length'] for r in rows).items())),'longer_than_31':sum(not r['fixed_31_supported'] for r in rows),'report':str(out)},indent=2))
if __name__=='__main__': main(*sys.argv[1:])
