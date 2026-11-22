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

def main(b,c,out):
 rows=read(b,'Selected panel','RNA_ID','UTR50','UTR15','benchmark')+read(c,'Unmodified controls','RNA ID','50-nt native UTR','15-nt native UTR','control')
 if Counter(r['group'] for r in rows)!=Counter({'benchmark':34,'control':27}): raise ValueError('Workbook counts differ from expected 34 and 27')
 with open(out,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(json.dumps({'counts':dict(Counter(r['group'] for r in rows)),'length_counts':dict(sorted(Counter(r['candidate_length'] for r in rows).items())),'longer_than_31':sum(not r['fixed_31_supported'] for r in rows),'report':str(out)},indent=2))
if __name__=='__main__': main(*sys.argv[1:])
