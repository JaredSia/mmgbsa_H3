import unittest, tempfile, csv, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from workbook_report import write_candidates
class CandidateTests(unittest.TestCase):
 def test_full_schema_and_unique_group_ids(self):
  rows=[{'group':'benchmark','source_id':'RNA_1','candidate_sequence':'G'*16+'A'*15,'native_terminal_A_count':0,'candidate_length':31},
        {'group':'control','source_id':'RNA_1','candidate_sequence':'G'*16+'A'*18,'native_terminal_A_count':3,'candidate_length':34}]
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'candidates.csv';write_candidates(rows,p)
   with p.open(newline='') as f:
    r=csv.DictReader(f);self.assertEqual(r.fieldnames,['variant_id','group','rna_sequence','modifications_json','replicate_id','notes']);data=list(r)
   self.assertNotEqual(data[0]['variant_id'],data[1]['variant_id']);self.assertEqual(data[1]['rna_sequence'],'G'*16+'A'*18);self.assertIn('extension required',data[1]['notes'])
if __name__=='__main__':unittest.main()
