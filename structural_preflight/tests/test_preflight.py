import unittest, tempfile, csv, json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from preflight import validate_variants, audit_pdb, K_RESTRAINT, normalized_length_analysis

class Validate(unittest.TestCase):
 def csvfile(self, rows):
  d=tempfile.TemporaryDirectory(); self.addCleanup(d.cleanup); p=Path(d.name)/'v.csv'
  with p.open('w',newline='') as f:
   w=csv.writer(f); w.writerow(['variant_id','group','rna_sequence','modifications_json','replicate_id','notes']); w.writerows(rows)
  return p
 def row(self,**kw):
  d={'variant_id':'b1','group':'benchmark','rna_sequence':'A'*31,'modifications_json':'[]','replicate_id':'1','notes':''}; d.update(kw); return list(d.values())
 def test_valid(self): self.assertEqual(len(validate_variants(self.csvfile([self.row()]),expected_counts=None)),1)
 def test_duplicate_id(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(),self.row()]),None)
 def test_bad_base(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(rna_sequence='A'*30+'T')]),None)
 def test_short(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(rna_sequence='A'*30)]),None)
 def test_tail(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(rna_sequence='A'*30+'G')]),None)
 def test_tail_authorized(self): self.assertEqual(len(validate_variants(self.csvfile([self.row(rna_sequence='A'*30+'G',notes='[allow_A15_changes]')]),None)),1)
 def test_bad_parent(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(modifications_json='[{"position": 2, "modification": "m6A", "parent_base": "C"}]')]),None)
 def test_duplicate_sites(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row(modifications_json='[{"position": 2,"modification":"m6A","parent_base":"A"},{"position": 2,"modification":"2OMe","parent_base":"A"}]')]),None)
 def test_modification_is_blocked(self):
  rows=validate_variants(self.csvfile([self.row(modifications_json='[{"position": 2,"modification":"m6A","parent_base":"A"}]')]),None)
  self.assertEqual(rows[0]['modification_status'],'requires_parameterization')
 def test_native_terminal_as_preserved(self):
  a=normalized_length_analysis('G'*47+'AAA','G'*12+'AAA')
  self.assertEqual(a['candidate_length'],34)
  self.assertEqual(a['candidate_sequence'],'G'*16+'A'*18)
  self.assertFalse(a['fixed_31_supported'])
 def test_utr_mismatch(self):
  with self.assertRaises(ValueError): normalized_length_analysis('G'*50,'C'*15)
 def test_units(self): self.assertEqual(K_RESTRAINT,4184.0)
 def test_invalid_box_and_contact(self):
  text=('CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n'
  'ATOM      1   CA ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
  "ATOM      2  C1'   A E  37       1.000   0.000   0.000  1.00  0.00           C\nEND\n")
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'tiny.pdb'; p.write_text(text); a=audit_pdb(p,None)
   self.assertTrue(a['cryst1']['placeholder_detected']); self.assertEqual(a['contacts']['pairs_lt_2A'],1)
 def test_counts(self):
  with self.assertRaises(ValueError): validate_variants(self.csvfile([self.row()]))
if __name__=='__main__': unittest.main()
