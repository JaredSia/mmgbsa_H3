import sys,unittest,tempfile,json
from unittest.mock import patch
from pathlib import Path
import openmm as mm
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import reference_runner as rr

class CheckpointTests(unittest.TestCase):
 def fixture(self,d):
  s=mm.System();s.addParticle(12);inte=mm.VerletIntegrator(.001);ctx=mm.Context(s,inte,mm.Platform.getPlatformByName('Reference'));ctx.setPositions([[.1,0,0]])
  state=Path(d)/'s_state.xml';state.write_text(mm.XmlSerializer.serialize(ctx.getState(getPositions=True)))
  pdb=Path(d)/'s.pdb';pdb.write_text('example atom-file payload')
  qa={'stereochemistry_outliers':[],'RNA_bad_link_count':0,'interchain_pairs_lt_2A':0,'interchain_minimum_A':2.5,'protein_CA_RMSD_protein_aligned_A':.3,'RNA_backbone_RMSD_protein_aligned_A':.2,'protein_CA_max_displacement_unaligned_A':.7,'RNA_backbone_max_displacement_unaligned_A':.8,'stereochemistry_penalty_kJ_mol':0.0,'reporter_threshold_observed':True,'artifact_sha256':{'state':rr.digest(state),'pdb':rr.digest(pdb)},'upstream_sha256':'expected'}
  (Path(d)/'s_diagnostics.json').write_text(json.dumps(qa));return ctx,inte,qa
 def test_verified_checkpoint_reused(self):
  with tempfile.TemporaryDirectory() as d:
   ctx,inte,qa=self.fixture(d)
   with patch.object(rr,'minimize_stage',side_effect=AssertionError('unexpected rerun')):
    result,pos=rr.advance_stage('s',ctx,None,None,None,Path(d),10,0,0,upstream_sha='expected')
   self.assertEqual(result['upstream_sha256'],'expected')
 def test_changed_upstream_invalidates_checkpoint(self):
  with tempfile.TemporaryDirectory() as d:
   ctx,inte,qa=self.fixture(d)
   with patch.object(rr,'minimize_stage',return_value=('rerun',None)) as mocked:
    result,_=rr.advance_stage('s',ctx,None,None,None,Path(d),10,0,0,upstream_sha='different')
   self.assertEqual(result,'rerun');self.assertTrue(mocked.called)
 def test_changed_state_invalidates_checkpoint(self):
  with tempfile.TemporaryDirectory() as d:
   ctx,inte,qa=self.fixture(d);(Path(d)/'s.pdb').write_text('tampered')
   with patch.object(rr,'minimize_stage',return_value=('rerun',None)) as mocked:
    result,_=rr.advance_stage('s',ctx,None,None,None,Path(d),10,0,0,upstream_sha='expected')
   self.assertEqual(result,'rerun');self.assertTrue(mocked.called)
if __name__=='__main__':unittest.main()
