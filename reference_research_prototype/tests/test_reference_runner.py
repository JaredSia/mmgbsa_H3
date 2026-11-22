import sys,unittest
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reference_runner import Tracker,stereo_barrier,validate_geometry,convergence_observed,minimize_with_retries

class ReferenceRunnerTests(unittest.TestCase):
 def test_stereo_barrier_is_zero_on_correct_side_and_positive_on_inversion(self):
  system=mm.System()
  for i in range(4):system.addParticle(12)
  system.addForce(stereo_barrier({'stereo':[{'indices':[0,1,2,3],'expected_sign':1}]}))
  ctx=mm.Context(system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'))
  ctx.setPositions([[0,0,0],[.14,0,0],[0,.14,0],[0,0,.14]])
  good=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
  ctx.setPositions([[0,0,0],[.14,0,0],[0,.14,0],[0,0,-.14]])
  bad=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
  self.assertEqual(good,0);self.assertGreater(bad,1e5)
 def test_openmm_reporter_records_real_iterations(self):
  system=mm.System();system.addParticle(1);f=mm.CustomExternalForce('x*x');f.addParticle(0,[]);system.addForce(f)
  ctx=mm.Context(system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'));ctx.setPositions([[1,0,0]])
  t=Tracker();mm.LocalEnergyMinimizer.minimize(ctx,1e-9,5,t)
  self.assertGreater(t.calls,0);self.assertIn('objective_gradient_rms_kJ_mol_nm',t.last)
 def test_minimizer_retry_wrapper_records_attempts(self):
  system=mm.System();system.addParticle(1);force=mm.CustomExternalForce('x*x');force.addParticle(0,[]);system.addForce(force)
  ctx=mm.Context(system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'));ctx.setPositions([[1,0,0]])
  tracker,attempts=minimize_with_retries(ctx,1e-9,10)
  self.assertGreaterEqual(len(attempts),1);self.assertGreater(tracker.calls,0)
 def test_reporter_threshold_requires_both_gradient_and_constraint(self):
  self.assertFalse(convergence_observed(None,10))
  self.assertTrue(convergence_observed({'objective_gradient_rms_kJ_mol_nm':9,'max_constraint_error':1e-6},10))
  self.assertFalse(convergence_observed({'objective_gradient_rms_kJ_mol_nm':11,'max_constraint_error':1e-6},10))
  self.assertFalse(convergence_observed({'objective_gradient_rms_kJ_mol_nm':9,'max_constraint_error':1e-4},10))
 def test_quality_guard_rejects_bad_structures(self):
  d={'stereochemistry_outliers':[],'RNA_bad_link_count':0,'interchain_pairs_lt_2A':0,'interchain_minimum_A':2.5,'protein_CA_RMSD_protein_aligned_A':.3,'RNA_backbone_RMSD_protein_aligned_A':.2,'protein_CA_max_displacement_unaligned_A':.7,'RNA_backbone_max_displacement_unaligned_A':.8}
  validate_geometry(d,True)
  for key,value in [('stereochemistry_outliers',[{'id':'test'}]),('RNA_bad_link_count',1),('interchain_pairs_lt_2A',1),('interchain_minimum_A',float('nan')),('protein_CA_RMSD_protein_aligned_A',1.1),('RNA_backbone_max_displacement_unaligned_A',3.1)]:
   bad=d.copy();bad[key]=value
   with self.assertRaises(RuntimeError):validate_geometry(bad,True)
if __name__=='__main__':unittest.main()
