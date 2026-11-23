import sys,unittest
from pathlib import Path
import numpy as np
import openmm as mm
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import reference_runner as rr

class FinalConstraintTests(unittest.TestCase):
 def test_preprojection_report_may_exceed_target_but_final_satisfy(self):
  last={'objective_gradient_rms_kJ_mol_nm':3.55,'max_constraint_error':2.84e-5}
  self.assertTrue(rr.convergence_observed(last,10.0,3e-8,0.00003,3.0))
 def test_final_constraint_failure_is_not_accepted(self):
  last={'objective_gradient_rms_kJ_mol_nm':3.55,'max_constraint_error':2.84e-5}
  self.assertFalse(rr.convergence_observed(last,10.0,2e-5,0.00003,3.0))
 def test_large_or_unknown_projection_is_not_accepted(self):
  last={'objective_gradient_rms_kJ_mol_nm':3.55,'max_constraint_error':0.0}
  self.assertFalse(rr.convergence_observed(last,10.0,1e-7,0.1,3.0))
  self.assertFalse(rr.convergence_observed(last,10.0,1e-7,None,3.0))
 def test_final_projected_force_must_meet_threshold(self):
  last={'objective_gradient_rms_kJ_mol_nm':3.55,'max_constraint_error':2.84e-5}
  self.assertFalse(rr.convergence_observed(last,10.0,3e-8,0.00003,11.0))
 def test_force_projection_on_two_equal_mass_atoms(self):
  from openmm import unit
  system=mm.System();system.addParticle(12.);system.addParticle(1.);system.addConstraint(0,1,.1)
  positions=np.array([[0.,0.,0.],[.1,0.,0.]])*unit.nanometer
  forces=np.array([[1.,0.,0.],[0.,0.,0.]])
  result=rr.projected_force_metrics(system,positions,forces,system.getDefaultPeriodicBoxVectors())
  self.assertAlmostEqual(result['rms_kJ_mol_nm'],np.sqrt(.5/6),places=7)
  self.assertLess(result['max_constraint_normal_kJ_mol_nm'],1e-6)
 def test_actual_constraint_error_is_measured(self):
  system=mm.System();system.addParticle(1);system.addParticle(1);system.addConstraint(0,1,.1)
  pos=np.array([[0.,0.,0.],[.100002,0.,0.]])
  m=rr.measure_constraints(system,pos)
  self.assertEqual(m['constraint_count'],1)
  self.assertAlmostEqual(m['max_relative_error'],2e-5,places=10)
  self.assertEqual(m['violations_above_1e_minus_5'],1)
 def test_real_openmm_final_projection_can_change_constraint_error(self):
  system=mm.System()
  for mass in (16.,1.,1.):system.addParticle(mass)
  for i,j,distance in [(0,1,.09572),(0,2,.09572),(1,2,.15139006545)]:system.addConstraint(i,j,distance)
  force=mm.CustomExternalForce('0.5*k*(x*x+y*y+z*z)+f*x')
  force.addPerParticleParameter('k');force.addPerParticleParameter('f')
  force.addParticle(0,[1000.,0.]);force.addParticle(1,[0.,300.]);force.addParticle(2,[0.,0.]);system.addForce(force)
  integrator=mm.VerletIntegrator(.001)
  context=mm.Context(system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
  context.setPositions([[.02,.03,0],[.12,.03,0],[-.01,.12,0]])
  tracker,attempts=rr.minimize_with_retries(context,10.,1000)
  result=attempts[-1]
  self.assertGreater(tracker.calls,0)
  self.assertGreater(tracker.last['max_constraint_error'],1e-5)
  self.assertLessEqual(result['final_constraints']['max_relative_error'],1e-5)
  self.assertTrue(result['screening_thresholds_observed'])
 def test_water_checkpoint_backup_preserves_originals(self):
  import tempfile,json
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);state=p/'stage0_water_state.xml';state.write_text('state-data');pdb=p/'stage0_water.pdb';pdb.write_text('pdb-data')
   diag=p/'stage0_water_diagnostics.json';diag.write_text(json.dumps({'artifact_sha256':{'state':rr.digest(state),'pdb':rr.digest(pdb)}}))
   rr.backup_previous_water(p)
   backup=p/'before_final_constraint_fix'
   self.assertEqual((backup/state.name).read_text(),'state-data')
   self.assertEqual((backup/pdb.name).read_text(),'pdb-data')
   self.assertEqual(state.read_text(),'state-data')
 def test_corrupt_water_checkpoint_does_not_migrate(self):
  import tempfile,json
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'stage0_water_state.xml').write_text('state');(p/'stage0_water.pdb').write_text('pdb')
   (p/'stage0_water_diagnostics.json').write_text(json.dumps({'artifact_sha256':{'state':'wrong','pdb':'wrong'}}))
   with self.assertRaises(ValueError):rr.backup_previous_water(p)
 def test_migration_only_for_known_unchanged_input(self):
  old={'runner_sha256':rr.PREVIOUS_RUNNER_SHA,'source_sha256':'same','seed':9181}
  new={'runner_sha256':'new','source_sha256':'same','seed':9181}
  self.assertTrue(rr.supported_upgrade(old,new))
  new['source_sha256']='different';self.assertFalse(rr.supported_upgrade(old,new))
  old['runner_sha256']='unrecognized';self.assertFalse(rr.supported_upgrade(old,new))
if __name__=='__main__':unittest.main()
