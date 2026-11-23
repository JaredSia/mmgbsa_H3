import unittest,sys,tempfile,json
from pathlib import Path
from unittest.mock import patch,MagicMock
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import reference_runner as rr

class PrecisionFixTests(unittest.TestCase):
 def test_cuda_requests_mixed_precision(self):
  platform=MagicMock(); context=MagicMock()
  with patch.object(rr.mm.Platform,'getPlatformByName',return_value=platform), patch.object(rr.mm,'Context',return_value=context) as create:
   found,name,errors=rr.context_for('system','integrator','CUDA')
  self.assertEqual(name,'CUDA');self.assertEqual(create.call_args.args[3],{'Precision':'mixed'})
 def test_final_tangent_gradient_controls_not_provisional_gradient(self):
  last={'objective_gradient_rms_kJ_mol_nm':364.29,'max_constraint_error':4.16e-5}
  self.assertTrue(rr.convergence_observed(last,10.0,1e-7,6.1e-5,5.78))
 def test_nonfinite_reporter_gradient_is_not_accepted(self):
  last={'objective_gradient_rms_kJ_mol_nm':float('nan'),'max_constraint_error':0.0}
  self.assertFalse(rr.convergence_observed(last,10.0,1e-7,6e-5,5.0))
 def test_final_constraint_still_must_pass(self):
  last={'objective_gradient_rms_kJ_mol_nm':364.29,'max_constraint_error':4.16e-5}
  self.assertFalse(rr.convergence_observed(last,10.0,1.68e-5,6.1e-5,5.78))
 def test_final_force_still_must_pass(self):
  last={'objective_gradient_rms_kJ_mol_nm':1.0,'max_constraint_error':1e-7}
  self.assertFalse(rr.convergence_observed(last,10.0,1e-7,6e-5,10.01))
 def test_final_projected_force_matches_constraint_jacobian(self):
  import numpy as np
  import openmm as mm
  from openmm import unit
  xyz=np.array([[0.,0.,0.],[.10,0.,0.],[.03,.09,0.]])
  forces=np.array([[12.,-3.,5.],[-2.,7.,1.],[4.,1.,-6.]])
  system=mm.System()
  for mass in (16.,1.,1.):system.addParticle(mass)
  pairs=[(0,1),(1,2),(0,2)]
  for i,j in pairs:system.addConstraint(i,j,float(np.linalg.norm(xyz[i]-xyz[j])))
  jac=np.zeros((3,9))
  for k,(i,j) in enumerate(pairs):
   direction=(xyz[i]-xyz[j])/np.linalg.norm(xyz[i]-xyz[j]);jac[k,3*i:3*i+3]=direction;jac[k,3*j:3*j+3]=-direction
  flat=forces.reshape(-1);expected=flat-jac.T@np.linalg.solve(jac@jac.T,jac@flat)
  result=rr.projected_force_metrics(system,xyz*unit.nanometer,forces,system.getDefaultPeriodicBoxVectors())
  self.assertAlmostEqual(result['rms_kJ_mol_nm'],float(np.sqrt(np.mean(expected**2))),places=7)
  self.assertLess(result['max_constraint_normal_kJ_mol_nm'],1e-6)
 def test_upgrade_from_current_known_runner(self):
  old={'runner_sha256':'e575871e39f132e635baa8f8d5131a43b957c0899ce4fc497ac26450ca622e69','source_sha256':'same','seed':9181}
  new={'runner_sha256':'new','source_sha256':'same','seed':9181}
  self.assertTrue(rr.supported_upgrade(old,new))
  new['source_sha256']='changed';self.assertFalse(rr.supported_upgrade(old,new))
 def test_new_backup_does_not_overwrite_previous_backup(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);prev=p/'before_final_constraint_fix';prev.mkdir();(prev/'sentinel').write_text('preserved')
   state=p/'stage0_water_state.xml';state.write_text('new-state');pdb=p/'stage0_water.pdb';pdb.write_text('new-pdb')
   (p/'stage0_water_diagnostics.json').write_text(json.dumps({'artifact_sha256':{'state':rr.digest(state),'pdb':rr.digest(pdb)}}))
   rr.backup_previous_water(p,tag='cuda_precision_fix')
   self.assertEqual((prev/'sentinel').read_text(),'preserved')
   self.assertEqual((p/'before_cuda_precision_fix'/state.name).read_text(),'new-state')
if __name__=='__main__':unittest.main()
