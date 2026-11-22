import unittest
import numpy as np
import openmm as mm
from openmm import unit, app

class OpenMMChecks(unittest.TestCase):
 def test_forcefields_load(self):
  app.ForceField('amber19/protein.ff19SB.xml','amber14/RNA.OL3.xml','amber19/tip3p.xml')
 def test_actual_minimization_and_restraint(self):
  system=mm.System(); system.addParticle(12.0)
  f=mm.CustomExternalForce('0.5*k*(x-x0)^2'); f.addGlobalParameter('k',4184.0); f.addPerParticleParameter('x0'); f.addParticle(0,[0.0]); system.addForce(f)
  integrator=mm.VerletIntegrator(0.001); ctx=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference')); ctx.setPositions([[0.1,0.,0.]])
  before=ctx.getState(getEnergy=True,getForces=True); e=before.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
  self.assertAlmostEqual(e,20.92,places=6)
  mm.LocalEnergyMinimizer.minimize(ctx,1e-6,100)
  state=ctx.getState(getEnergy=True,getPositions=True); after=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
  self.assertTrue(np.isfinite(after));self.assertLess(after,1e-6)
if __name__=='__main__': unittest.main()
