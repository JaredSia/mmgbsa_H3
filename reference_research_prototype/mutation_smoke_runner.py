"""Test independent ChimeraX A->G/C/U substitutions; not a batch pipeline."""
from pathlib import Path
import json,os,random,shutil,subprocess,sys
import numpy as np
from Bio.PDB import PDBParser
import openmm as mm
from openmm import app,unit
from reference_runner import FF_NAMES,RNA_BB,selections,diagnostics,validate_geometry
from preflight import BASE

def run(reference,destination,exe=None):
 reference=Path(reference).resolve();out=Path(destination).resolve();out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise FileExistsError('Mutation test folder must be empty')
 shutil.copyfile(reference,out/'reference.pdb')
 scripts={}
 for target in ('G','C','U'):
  f=out/f'A_to_{target}.cxc';f.write_text(f'open reference.pdb\nswapna /E:40 {target} preserve true\nsave A_to_{target}_raw.pdb models #1\nexit\n');scripts[target]=f
 if not exe:
  (out/'status.json').write_text(json.dumps({'status':'scripts_generated_only','note':'ChimeraX was not run; run each .cxc from this directory using ChimeraX 1.12 and then validate outputs.'},indent=2))
  print('ChimeraX unavailable: generated scripts only, no structures created.',flush=True);return
 version=subprocess.run([str(exe),'--nogui','--version'],capture_output=True,text=True)
 if version.returncode or 'version: 1.12' not in version.stdout+version.stderr:raise RuntimeError('This smoke test was checked only with ChimeraX 1.12: '+version.stdout+version.stderr)
 P=PDBParser(QUIET=True);ref=P.get_structure('reference',str(reference));ref_protein=[r.resname for r in ref[0]['A']]
 def heavy(s,ch):return {(r.id[1],a.name):a.coord for r in s[0][ch] for a in r if a.element not in ('H','D')}
 ref_p=heavy(ref,'A');ref_r=heavy(ref,'E');ff=app.ForceField(*FF_NAMES);results=[]
 for target in ('G','C','U'):
  try:
   result=subprocess.run([str(exe),'--nogui','--exit',str(scripts[target])],cwd=out,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
   (out/f'A_to_{target}_chimerax.log').write_text(result.stdout)
   p=out/f'A_to_{target}_raw.pdb'
   if result.returncode or not p.exists():raise RuntimeError('ChimeraX failed; inspect log')
   mutant=P.get_structure(target,str(p));rna=mutant[0]['E'];res=rna[40]
   if [r.resname for r in mutant[0]['A']]!=ref_protein:raise ValueError('protein sequence changed')
   if ''.join(r.resname for r in rna)!='AAA'+target+'A'*27:raise ValueError('wrong RNA sequence')
   sugar=RNA_BB;expected=sugar|BASE[target];got={a.name for a in res if a.element not in ('H','D')}
   if got!=expected:raise ValueError('base or backbone heavy atoms do not match expected set')
   p_atoms=heavy(mutant,'A');r_atoms=heavy(mutant,'E');keys=[k for k in ref_r if k[1] in sugar]
   if p_atoms.keys()!=ref_p.keys() or any(k not in r_atoms for k in keys):raise ValueError('protein or RNA backbone atom missing')
   protein_max=max(float(np.linalg.norm(ref_p[k]-p_atoms[k])) for k in ref_p);backbone_max=max(float(np.linalg.norm(ref_r[k]-r_atoms[k])) for k in keys)
   if protein_max>.002 or backbone_max>.002:raise ValueError('protein or RNA backbone coordinates changed')
   bond='N9' if target=='G' else 'N1';bond_length=float(np.linalg.norm(res["C1'"].coord-res[bond].coord))
   if not 1.3<=bond_length<=1.6:raise ValueError('bad glycosidic attachment distance')
   pdb=app.PDBFile(str(p));mod=app.Modeller(pdb.topology,pdb.positions);random.seed(391);mod.addHydrogens(ff,pH=7.4,platform=mm.Platform.getPlatformByName('CPU'))
   unmatched=ff.getUnmatchedResidues(mod.topology)
   if unmatched:raise ValueError('unmatched residues: '+str(unmatched))
   def heavy_coords(top,positions):
    v=np.asarray(positions.value_in_unit(unit.nanometer));return {(a.residue.chain.id,a.residue.id,a.name):v[a.index] for a in top.atoms() if a.element and a.element.symbol not in ('H','D')}
   h0=heavy_coords(pdb.topology,pdb.positions);h1=heavy_coords(mod.topology,mod.positions)
   if h0.keys()!=h1.keys() or any(np.linalg.norm(h0[k]-h1[k])>1e-7 for k in h0):raise ValueError('adding hydrogens changed a pre-existing heavy atom')
   sel=selections(mod.topology);xyz=np.asarray(mod.positions.value_in_unit(unit.nanometer));d0=diagnostics(xyz,xyz,sel);validate_geometry(d0,False)
   system=ff.createSystem(mod.topology,nonbondedMethod=app.NoCutoff,constraints=None,removeCMMotion=False)
   force=mm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)');force.addGlobalParameter('k',4184.)
   for k in ('x0','y0','z0'):force.addPerParticleParameter(k)
   for i in sel['protein_ca']+sel['rna_backbone']:force.addParticle(i,xyz[i])
   system.addForce(force);inte=mm.VerletIntegrator(.001);ctx=mm.Context(system,inte,mm.Platform.getPlatformByName('CPU'),{'Threads':'2'});ctx.setPositions(mod.positions)
   before=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);mm.LocalEnergyMinimizer.minimize(ctx,10,25)
   st=ctx.getState(getPositions=True,getEnergy=True);post=np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer));after=st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
   if not np.isfinite(post).all() or not np.isfinite(after):raise ValueError('nonfinite test results')
   d=diagnostics(post,xyz,sel);validate_geometry(d,False)
   result={'target':target,'status':'smoke_test_passed','protein_heavy_max_displacement_A':protein_max,'rna_backbone_max_displacement_A':backbone_max,'glycosidic_bond_A':bond_length,'unmatched_residues':0,'post_test_stereo_outliers':len(d['stereochemistry_outliers']),'dry_energy_before_kJ_mol':float(before),'dry_energy_after_kJ_mol':float(after),'iterations_cap':25,'limitation':'short vacuum smoke test; neither converged nor a variant for scoring'}
   del ctx,inte
  except Exception as e:result={'target':target,'status':'failed','error':str(e)}
  results.append(result);print(json.dumps(result),flush=True)
  (out/'status.json').write_text(json.dumps({'status':'in_progress','results':results},indent=2))
 status='smoke_tests_finished' if all(r['status']=='smoke_test_passed' for r in results) else 'smoke_tests_failed'
 (out/'status.json').write_text(json.dumps({'status':status,'results':results},indent=2))
 if status!='smoke_tests_finished':raise RuntimeError('One or more smoke tests failed; inspect status.json')
if __name__=='__main__':run(*sys.argv[1:])
