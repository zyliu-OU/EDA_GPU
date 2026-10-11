"""CPU/GPU optimization + common-geometry direct GKS-EDA paper comparison.
Example: python scripts/compare_tryptamine.py --xyz original.xyz --metadata neutral.json
  --case A_H2O --output runs/tryptamine_comparison
No fitted components. A fresh directory is required. Runs CPU/GPU sequentially.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from eda_gpu.model import load, xyz_text
from eda_gpu.engine import Settings, optimize, hardware
from eda_gpu.gks import decompose

EH_TO_KCAL=627.5094740631

def aligned_rmsd(a,b):
    a=a-a.mean(axis=0);b=b-b.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
    return float(np.sqrt(np.mean(np.sum((a@u@correction@vt-b)**2,axis=1))))

def table(result):
    c=result['components_hartree']
    return dict(electrostatic=c['electrostatic']*EH_TO_KCAL,
        pauli=(c['exchange']+c['repulsion'])*EH_TO_KCAL,
        polarization=c['polarization']*EH_TO_KCAL,
        generalized_correlation=c['generalized_correlation']*EH_TO_KCAL,
        total=result['counterpoise_interaction_energy_hartree']*EH_TO_KCAL)

def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xyz',required=True);p.add_argument('--metadata',required=True)
    p.add_argument('--case',choices=['A_H2O','A_cation_H2O','H_cation_H2O'],required=True)
    p.add_argument('--output',required=True);p.add_argument('--grid-level',type=int,default=3)
    p.add_argument('--cart',action='store_true');p.add_argument('--maxsteps',type=int,default=100)
    p.add_argument('--gpu-optimized',help='reuse converged GPU optimized.xyz, with matching result.json')
    p.add_argument('--cpu-optimized',help='reuse converged CPU optimized.xyz with matching result.json')
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
    s,c,m=load(a.xyz,a.metadata)
    reference=json.loads((Path(__file__).resolve().parents[1]/'references/tryptamine_water_2017.json').read_text())
    ref=reference['cases'][a.case]
    if (m['charge'],m['spin'])!=(ref['charge'],ref['spin']):raise ValueError('case charge/spin mismatch')
    settings=Settings(xc='m06-2x',basis='aug-cc-pvdz',grid_level=a.grid_level,density_fitting=False,conv_tol=1e-10,cart=a.cart)
    results=dict(hardware=hardware(),paper_reference=reference,geometry_provenance=m.get('geometry_provenance','user supplied; original provenance not verified'),case=a.case,scientifically_validated=False,optimizations={},common_geometry_eda={})
    coords={}
    for backend in ('gpu','cpu'):
        dest=out/(backend+'_optimization');dest.mkdir()
        reuse=a.gpu_optimized if backend=='gpu' else a.cpu_optimized
        if reuse:
            path=Path(reuse);r=json.loads(path.with_name('result.json').read_text())
            # Ensure reused geometry really comes from identical settings and seed.
            for key,value in vars(replace(settings,backend=backend)).items():
                if r['settings'][key]!=value:raise ValueError('reused optimizer settings mismatch: '+key)
            if not np.allclose(r['initial_coordinates_angstrom'],c,rtol=0,atol=1e-10):raise ValueError('reused initial geometry mismatch')
            _,xyz,_=load(path,a.metadata)
        else:
            r=optimize(s,c,m,replace(settings,backend=backend),dest,a.maxsteps)
            _,xyz,_=load(dest/'optimized.xyz',a.metadata)
        if not r['optimization_converged']:raise RuntimeError(backend+' optimization not converged')
        write(dest/'result.json',r);(dest/'optimized.xyz').write_text(xyz_text(s,xyz))
        results['optimizations'][backend]=r;coords[backend]=xyz
        write(out/'progress.json',results)
    results['aligned_cpu_gpu_geometry_rmsd_angstrom']=aligned_rmsd(coords['cpu'],coords['gpu'])
    # Direct EDA at exactly the SAME GPU-optimized geometry, avoiding confounding
    # geometry and backend differences. CPU optimization result archived separately.
    for backend in ('cpu','gpu'):
        dest=out/(backend+'_common_geometry_eda');dest.mkdir()
        r=decompose(s,coords['gpu'],m,replace(settings,backend=backend,conv_tol=1e-12),dest)
        write(dest/'result.json',r);results['common_geometry_eda'][backend]=dict(result=r,components_kcal_mol=table(r))
        write(out/'progress.json',results)
    cpu=results['common_geometry_eda']['cpu']['components_kcal_mol'];gpu=results['common_geometry_eda']['gpu']['components_kcal_mol']
    results['gpu_minus_cpu_kcal_mol']={k:gpu[k]-cpu[k] for k in cpu}
    results['gpu_minus_paper_kcal_mol']={k:gpu[k]-ref[k] for k in gpu}
    results['component_agreement_tolerance_hartree']=1e-7
    results['cpu_gpu_eda_passed']=max(abs(v) for v in results['gpu_minus_cpu_kcal_mol'].values())<1e-7*EH_TO_KCAL
    results['hydrogen_bond_distances_angstrom']={key:float(np.linalg.norm(coords['gpu'][i]-coords['gpu'][j])) for key,(i,j) in m.get('comparison_atoms',{}).items()}
    results['paper_comparison_status']='exploratory: original Cartesian coordinates and numerical settings not verified; no D0, Hessian, intra-EDA or individual H-bond decomposition'
    write(out/'result.json',results)
    return 0 if results['cpu_gpu_eda_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
