import json, sys, glob, re
from pathlib import Path
import numpy as np
from itertools import combinations
cov = {'H':0.31,'C':0.76,'N':0.71,'O':0.66}
def rx(p):
    L=Path(p).read_text().split('\n'); n=int(L[0]); s=[];c=[]
    for l in L[2:2+n]: t=l.split(); s.append(t[0]); c.append([float(x) for x in t[1:4]])
    return s,np.array(c)
rows={}
for run, seeddir in [('runs/hplc_patch_cpcm_def2svp','examples/hplc'), ('runs/hplc_patch_sampling_cpcm_def2svp','examples/hplc/patch_sampling')]:
    np_=json.load(open(f'{run}/nonpolar.json'))['cases']
    for d in sorted(glob.glob(f'{run}/patch_*_[0-9]*')):
        name=Path(d).name
        meta=json.load(open(next(p for p in [f'{seeddir}/{name}.json', f'examples/hplc/patch/{name}.json', *glob.glob(f'examples/hplc/**/{name}.json',recursive=True)] if Path(p).exists())))
        o=json.load(open(f'{d}/optimization/result.json')); e=json.load(open(f'{d}/eda/result.json'))
        s,X=rx(f'{d}/optimization/optimized.xyz')
        bonds=[tuple(b[:2]) for b in meta['bonds']]
        bset=set(map(lambda b: tuple(sorted(b)), bonds))
        r=np.linalg.norm(X[:,None]-X[None],axis=-1)
        broken=[b for b in bonds if r[b]>1.25*(cov[s[b[0]]]+cov[s[b[1]]])]
        new=[(i,j) for i,j in combinations(range(len(s)),2) if (i,j) not in bset and r[i,j]<1.15*(cov[s[i]]+cov[s[j]]) ]
        nh=78
        C1=X[[0,26,52]]; cen=C1.mean(0); nrm=np.cross(C1[1]-C1[0],C1[2]-C1[0]); nrm/=np.linalg.norm(nrm)
        if np.dot(X[[7,33,59]].mean(0)-cen,nrm)<0: nrm=-nrm
        N=[i for i in range(nh,len(s)) if s[i]=='N'][0]; v=X[N]-cen; h=v@nrm; lat=np.linalg.norm(v-h*nrm)
        guestmin=min(((X[i]-cen)@nrm) for i in range(nh,len(s)) if s[i]!='H')
        pc=[i for i in range(nh) if s[i]=='C']; gc=[i for i in range(nh,len(s)) if s[i]=='C']
        cc=int((r[np.ix_(gc,pc)]<4.5).sum())
        comp={k:v*627.509474 for k,v in e['components_hartree'].items()}
        npc=next(v for k,v in np_.items() if k.endswith('/'+name)); sp=next(iter(npc['splits'].values()))
        rows[name]=dict(conv=o['optimization_converged'], scf=e['scf_converged'], frz=o['frozen_max_displacement_angstrom'], clos=abs(e['total_closure_residual_hartree']),
          broken=len(broken), new=len(new), E=o['energy_hartree'], evals=len(o['evaluations']), wall=o['seconds']/3600, h=h, r=lat, gmin=guestmin, cc=cc,
          elec=comp['electrostatic'], pauli=comp['exchange']+comp['repulsion'], pol=comp['polarization'], gc=comp['generalized_correlation'], des=comp['desolvation'],
          tot=sum(comp.values()), cdsw=sp['dG_cds_water'], epc=sp['eda_plus_cds_water'], cdsm=sp['dG_cds_acetonitrile'], np=sp['buried_nonpolar_A2'], run=run)
json.dump(rows, open(sys.argv[1],'w'), indent=1)
for k,v in rows.items(): print(k, {a:(round(b,3) if isinstance(b,float) else b) for a,b in v.items() if a!='run'})
