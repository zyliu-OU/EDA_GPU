"""Construct a reproducible approximate seed; NOT published Cartesian coordinates.
Requires optional rdkit==2026.9.1. DOI 10.1039/C7CP03491A, Table 1/SI Fig S1.
"""
from pathlib import Path
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from eda_gpu.model import xyz_text

out=Path('examples/tryptamine_water');out.mkdir(exist_ok=True)
mol=Chem.AddHs(Chem.MolFromSmiles('NCCc1c[nH]c2ccccc12'))
ids=AllChem.EmbedMultipleConfs(mol,numConfs=100,randomSeed=2017)
energies=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=4,maxIters=500)
# SMILES atoms 0 terminal N, 1 terminal CH2, 2 ring-side CH2,
# 3 indole C3, 4 indole C2, 5 pyrrolic N. All 24 atoms retained.
h2=next(a.GetIdx() for a in mol.GetAtomWithIdx(4).GetNeighbors() if a.GetAtomicNum()==1)
hside=[a.GetIdx() for a in mol.GetAtomWithIdx(2).GetNeighbors() if a.GetAtomicNum()==1]
best=None
for conf_id,(_,energy) in zip(ids,energies):
    c=np.array(mol.GetConformer(conf_id).GetPositions());n=c[0];h=c[h2]
    d=np.linalg.norm(h-n);r=2.88;q=2.39
    if not abs(r-q)<d<r+q:continue
    u=(h-n)/d;t=(r*r-q*q+d*d)/(2*d);radius=np.sqrt(r*r-t*t)
    v=np.cross(u,c[3]-n)
    if np.linalg.norm(v)<1e-8:continue
    v/=np.linalg.norm(v);w=np.cross(u,v)
    for theta in np.linspace(0,2*np.pi,180,endpoint=False):
        oxygen=n+t*u+radius*(v*np.cos(theta)+w*np.sin(theta))
        distances=np.linalg.norm(c-oxygen,axis=1)
        side=min(hside,key=lambda i:distances[i])
        if any(distances[i]<2.5 for i,a in enumerate(mol.GetAtoms()) if a.GetAtomicNum()>1 and i!=0):continue
        score=abs(distances[side]-2.612)*20+energy
        if best is None or score<best[0]:best=(score,c,oxygen,side,energy)
if best is None:raise RuntimeError('no nonclashing bridge seed')
_,c,o,side,energy=best
u=(c[0]-o)/np.linalg.norm(c[0]-o)
v=np.cross(u,c[h2]-o);v/=np.linalg.norm(v)
# Free hydrogen points away from the ring/sidechain.
angle=np.deg2rad(104.5);free=np.cos(angle)*u+np.sin(angle)*v
if np.min(np.linalg.norm(c-(o+.96*free),axis=1))<1.5:free=np.cos(angle)*u-np.sin(angle)*v
coords=np.vstack((c,o,o+.96*u,o+.96*free));symbols=[a.GetSymbol() for a in mol.GetAtoms()]+['O','H','H']
(out/'constructed_bridge.xyz').write_text(xyz_text(symbols,coords,'Constructed RDKit/SI-pattern seed; NOT paper Cartesian coordinates'))
for charge in (0,1):
    meta=dict(identity='tryptamine + water (constructed bridge seed)',protonation='neutral free base; radical ionization for charge +1; no proton added',charge=charge,spin=charge,
        fragments=[dict(name='tryptamine',atoms=list(range(24)),charge=charge,spin=charge),dict(name='water',atoms=[24,25,26],charge=0,spin=0)],
        bonds=[[b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondTypeAsDouble()] for b in mol.GetBonds()]+[[24,25,1],[24,26,1]],
        geometry_provenance='Constructed from RDKit ETKDG/MMFF and the published hydrogen-bond pattern; not extracted from original calculation',
        comparison_atoms=dict(hbond1=[25,0],hbond2=[24,h2],hbond3=[24,side]))
    (out/('neutral.json' if charge==0 else 'cation.json')).write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(dict(atoms=len(symbols),mmff_energy=energy,hbond2_atom=h2,hbond3_atom=side)))
