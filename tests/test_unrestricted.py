from dataclasses import replace
import numpy as np
import pytest
from eda_gpu.engine import Settings, singlepoint, make_mf
from eda_gpu.gks import Algebra, decompose, hf_like
from eda_gpu.model import load
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def radical():
    s,c,m=load(ROOT/'examples/water_dimer.xyz',ROOT/'examples/water_dimer.json')
    m['charge']=1;m['spin']=1;m['fragments'][0]['charge']=1;m['fragments'][0]['spin']=1
    return s,c,m

@pytest.mark.scientific
def test_unrestricted_exact_exchange_and_gradient(tmp_path):
    s,c,m=radical();settings=Settings(xc='m06-2x',basis='sto-3g',grid_level=2,density_fitting=False,conv_tol=1e-12)
    r=decompose(s,c,m,settings,tmp_path)
    assert r['stage']=='experimental_complete_unrestricted_gks_eda'
    assert abs(r['total_closure_residual_hartree'])<1e-8
    arrays=np.load(tmp_path/'reference_densities.npz');dm=arrays['complex_density']
    mf=make_mf(s,c,m,settings);eri=mf.mol.intor('int2e',aosym='s1');h=mf.get_hcore();total=dm.sum(axis=0)
    expected=np.einsum('ij,ji->',total,h)+.5*np.einsum('ij,kl,ijkl->',total,total,eri)
    expected-=.5*sum(np.einsum('ij,kl,ikjl->',d,d,eri) for d in dm)
    expected+=mf.mol.energy_nuc()
    assert abs(expected-r['reference_states']['complex']['hf_like_total_hartree'])<1e-9
    # Independent spin-resolved cross exchange contracts four-index ERIs.
    da,db=arrays['fragment_0_density'],arrays['fragment_1_density']/2
    ex=-sum(np.einsum('ij,kl,ikjl->',da[i],db,eri) for i in range(2))
    assert abs(ex-r['components_hartree']['exchange'])<1e-9
    overlap=arrays['overlap'];w,v=np.linalg.eigh(overlap)
    sqrt=(v*np.sqrt(w))@v.T;invsqrt=(v/np.sqrt(w))@v.T
    frozen=[]
    for spin in ('alpha','beta'):
        ca=arrays['fragment_0_occupied_'+spin];cb=arrays['fragment_1_occupied']
        q,_=np.linalg.qr(sqrt@np.concatenate((ca,cb),axis=1))
        coeff=invsqrt@q;frozen.append(coeff@coeff.T)
    def energy(d):
        t=np.sum(d,axis=0)
        return np.einsum('ij,ji->',t,h)+.5*np.einsum('ij,kl,ijkl->',t,t,eri)-.5*sum(np.einsum('ij,kl,ikjl->',x,x,eri) for x in d)+mf.mol.energy_nuc()
    assert np.max(np.abs(np.array(frozen)-arrays['asn_density']))<1e-10
    as_energy=energy(da+np.stack((db,db)))
    assert abs(energy(frozen)-as_energy-r['components_hartree']['repulsion'])<1e-9
    assert abs(energy(dm)-energy(frozen)-r['components_hartree']['polarization'])<1e-9
    # Analytic open-shell meta-GGA gradient checked against displaced energies.
    sp=singlepoint(s,c,m,settings);step=1e-3;delta=np.zeros_like(c);delta[0,2]=step*.529177210903
    energies=[]
    for xyz in (c+delta,c-delta):
        f=make_mf(s,xyz,m,settings);energies.append(float(f.kernel()));assert f.converged
    fd=(energies[0]-energies[1])/(2*step)
    assert abs(fd-sp['gradient_hartree_per_bohr'][0][2])<2e-4

@pytest.mark.gpu
@pytest.mark.scientific
def test_unrestricted_gpu_cpu_agreement(tmp_path):
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount()<1: pytest.skip('no CUDA device')
    except (ImportError,RuntimeError,OSError):
        pytest.skip('GPU unavailable')
    s,c,m=radical();settings=Settings(xc='m06-2x',basis='sto-3g',grid_level=2,density_fitting=False,conv_tol=1e-12)
    cpu_dir=tmp_path/'cpu';gpu_dir=tmp_path/'gpu';cpu_dir.mkdir();gpu_dir.mkdir()
    cpu=decompose(s,c,m,settings,cpu_dir)
    gpu=decompose(s,c,m,replace(settings,backend='gpu'),gpu_dir)
    for k,v in cpu['components_hartree'].items():
        assert abs(v-gpu['components_hartree'][k])<2e-7,k
    a=singlepoint(s,c,m,settings);b=singlepoint(s,c,m,replace(settings,backend='gpu'))
    assert abs(a['energy_hartree']-b['energy_hartree'])<1e-7
    assert np.max(np.abs(np.array(a['gradient_hartree_per_bohr'])-b['gradient_hartree_per_bohr']))<1e-6
