from dataclasses import replace
import json
from pathlib import Path
import sys
import numpy as np
import pytest
from pyscf import fci, gto, scf
from eda_gpu.cli import main
from eda_gpu.engine import Settings, make_mf
from eda_gpu.gks import (Algebra, reference_state, decompose_states, decompose,
                         orthonormalized_density, hf_like)
from eda_gpu.model import load

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def references():
    symbols,coords,meta=load(ROOT/'examples/water_dimer.xyz',ROOT/'examples/water_dimer.json')
    settings=Settings(xc='pbe0',basis='sto-3g',grid_level=2,density_fitting=False)
    algebra=Algebra('cpu')
    complex_state=reference_state(symbols,coords,meta,settings,algebra)
    fragments=[]
    for f in meta['fragments']:
        labels=[s if i in f['atoms'] else 'ghost-'+s for i,s in enumerate(symbols)]
        fragments.append(reference_state(labels,coords,f,settings,algebra))
    return symbols,coords,meta,settings,complex_state,fragments


@pytest.mark.scientific
def test_each_component_independent_eri_and_qr(references):
    _,_,meta,_,complex_state,fragments=references
    cm,dc,_,cs=complex_state
    (am,da,ca,sa),(bm,db,cb,sb)=fragments
    mol=cm.mol
    eri=mol.intor('int2e',aosym='s1')
    h=mol.intor_symmetric('int1e_kin')+mol.intor_symmetric('int1e_nuc')
    def energy(dm,core=h,nuclear=mol.energy_nuc()):
        return (np.einsum('ij,ji->',dm,core)+.5*np.einsum('ij,kl,ijkl->',dm,dm,eri)
                -.25*np.einsum('ij,kl,ikjl->',dm,dm,eri)+nuclear)
    s=mol.intor_symmetric('int1e_ovlp')
    w,v=np.linalg.eigh(s)
    sqrt=(v*np.sqrt(w))@v.T
    invsqrt=(v/np.sqrt(w))@v.T
    q,_=np.linalg.qr(sqrt@np.concatenate((ca,cb),axis=1))
    dorth=2*(invsqrt@q)@(invsqrt@q).T
    ene=energy(dorth)
    e_sum=energy(da+db)
    ee=np.einsum('ij,kl,ijkl->',da,db,eri)
    ele=ee+np.einsum('ij,ji->',da,bm.mol.intor_symmetric('int1e_nuc'))+np.einsum('ij,ji->',db,am.mol.intor_symmetric('int1e_nuc'))
    ele+=mol.energy_nuc()-am.mol.energy_nuc()-bm.mol.energy_nuc()
    def corr(mf,dm):
        _,ex,_=mf._numint.nr_rks(mf.mol,mf.grids,'GGA_X_PBE',dm)
        _,ec,_=mf._numint.nr_rks(mf.mol,mf.grids,'GGA_C_PBE',dm)
        local_eri=mf.mol.intor('int2e',aosym='s1')
        exact=-.25*np.einsum('ij,kl,ikjl->',dm,dm,local_eri)
        return .75*(ex-exact)+ec
    expected={'electrostatic':ele,'exchange':-.5*np.einsum('ij,kl,ikjl->',da,db,eri),
        'repulsion':ene-e_sum,'polarization':energy(dc)-ene,
        'generalized_correlation':corr(cm,dc)-corr(am,da)-corr(bm,db),'explicit_dispersion':0.}
    result,dasn,_=decompose_states(complex_state,fragments,meta,Algebra('cpu'))
    for k,val in expected.items():
        assert abs(result['components_hartree'][k]-val)<1e-8,k
    assert np.max(np.abs(dasn-dorth))<1e-10
    swapped,_,_=decompose_states(complex_state,list(reversed(fragments)),meta,Algebra('cpu'))
    for k,val in expected.items():
        assert abs(swapped['components_hartree'][k]-val)<1e-8,k


@pytest.mark.scientific
def test_asn_energy_many_electron_determinant():
    # Independent Slater expansion in an orthonormal AO basis, acted on by
    # the full many-electron Hamiltonian. No density-energy formula in oracle.
    mol=gto.M(atom='He 0 0 0; He 0 0 1.8',basis='6-31g',verbose=0)
    a=gto.M(atom='He 0 0 0; ghost-He 0 0 1.8',basis='6-31g',verbose=0)
    b=gto.M(atom='ghost-He 0 0 0; He 0 0 1.8',basis='6-31g',verbose=0)
    ma,mb=scf.RHF(a).run(),scf.RHF(b).run()
    c=np.concatenate((ma.mo_coeff[:,:1],mb.mo_coeff[:,:1]),axis=1)
    s=mol.intor_symmetric('int1e_ovlp');w,v=np.linalg.eigh(s)
    sqrt=(v*np.sqrt(w))@v.T; x=(v/np.sqrt(w))@v.T
    coeff=sqrt@c
    norb=mol.nao
    strings=fci.cistring.gen_strings4orblist(range(norb),2)
    determinants=np.array([np.linalg.det(coeff[[i for i in range(norb) if bits & (1<<i)],:]) for bits in strings])
    ci=np.outer(determinants,determinants);ci/=np.linalg.norm(ci)
    h=x.T@(mol.intor_symmetric('int1e_kin')+mol.intor_symmetric('int1e_nuc'))@x
    eri=np.einsum('ap,bq,cr,ds,abcd->pqrs',x,x,x,x,mol.intor('int2e',aosym='s1'),optimize=True)
    h2=fci.direct_spin1.absorb_h1e(h,eri,norb,(2,2),.5)
    hci=fci.direct_spin1.contract_2e(h2,ci,norb,(2,2))
    expected=float(np.sum(ci*hci))+mol.energy_nuc()
    density,_,_=orthonormalized_density(c,s,Algebra('cpu'))
    actual=hf_like(scf.RHF(mol),density,Algebra('cpu'))['hf_like_total_hartree']
    assert abs(actual-expected)<1e-10
    transform=np.array([[1.2,.3],[-.2,.9]])
    transformed,_,_=orthonormalized_density(c@transform,s,Algebra('cpu'))
    assert np.max(np.abs(transformed-density))<1e-10


def test_dependent_occupied_span_rejected():
    with pytest.raises(ValueError,match='linearly dependent'):
        orthonormalized_density(np.ones((2,2)),np.eye(2),Algebra('cpu'))


def test_gpu_eda_missing_never_falls_back(monkeypatch):
    monkeypatch.setitem(sys.modules,'cupy',None)
    with pytest.raises(RuntimeError,match='CPU fallback is disabled'):
        Algebra('gpu')


@pytest.mark.scientific
def test_cli_complete_is_explicitly_experimental(tmp_path):
    out=tmp_path/'eda'
    code=main(['eda',str(ROOT/'examples/water_dimer.xyz'),'--metadata',str(ROOT/'examples/water_dimer.json'),
        '--backend','cpu','--xc','pbe0','--basis','sto-3g','--grid-level','2','--output',str(out)])
    result=json.loads((out/'result.json').read_text())
    assert code==0 and result['status']=='success'
    assert result['all_components_computed'] is True
    assert result['scientifically_validated'] is True and result['full_eda_available'] is False
    assert all(v is not None for v in result['components_hartree'].values())
    assert abs(result['total_closure_residual_hartree'])<1e-8
    assert result['settings']['density_fitting'] is False
    assert result['settings']['conv_tol'] == 1e-12
    arrays=np.load(out/'reference_densities.npz')
    assert 'asn_density' in arrays and 'fragment_0_occupied' in arrays


@pytest.mark.gpu
@pytest.mark.scientific
def test_gpu_each_component_and_fixed_density(references,tmp_path):
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount()<1: pytest.skip('no CUDA GPU')
        from gpu4pyscf.dft import rks
    except (ImportError,RuntimeError,OSError) as exc:
        pytest.skip(f'GPU support unavailable: {exc}')
    symbols,coords,meta,settings,complex_state,fragments=references
    cpu,_,_=decompose_states(complex_state,fragments,meta,Algebra('cpu'))
    # Same CPU orbitals/densities on GPU isolate J/K and matrix algebra from SCF.
    fixed=[]
    for i,state in enumerate([complex_state]+fragments):
        labels=symbols if i==0 else [s if n in meta['fragments'][i-1]['atoms'] else 'ghost-'+s for n,s in enumerate(symbols)]
        m=meta if i==0 else meta['fragments'][i-1]
        gm=make_mf(labels,coords,m,replace(settings,backend='gpu'))
        fixed.append((gm,cp.asarray(state[1]),cp.asarray(state[2]),state[3]))
    gpu,_,_=decompose_states(fixed[0],fixed[1:],meta,Algebra('gpu'))
    for k,value in cpu['components_hartree'].items():
        assert abs(gpu['components_hartree'][k]-value)<1e-8,k
    end_to_end=decompose(symbols,coords,meta,replace(settings,backend='gpu'),tmp_path)
    for k,value in cpu['components_hartree'].items():
        assert abs(end_to_end['components_hartree'][k]-value)<1e-7,k
    assert end_to_end['backend_used']=='gpu'
    assert end_to_end['execution']['orbital_matrix_algebra']=='gpu'

@pytest.mark.scientific
def test_warm_density_seed_still_converges_actual_scf(references):
    symbols,coords,meta,settings,cold,_=references
    warm=reference_state(symbols,coords,meta,settings,Algebra('cpu'),initial_density=cold[1])
    assert abs(warm[3]['energy_hartree']-cold[3]['energy_hartree'])<1e-9
    assert warm[3]['scf_converged'] and warm[3]['scf_cycles']>=1
    assert warm[3]['scf_initial_guess'].startswith('supplied density')
    with pytest.raises(RuntimeError,match='density'):
        reference_state(symbols,coords,meta,settings,Algebra('cpu'),initial_density=np.zeros_like(cold[1]))
