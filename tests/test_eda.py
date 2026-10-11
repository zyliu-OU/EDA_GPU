from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import pytest
from pyscf import dft, scf
from eda_gpu.eda import calculate_state, frozen_electrostatics, prepare, validate_scope
from eda_gpu.engine import Settings
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def states():
    symbols,coords,m = load(ROOT/"examples/water_dimer.xyz",ROOT/"examples/water_dimer.json")
    settings = Settings(xc="pbe0",basis="sto-3g",density_fitting=False,grid_level=2)
    cm,cd,cs = calculate_state(symbols,coords,m,settings)
    fragments=[]
    for f in m["fragments"]:
        atoms=set(f["atoms"])
        labels=[s if i in atoms else "ghost-"+s for i,s in enumerate(symbols)]
        fragments.append(calculate_state(labels,coords,f,settings))
    return symbols,coords,m,settings,(cm,cd,cs),fragments


@pytest.mark.scientific
def test_frozen_electrostatics_independent_integrals(states):
    _,_,meta,_,(cm,cd,cs),fragments=states
    (am,da,_),(bm,db,_)=fragments
    # Independent dense ERI contraction ONLY for tiny STO-3G validation.
    eri=cm.mol.intor("int2e",aosym="s1")
    coulomb=float(np.einsum("ij,kl,ijkl->",da,db,eri))
    va=am.mol.intor_symmetric("int1e_nuc")
    vb=bm.mol.intor_symmetric("int1e_nuc")
    nn=cm.mol.energy_nuc()-am.mol.energy_nuc()-bm.mol.energy_nuc()
    expected=coulomb+np.einsum("ij,ji->",da,vb)+np.einsum("ij,ji->",db,va)+nn
    result=frozen_electrostatics(cm.mol,da,db,meta["fragments"][0]["atoms"],meta["fragments"][1]["atoms"])
    assert abs(result["electrostatic_hartree"]-expected)<1e-10
    expected_exchange=-.5*np.einsum("ij,kl,ikjl->",da,db,eri)
    assert abs(result["bare_cross_hf_exchange_hartree"]-expected_exchange)<1e-10
    swapped=frozen_electrostatics(cm.mol,db,da,meta["fragments"][1]["atoms"],meta["fragments"][0]["atoms"])
    assert abs(swapped["electrostatic_hartree"]-expected)<1e-10
    assert abs(swapped["bare_cross_hf_exchange_hartree"]-expected_exchange)<1e-10


@pytest.mark.scientific
def test_generalized_correlation_explicit_pbe0_equation(states):
    _,_,_,_,complex_state,fragments=states
    for mf,dm,summary in [complex_state]+fragments:
        _,ex,_=mf._numint.nr_rks(mf.mol,mf.grids,"GGA_X_PBE",dm)
        _,ec,_=mf._numint.nr_rks(mf.mol,mf.grids,"GGA_C_PBE",dm)
        eri=mf.mol.intor("int2e",aosym="s1")
        ex_hf=-.25*np.einsum("ij,kl,ikjl->",dm,dm,eri)
        expected=.75*(ex-ex_hf)+ec
        assert abs(expected-summary["generalized_correlation_hartree"])<1e-8
        assert abs(summary["electron_count"]-mf.mol.nelectron)<1e-8


@pytest.mark.scientific
def test_prepare_never_claims_full_decomposition(states,tmp_path):
    symbols,coords,meta,settings,_,_=states
    r=prepare(symbols,coords,meta,settings,tmp_path)
    assert r["full_eda_available"] is False
    assert r["total_closure_residual_hartree"] is None
    for name in ("exchange","repulsion","polarization"):
        assert r["partial_components_hartree"][name] is None
    refs=r["reference_states"]
    expected=refs["complex"]["energy_hartree"]-sum(f["energy_hartree"] for f in refs["fragments"])
    assert abs(r["counterpoise_interaction_energy_hartree"]-expected)<1e-12
    assert refs["fragments"][0]["ghost_atoms"] == [3,4,5]
    arrays=np.load(tmp_path/"reference_densities.npz")
    assert arrays["complex_density"].shape == arrays["fragment_0_density"].shape


@pytest.mark.parametrize("change",["gpu","density_fitting","functional","open_shell","fragments"])
def test_eda_scope_errors(states,change):
    symbols,coords,meta,settings,_,_=states
    m=json.loads(json.dumps(meta))
    if change=="gpu": settings=replace(settings,backend="gpu")
    if change=="density_fitting": settings=replace(settings,density_fitting=True)
    if change=="functional": settings=replace(settings,xc="cam-b3lyp")
    if change=="open_shell": m["fragments"][0]["spin"]=2
    if change=="fragments": m["fragments"]=m["fragments"][:1]
    with pytest.raises(ValueError): validate_scope(symbols,coords,m,settings)
