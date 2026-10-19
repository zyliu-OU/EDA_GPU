from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from eda_gpu.engine import Settings, make_mf
from eda_gpu.gks import decompose, validate_gks_scope
from eda_gpu.model import load
from eda_gpu.solvent import surface_nuclear_potential

ROOT = Path(__file__).resolve().parents[1]
FAST = Settings(backend="cpu", xc="pbe0", basis="sto-3g", density_fitting=False, grid_level=1, conv_tol=1e-11,
                conv_tol_grad=1e-6, threads=2, cpu_core_limit=0, compute_duty_cycle=1.0)


def dimer():
    return load(ROOT / "examples/water_dimer.xyz", ROOT / "examples/water_dimer.json")


def test_solvent_scope_and_dielectric_are_checked():
    s, c, m = dimer()
    with pytest.raises(ValueError, match="solvent must be"):
        validate_gks_scope(s, c, m, replace(FAST, solvent="smd"))
    with pytest.raises(ValueError, match="dielectric"):
        make_mf(s, c, m, replace(FAST, solvent="cpcm", solvent_eps=1.0))


def test_monomer_uses_complex_cavity_with_own_nuclei_only():
    s, c, m = dimer()
    complex_mf = make_mf(s, c, m, replace(FAST, solvent="cpcm"))
    f = m["fragments"][0]
    labels = [x if i in f["atoms"] else "ghost-" + x for i, x in enumerate(s)]
    monomer = make_mf(labels, c, f, replace(FAST, solvent="cpcm"), cavity_mol=complex_mf.mol)
    pcm = monomer.with_solvent
    assert pcm.mol.natm == 6 and len(pcm.surface["grid_coords"]) == len(complex_mf.with_solvent.surface["grid_coords"])
    own = surface_nuclear_potential(complex_mf.mol, pcm.surface, monomer.mol.atom_charges())
    assert np.allclose(pcm.v_grids_n, own, atol=1e-12)
    assert not np.allclose(pcm.v_grids_n, complex_mf.with_solvent.v_grids_n)


def test_gks_eda_sol_closure_and_vacuum_limit(tmp_path):
    s, c, m = dimer()
    for d in ("gas", "sol", "eps1"):
        (tmp_path / d).mkdir()
    gas = decompose(s, c, m, FAST, tmp_path / "gas")
    sol = decompose(s, c, m, replace(FAST, solvent="cpcm"), tmp_path / "sol")
    near = decompose(s, c, m, replace(FAST, solvent="cpcm", solvent_eps=1.000001), tmp_path / "eps1")
    assert sol["stage"].endswith("_sol") and "desolvation" in sol["components_hartree"]
    assert gas["scientifically_validated"] is True and sol["scientifically_validated"] is False
    assert abs(sum(sol["components_hartree"].values()) - sol["counterpoise_interaction_energy_hartree"]) < 1e-8
    assert abs(sol["components_hartree"]["desolvation"]) > 1e-5
    for k, v in gas["components_hartree"].items():
        assert abs(near["components_hartree"][k] - v) < 1e-6
    assert abs(near["components_hartree"]["desolvation"]) < 1e-6
