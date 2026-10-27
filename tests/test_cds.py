from pathlib import Path
import numpy as np
import pytest
from eda_gpu.cds import association_cds, cds_kcal
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
METHANE = [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]]


def methane_dimer(gap):
    one = np.array(METHANE)
    return ['C', 'H', 'H', 'H', 'H'] * 2, np.r_[one, one + [0, 0, gap]]


def test_separated_fragments_have_no_association_term():
    s, c = methane_dimer(30.0)
    assert abs(association_cds(s, c, [range(5), range(5, 10)])) < 1e-10


def test_methane_contact_is_favourable_in_water():
    s, c = methane_dimer(3.8)
    assert association_cds(s, c, [list(range(5)), list(range(5, 10))]) < -0.1


def test_charged_fragments_and_solvents():
    s, c, m = load(ROOT / "examples/water_dimer.xyz", ROOT / "examples/water_dimer.json")
    groups = [f["atoms"] for f in m["fragments"]]
    assert np.isfinite([association_cds(s, c, groups, solvent) for solvent in ("water", "acetonitrile")]).all()
    assert np.isfinite(cds_kcal(["N", "H", "H", "H", "H"], METHANE))  # odd-electron atom set (NH4 without charge)
    with pytest.raises(ValueError, match="overlap"):
        association_cds(s, c, [groups[0], groups[0]])
    with pytest.raises(ValueError, match="unknown"):
        cds_kcal(s, c, "plasma")
