import json
from pathlib import Path
import pytest
from eda_gpu.model import load
from eda_gpu.engine import Settings, singlepoint

ROOT = Path(__file__).resolve().parents[1]


def test_valid_partition():
    s, r, m = load(ROOT / "examples/water_dimer.xyz", ROOT / "examples/water_dimer.json")
    assert len(s) == 6 and r.shape == (6,3)
    assert m["charge"] == 0


@pytest.mark.parametrize("case", ["spin", "duplicate", "missing", "charge", "noninteger", "frozen"])
def test_bad_metadata(tmp_path, case):
    m = json.loads((ROOT / "examples/water.json").read_text())
    if case == "spin": m["spin"] = 1
    if case == "duplicate": m["fragments"][0]["atoms"] = [0,1,1,2]
    if case == "missing": m["fragments"][0]["atoms"] = [0,1]
    if case == "charge": m["fragments"][0]["charge"] = 2
    if case == "noninteger": m["charge"] = 0.0
    if case == "frozen": m["frozen_atoms"] = [3]
    f = tmp_path / "metadata.json"
    f.write_text(json.dumps(m))
    with pytest.raises(ValueError): load(ROOT / "examples/water.xyz", f)


@pytest.mark.scientific
def test_scf_failure_is_explicit():
    s,r,m = load(ROOT / "examples/water.xyz", ROOT / "examples/water.json")
    with pytest.raises(RuntimeError, match="SCF did not converge"):
        singlepoint(s,r,m,Settings(max_cycle=1))
