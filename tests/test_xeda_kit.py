from pathlib import Path
import importlib.util
import pytest
from pyscf.data.nist import BOHR
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("xeda_matched", ROOT / "scripts/xeda_matched.py")
xeda = importlib.util.module_from_spec(spec)
spec.loader.exec_module(xeda)

TABLE = """ 1   EDA RESULTS
 3    ALL BASIS SET                                       HARTREE          KCAL/MOL
 5    ELECTROSTATIC ENERGY              ES=              -0.026946          -16.91
 6    EXCHANGE ENERGY                   EX=              -0.055228          -34.66
 7    REPULSION ENERGY                 REP=               0.094078           59.03
 8    POLARIZATION ENERGY              POL=              -0.002168           -1.36
 9    ELECTRON CORRELATION            CORR=              -0.014048           -8.82
10    TOTAL INTERACTION ENERGY           E=              -0.004311           -2.71
"""


def test_parse_log_reads_table_and_ignores_unrelated_assignments(tmp_path):
    log = tmp_path / "case.log"
    log.write_text(" ITER  DEL-E=  -0.5000001  E= -152.1\n" + TABLE)
    assert xeda.parse_log(log) == {"electrostatic": -0.026946, "exchange": -0.055228, "repulsion": 0.094078,
                                   "polarization": -0.002168, "generalized_correlation": -0.014048, "total": -0.004311}


def test_gamess_input_uses_bohr_fragments_and_explicit_representation():
    s, c, m = load(ROOT / "examples/water_dimer.xyz", ROOT / "examples/water_dimer.json")
    text = xeda.gamess_input("t", s, c, m, "PBE0", "GBASIS=CCD", cart=False)
    assert "MATOM(1)=3 3" in text and "EDATYP=GKS SUPBAS=.T." in text and "ISPHER=1" in text and "UNITS=BOHR" in text
    oxygen_b = [line for line in text.splitlines() if line.startswith(" O")][1].split()
    assert float(oxygen_b[2]) == pytest.approx(c[3][0] / BOHR, abs=1e-11)


def test_gamess_input_rejects_noncontiguous_fragments():
    s, c, m = load(ROOT / "examples/water_dimer.xyz", ROOT / "examples/water_dimer.json")
    m["fragments"][0]["atoms"], m["fragments"][1]["atoms"] = [0, 1, 3], [2, 4, 5]
    with pytest.raises(ValueError, match="contiguous"):
        xeda.gamess_input("t", s, c, m, "PBE0", "GBASIS=CCD", cart=False)
