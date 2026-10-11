from pathlib import Path
import json
import sys
import pytest
from eda_gpu.cli import main
from eda_gpu.engine import Settings, singlepoint
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]


def test_gpu_missing_never_falls_back(monkeypatch):
    monkeypatch.setitem(sys.modules, "cupy", None)
    s,r,m = load(ROOT / "examples/water.xyz", ROOT / "examples/water.json")
    with pytest.raises(RuntimeError, match="CPU fallback is disabled"):
        singlepoint(s,r,m,Settings(backend="gpu"))


def test_cli_failed_scf_records_failure(tmp_path):
    out = tmp_path / "failed"
    code = main(["singlepoint", str(ROOT / "examples/water.xyz"), "--metadata", str(ROOT / "examples/water.json"), "--backend", "cpu", "--max-cycle", "1", "--output", str(out)])
    result = json.loads((out / "result.json").read_text())
    assert code == 1 and result["status"] == "failed"
    assert "SCF did not converge" in result["error"]
    before = (out / "result.json").read_text()
    assert main(["singlepoint", "missing.xyz", "--metadata", "missing.json", "--backend", "cpu", "--output", str(out)]) == 2
    assert (out / "result.json").read_text() == before


@pytest.mark.gpu
def test_real_gpu_energy_gradient():
    # Honest skip for CPU CI; actual build host must execute this marked test.
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1: pytest.skip("no CUDA GPU")
        from gpu4pyscf.dft import rks
    except (ImportError, RuntimeError, OSError) as e:
        pytest.skip(f"GPU support unavailable: {e}")
    s,r,m = load(ROOT / "examples/water.xyz", ROOT / "examples/water.json")
    result = singlepoint(s,r,m,Settings(backend="gpu"))
    assert result["backend_used"] == "gpu" and result["scf_converged"]
    assert len(result["gradient_hartree_per_bohr"]) == 3


@pytest.mark.scientific
def test_optimizer_iteration_limit_is_not_success(tmp_path):
    out = tmp_path / "limited"
    code = main(["optimize", str(ROOT / "examples/water.xyz"), "--metadata", str(ROOT / "examples/water.json"), "--backend", "cpu", "--maxsteps", "1", "--output", str(out)])
    result = json.loads((out / "result.json").read_text())
    assert code == 1 and result["status"] == "not_converged"
    assert result["optimization_converged"] is False
    assert result["scf_converged"] is True
    assert (out / "trajectory.xyz").exists()
