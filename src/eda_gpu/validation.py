"""Scientific checks using equivalent settings, never guessed reference energies."""
from dataclasses import replace, asdict
from pathlib import Path
import json
import numpy as np
from pyscf.lib import param
from .engine import Settings, singlepoint, optimize
from .model import load


def validate(out, gpu=False):
    root = Path(__file__).resolve().parents[2]
    examples = root / "examples"
    if not examples.exists():
        raise RuntimeError("validation examples require an editable repository installation")
    s = Settings()
    results, checks = {}, {}
    for name in ("water", "water_dimer"):
        symbols, coords, meta = load(examples / f"{name}.xyz", examples / f"{name}.json")
        cpu = singlepoint(symbols, coords, meta, s)
        results[name] = {"cpu": cpu}
        if gpu:
            gr = singlepoint(symbols, coords, meta, replace(s, backend="gpu"))
            de = abs(cpu["energy_hartree"] - gr["energy_hartree"])
            dg = float(np.max(np.abs(np.array(cpu["gradient_hartree_per_bohr"]) - gr["gradient_hartree_per_bohr"])))
            results[name].update(gpu=gr, energy_difference_hartree=de, gradient_max_difference_hartree_per_bohr=dg)
            checks[f"{name}_cpu_gpu"] = de < 1e-7 and dg < 1e-6
        for backend in (["cpu", "gpu"] if gpu else ["cpu"]):
            folder = Path(out) / f"{name}_{backend}_optimization"
            folder.mkdir()
            opt = optimize(symbols, coords, meta, replace(s, backend=backend), folder)
            (folder / "result.json").write_text(json.dumps(opt, indent=2))
            results[name][f"{backend}_optimization"] = {k: v for k, v in opt.items() if k != "evaluations"}
            checks[f"{name}_{backend}_optimization"] = opt["optimization_converged"] and opt["scf_converged"]
    symbols, coords, meta = load(examples / "water.xyz", examples / "water.json")
    analytic = np.array(results["water"]["cpu"]["gradient_hartree_per_bohr"])
    fd = np.zeros_like(coords)
    h = 1e-3 # bohr; all nine Cartesian components, central differences
    for i in range(len(symbols)):
        for j in range(3):
            plus, minus = coords.copy(), coords.copy()
            plus[i,j] += h * param.BOHR
            minus[i,j] -= h * param.BOHR
            ep = singlepoint(symbols, plus, meta, s)["energy_hartree"]
            em = singlepoint(symbols, minus, meta, s)["energy_hartree"]
            fd[i,j] = (ep - em) / (2*h)
    error = float(np.max(np.abs(fd-analytic)))
    results["finite_difference"] = {"step_bohr": h, "gradient_hartree_per_bohr": fd.tolist(), "max_error_hartree_per_bohr": error}
    checks["finite_difference"] = error < 2e-5
    _, _, frozen = load(examples / "water.xyz", examples / "water_frozen.json")
    for backend in (["cpu", "gpu"] if gpu else ["cpu"]):
        folder = Path(out) / f"frozen_{backend}"
        folder.mkdir()
        opt = optimize(symbols, coords, frozen, replace(s, backend=backend), folder)
        (folder / "result.json").write_text(json.dumps(opt,indent=2))
        results[f"frozen_{backend}"] = {k: v for k, v in opt.items() if k != "evaluations"}
        checks[f"frozen_{backend}"] = opt["optimization_converged"] and opt["frozen_max_displacement_angstrom"] < 1e-5
    return {"passed": all(checks.values()), "checks": checks, "results": results,
            "settings": asdict(s), "gpu_validation": "executed" if gpu else "not requested; not validated",
            "tolerances": {"energy_hartree": 1e-7, "gradient_hartree_per_bohr": 1e-6, "finite_difference_hartree_per_bohr": 2e-5, "frozen_angstrom": 1e-5}}
