import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
from .model import load, xyz_text
from .engine import Settings, hardware, singlepoint, optimize


def main(argv=None):
    p = argparse.ArgumentParser(prog="eda-gpu", description="Molecular DFT and experimental complete CPU/GPU GKS-EDA; external validation pending.")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("inspect")
    v = sub.add_parser("validate")
    v.add_argument("--output", default="runs/validation")
    v.add_argument("--gpu", action="store_true", help="require actual GPU comparisons")
    for name in ("singlepoint", "optimize", "eda-prepare", "eda"):
        q = sub.add_parser(name, help=("experimental CPU fragment references and partial components" if name == "eda-prepare" else ("experimental complete restricted/unrestricted GKS-EDA; external validation pending" if name == "eda" else None)))
        q.add_argument("xyz")
        q.add_argument("--metadata", required=True)
        q.add_argument("--backend", choices=(["cpu"] if name == "eda-prepare" else ["cpu", "gpu"]), required=True)
        q.add_argument("--output", required=True)
        q.add_argument("--xc", default="b3lyp")
        q.add_argument("--basis", default="def2-svp")
        q.add_argument("--cart", action="store_true", help="use Cartesian basis functions; default spherical")
        q.add_argument("--auxbasis", default="def2-universal-jkfit")
        q.add_argument("--grid-level", type=int, default=3)
        q.add_argument("--conv-tol", type=float, default=(1e-12 if name == "eda" else 1e-10))
        q.add_argument("--max-cycle", type=int, default=100)
        q.add_argument("--threads", type=int, default=4)
        q.add_argument("--no-density-fitting", action="store_true")
        if name == "optimize":
            q.add_argument("--maxsteps", type=int, default=100)
            q.add_argument("--constraints")
    a = p.parse_args(argv)
    if a.command == "inspect":
        print(json.dumps(hardware(), indent=2))
        return 0
    out = Path(a.output)
    # Avoid overwriting existing calculations.
    try:
        out.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"output directory already exists: {out}", file=sys.stderr)
        return 2
    result = {"schema_version": 1, "command": a.command, "hardware": hardware(), "status": "failed"}
    try:
        if a.command == "validate":
            from .validation import validate
            result.update(validate(out, a.gpu))
        else:
            symbols, coords, metadata = load(a.xyz, a.metadata)
            settings = Settings(backend=a.backend, xc=a.xc, basis=a.basis, auxbasis=a.auxbasis,
                                grid_level=a.grid_level, conv_tol=a.conv_tol, max_cycle=a.max_cycle,
                                threads=a.threads, cart=a.cart, density_fitting=(False if a.command in ("eda-prepare", "eda") else not a.no_density_fitting))
            result.update(metadata=metadata, symbols=symbols, initial_coordinates_angstrom=coords.tolist(), settings=asdict(settings), backend_requested=a.backend)
            (out / "initial.xyz").write_text(xyz_text(symbols, coords))
            if a.command == "eda":
                from .gks import decompose
                calc = decompose(symbols, coords, metadata, settings, out)
            elif a.command == "eda-prepare":
                from .eda import prepare
                calc = prepare(symbols, coords, metadata, settings, out)
            else:
                calc = optimize(symbols, coords, metadata, settings, out, a.maxsteps, a.constraints) if a.command == "optimize" else singlepoint(symbols, coords, metadata, settings)
            result.update(calc)
        result["status"] = "success" if result.get("optimization_converged", result.get("passed", True)) else "not_converged"
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"
    (out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "result": str(out / "result.json"), "error": result.get("error"), **({"full_eda_available": False, "all_components_computed": result.get("all_components_computed", False)} if a.command in ("eda-prepare", "eda") else {})}, indent=2))
    return 0 if result["status"] == "success" else 1
