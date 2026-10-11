import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


def main(argv=None):
    p = argparse.ArgumentParser(prog="eda-gpu", description="Molecular DFT and experimental complete CPU/GPU GKS-EDA; external validation pending.")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("inspect")
    for name in ("monitor","status"):
        q=sub.add_parser(name,help="read live calculation progress and results")
        q.add_argument("run_folder")
        q.add_argument("--log",help="optional live log for runs started before progress recording")
        q.add_argument("--interval",type=float,default=1.0)
        q.add_argument("--json",action="store_true")
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
        q.add_argument("--scf-gradient-tol",type=float,default=None)
        q.add_argument("--max-cycle", type=int, default=100)
        q.add_argument("--threads", type=int, default=16)
        q.add_argument("--cpu-core-limit", type=int, default=16, help="process affinity cap on Linux; 0 disables affinity cap")
        q.add_argument("--compute-duty-cycle", type=float, default=0.5, help="process workload pacing target, default 0.5; not an instantaneous GPU utilization cap")
        q.add_argument("--gpu-temperature-limit", type=float, default=75.0)
        q.add_argument("--cpu-temperature-limit", type=float, default=85.0)
        q.add_argument("--no-density-fitting", action="store_true")
        if name == "eda":
            q.add_argument("--reference-seed",help="NPZ reference densities at this exact geometry/basis; actual SCFs still run")
        if name == "optimize":
            q.add_argument("--maxsteps", type=int, default=100)
            q.add_argument("--constraints")
    a = p.parse_args(argv)
    if a.command == "inspect":
        from .engine import hardware
        print(json.dumps(hardware(), indent=2))
        return 0
    if a.command in ("monitor","status"):
        from .progress import monitor
        return monitor(a.run_folder,a.log,a.interval,once=a.command=="status",json_output=a.json)
    from .model import load, xyz_text
    from .engine import Settings, hardware, singlepoint, optimize
    out = Path(a.output)
    # Avoid overwriting existing calculations.
    try:
        out.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"output directory already exists: {out}", file=sys.stderr)
        return 2
    from .citations import information, text
    print(text(),flush=True)
    result = {**information(),"schema_version": 1, "command": a.command, "hardware": hardware(), "status": "failed"}
    try:
        if a.command == "validate":
            from .validation import validate
            result.update(validate(out, a.gpu))
        else:
            symbols, coords, metadata = load(a.xyz, a.metadata)
            settings = Settings(backend=a.backend, xc=a.xc, basis=a.basis, auxbasis=a.auxbasis,
                                grid_level=a.grid_level, conv_tol=a.conv_tol,conv_tol_grad=a.scf_gradient_tol, max_cycle=a.max_cycle,
                                threads=a.threads, cpu_core_limit=a.cpu_core_limit,cart=a.cart, compute_duty_cycle=a.compute_duty_cycle,
                                gpu_temperature_limit=a.gpu_temperature_limit,cpu_temperature_limit=a.cpu_temperature_limit, density_fitting=(False if a.command in ("eda-prepare", "eda") else not a.no_density_fitting))
            result.update(information(metadata))
            result.update(metadata=metadata, symbols=symbols, initial_coordinates_angstrom=coords.tolist(), settings=asdict(settings), backend_requested=a.backend)
            (out / "initial.xyz").write_text(xyz_text(symbols, coords))
            if a.command == "eda":
                from .gks import decompose
                seeds=None
                if a.reference_seed:
                    import numpy as np
                    with np.load(a.reference_seed,allow_pickle=False) as archive:
                        seeds={key:archive[key] for key in ("complex_density","fragment_0_density","fragment_1_density")}
                    result["reference_seed_file"]=a.reference_seed
                calc = decompose(symbols, coords, metadata, settings, out,initial_densities=seeds)
            elif a.command == "eda-prepare":
                from .eda import prepare
                calc = prepare(symbols, coords, metadata, settings, out)
            else:
                calc = optimize(symbols, coords, metadata, settings, out, a.maxsteps, a.constraints) if a.command == "optimize" else singlepoint(symbols, coords, metadata, settings,out=out)
            result.update(calc)
        result["status"] = "success" if result.get("optimization_converged", result.get("passed", True)) else "not_converged"
    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}"
    from .progress import atomic_json, ProgressWriter
    atomic_json(out / "result.json",result)
    ProgressWriter(out).update("finished",status=result["status"])
    print(json.dumps({"status": result["status"], "result": str(out / "result.json"), "error": result.get("error"), **({"full_eda_available": False, "all_components_computed": result.get("all_components_computed", False)} if a.command in ("eda-prepare", "eda") else {})}, indent=2))
    return 0 if result["status"] == "success" else 1
