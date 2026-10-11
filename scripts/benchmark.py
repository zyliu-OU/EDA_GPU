"""Small-system energy+gradient timings; run after sourcing gpu-env.sh."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import resource
import statistics
import time
from eda_gpu.engine import Settings, hardware, singlepoint
from eda_gpu.model import load


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--repeats", type=int, default=3)
    a = p.parse_args()
    if a.repeats < 2: p.error("at least two repeats required")
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    result = {"hardware": hardware(), "threads": 4, "warmup": "one full energy+gradient call per backend/system excluded", "memory_policy": "process cumulative ru_maxrss and CuPy pool allocation snapshots; not per-operation peak GPU memory", "results": {}}
    for name in ("water", "water_dimer"):
        symbols, coords, m = load(root / f"examples/{name}.xyz", root / f"examples/{name}.json")
        result["results"][name] = {}
        for backend in ("cpu", "gpu"):
            s = replace(Settings(), backend=backend)
            singlepoint(symbols, coords, m, s)
            if backend == "gpu":
                import cupy as cp
                cp.cuda.Device().synchronize()
            timings, energies = [], []
            for _ in range(a.repeats):
                start = time.perf_counter()
                r = singlepoint(symbols, coords, m, s)
                if backend == "gpu": cp.cuda.Device().synchronize()
                timings.append(time.perf_counter()-start)
                energies.append(r["energy_hartree"])
            row = {"seconds": timings, "median_seconds": statistics.median(timings), "energies_hartree": energies,
                   "process_cumulative_maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   "settings": r["settings"]}
            if backend == "gpu":
                row["cupy_pool_used_bytes_after_calls"] = cp.get_default_memory_pool().used_bytes()
                row["cupy_pool_reserved_bytes_after_calls"] = cp.get_default_memory_pool().total_bytes()
                free, total = cp.cuda.runtime.memGetInfo()
                row["device_used_bytes_after_calls"] = total-free
            result["results"][name][backend] = row
    (out / "benchmark.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({name:{b:r["median_seconds"] for b,r in rows.items()} for name,rows in result["results"].items()},indent=2))


if __name__ == "__main__": main()
