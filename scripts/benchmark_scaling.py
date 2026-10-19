#!/usr/bin/env python
"""Repeated CPU/GPU timings across system sizes for the two production paths.

singlepoint: density-fitted energy + analytic gradient (the optimizer's evaluation).
eda:         direct-integral counterpoise GKS-EDA (all reference and intermediate SCFs).
One warm-up call per case/backend is excluded; every repetition reruns all SCFs.
Peak GPU memory is sampled from nvidia-smi every 0.25 s during each call (whole device,
so other processes are included; the idle baseline is recorded). Run after sourcing
scripts/gpu-env.sh. Timings describe this host only.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import statistics
import subprocess
import threading
import time
from eda_gpu.engine import Settings, hardware, singlepoint
from eda_gpu.gks import decompose
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
TRYP = (ROOT / 'reports/tryptamine_neutral_gpu_2026-10-11/initial.xyz', ROOT / 'examples/tryptamine_water/neutral.json')
CASES = [  # name, xyz, metadata, xc, basis, auxbasis
    ('water_dimer_def2svp', ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json', 'b3lyp', 'def2-svp', 'def2-universal-jkfit'),
    ('tryptamine_water_def2svp', *TRYP, 'm06-2x', 'def2-svp', 'def2-universal-jkfit'),
    ('tryptamine_water_augccpvdz', *TRYP, 'm06-2x', 'aug-cc-pvdz', 'aug-cc-pvdz-jkfit'),
]


class GpuMemorySampler:
    def __init__(self):
        self.peak, self._stop = 0, threading.Event()

    def _read(self):
        r = subprocess.run(['nvidia-smi', '--query-gpu=memory.used', '--format=csv,noheader,nounits'], capture_output=True, text=True)
        return int(r.stdout.split()[0]) if r.returncode == 0 else 0

    def __enter__(self):
        self.peak = self._read()
        def loop():
            while not self._stop.wait(0.25): self.peak = max(self.peak, self._read())
        self._thread = threading.Thread(target=loop, daemon=True); self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set(); self._thread.join()


def synchronize(backend):
    if backend == 'gpu':
        import cupy
        cupy.cuda.Device().synchronize()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--tasks', nargs='+', default=['singlepoint', 'eda'], choices=['singlepoint', 'eda'])
    p.add_argument('--only', nargs='+')
    a = p.parse_args()
    if a.repeats < 2: p.error('at least two repeats required')
    out = Path(a.output); out.mkdir(parents=True, exist_ok=False)
    idle = GpuMemorySampler()._read()
    report = {'hardware': hardware(), 'threads': a.threads, 'repeats': a.repeats, 'idle_gpu_memory_mib': idle,
              'warmup': 'one full call per case/task/backend excluded', 'cases': []}
    for name, xyz, meta_path, xc, basis, aux in CASES:
        if a.only and name not in a.only: continue
        symbols, coords, meta = load(xyz, meta_path)
        for task in a.tasks:
            row = {'system': name, 'task': task, 'natoms': len(symbols), 'backends': {}}
            for backend in ('cpu', 'gpu'):
                settings = Settings(backend=backend, xc=xc, basis=basis, auxbasis=aux, density_fitting=task == 'singlepoint',
                                    grid_level=3, conv_tol=1e-9 if task == 'eda' else 1e-10, conv_tol_grad=1e-6 if task == 'eda' else None,
                                    max_cycle=150, threads=a.threads, compute_duty_cycle=1.0)
                timings, peaks, values = [], [], []
                for rep in range(a.repeats + 1):
                    run = out / f'{name}_{task}_{backend}_{rep}'; run.mkdir()
                    synchronize(backend)
                    with GpuMemorySampler() as mem:
                        start = time.perf_counter()
                        r = singlepoint(symbols, coords, meta, settings) if task == 'singlepoint' else decompose(symbols, coords, meta, settings, run)
                        synchronize(backend)
                        elapsed = time.perf_counter() - start
                    if rep == 0: continue  # warm-up
                    timings.append(elapsed); peaks.append(mem.peak)
                    values.append(r['energy_hartree'] if task == 'singlepoint' else r['counterpoise_interaction_energy_hartree'])
                    print(f'{name} {task} {backend} rep {rep}: {elapsed:.1f} s, peak GPU {mem.peak} MiB', flush=True)
                row['backends'][backend] = {'seconds': timings, 'median_seconds': statistics.median(timings),
                                            'peak_device_memory_mib': max(peaks), 'values_hartree': values, 'settings': asdict(settings)}
            c, g = row['backends']['cpu'], row['backends']['gpu']
            row['cpu_over_gpu_median'] = c['median_seconds'] / g['median_seconds']
            row['cpu_gpu_value_difference_hartree'] = max(abs(x - y) for x, y in zip(c['values_hartree'], g['values_hartree']))
            report['cases'].append(row)
            print(f'== {name} {task}: CPU {c["median_seconds"]:.1f} s, GPU {g["median_seconds"]:.1f} s, ratio {row["cpu_over_gpu_median"]:.2f}', flush=True)
            (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
