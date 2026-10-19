#!/usr/bin/env python
"""Resumable optimize -> GKS-EDA(sol) campaign over labelled complex seeds.

For each <seed>.xyz/.json: GPU optimization (density fitting, PCM, frozen atoms from
metadata), then direct-integral counterpoise GKS-EDA(sol) at the optimized geometry with the
same functional, basis and solvent. Completed stages are skipped on restart; a failure is
recorded and the campaign continues. Optimized structures are stationary points, not
established minima, and EDA values are electronic descriptors, not retention free energies.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import traceback
from eda_gpu.engine import Settings, hardware, optimize
from eda_gpu.gks import decompose
from eda_gpu.model import load

KCAL = 627.509474


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('seeds', nargs='+', help='seed stems (without .xyz/.json)')
    p.add_argument('--output', required=True)
    p.add_argument('--xc', default='m06-2x'); p.add_argument('--basis', default='def2-svp')
    p.add_argument('--auxbasis', default='def2-universal-jkfit')
    p.add_argument('--solvent', default='cpcm', choices=['none', 'cpcm', 'iefpcm']); p.add_argument('--solvent-eps', type=float, default=78.3553)
    p.add_argument('--maxsteps', type=int, default=200); p.add_argument('--threads', type=int, default=16)
    p.add_argument('--skip-eda', action='store_true')
    a = p.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    common = dict(backend='gpu', xc=a.xc, basis=a.basis, auxbasis=a.auxbasis, grid_level=3, max_cycle=200, threads=a.threads,
                  compute_duty_cycle=1.0, gpu_temperature_limit=80.0, solvent=a.solvent, solvent_eps=a.solvent_eps)
    opt_settings = Settings(density_fitting=True, conv_tol=1e-10, **common)
    eda_settings = Settings(density_fitting=False, conv_tol=1e-9, conv_tol_grad=1e-6, **common)
    summary_path = out / 'summary.json'
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {'cases': {}}
    summary.update(hardware=hardware(), optimization_settings=asdict(opt_settings), eda_settings=asdict(eda_settings))
    for stem in a.seeds:
        name = Path(stem).name
        case = summary['cases'].setdefault(name, {'seed': stem})
        symbols, coords, meta = load(stem + '.xyz', stem + '.json')
        opt_dir, eda_dir = out / name / 'optimization', out / name / 'eda'
        try:
            if case.get('optimization') != 'converged':
                if opt_dir.exists(): opt_dir.rename(opt_dir.with_name(f'optimization_failed_{len(list(opt_dir.parent.glob("optimization_failed_*")))}'))
                opt_dir.mkdir(parents=True)
                print(f'{name}: optimizing ({len(symbols)} atoms)', flush=True)
                r = optimize(symbols, coords, meta, opt_settings, opt_dir, maxsteps=a.maxsteps)
                write(opt_dir / 'result.json', r)
                case.update(optimization='converged' if r['optimization_converged'] else 'not converged',
                            optimized_energy_hartree=r['energy_hartree'], evaluations=len(r['evaluations']),
                            frozen_max_displacement_angstrom=r['frozen_max_displacement_angstrom'])
                write(summary_path, summary)
                if not r['optimization_converged']: raise RuntimeError('optimization not converged')
            if not a.skip_eda and case.get('eda') != 'success':
                _, xyz, _ = load(opt_dir / 'optimized.xyz', stem + '.json')
                if eda_dir.exists(): eda_dir.rename(eda_dir.with_name(f'eda_failed_{len(list(eda_dir.parent.glob("eda_failed_*")))}'))
                eda_dir.mkdir(parents=True)
                print(f'{name}: GKS-EDA(sol)', flush=True)
                r = decompose(symbols, xyz, meta, eda_settings, eda_dir)
                write(eda_dir / 'result.json', r)
                comps = dict(r['components_hartree']); comps['total'] = r['counterpoise_interaction_energy_hartree']
                case.update(eda='success', components_kcal_mol={k: v * KCAL for k, v in comps.items()},
                            closure_residual_hartree=r['total_closure_residual_hartree'])
                print(f'{name}: ' + '  '.join(f'{k}={v:+.2f}' for k, v in case['components_kcal_mol'].items()), flush=True)
        except Exception as exc:
            case['error'] = f'{type(exc).__name__}: {exc}'; case['traceback'] = traceback.format_exc()
            print(f'{name}: FAILED {case["error"]}', flush=True)
        else:
            case.pop('error', None); case.pop('traceback', None)
        write(summary_path, summary)


if __name__ == '__main__':
    main()
