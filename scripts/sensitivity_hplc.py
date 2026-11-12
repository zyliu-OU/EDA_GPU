#!/usr/bin/env python
"""Solvent-model sensitivity: GKS-EDA(sol) single points at existing optimized geometries.

For every case directory (<run>/<case>/optimization/optimized.xyz plus the seed metadata named in
<run>/summary.json) the EDA is repeated with the run's own recorded eda_settings, changing only the
fields of each variant (e.g. solvent_vdw_scale=1.4, solvent_eps=35.688, or basis aug-cc-pVDZ). Geometries are NOT
re-optimized: this isolates the response of the EDA terms to the solvent model at fixed structure.
Resumable; GPU/CPU temperature guards come from the recorded settings.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import traceback
from eda_gpu.engine import Settings
from eda_gpu.gks import decompose
from eda_gpu.model import load

KCAL = 627.509474
VARIANTS = {'vdw14': dict(solvent_vdw_scale=1.4), 'acn': dict(solvent_eps=35.688), 'ref': {},
            'augdz': dict(basis='aug-cc-pvdz', auxbasis='aug-cc-pvdz-jkfit')}  # EDA is direct-integral; auxbasis unused


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs', nargs='+', help='run directories with summary.json and <case>/optimization/optimized.xyz')
    p.add_argument('--variants', nargs='+', default=['vdw14', 'acn'], choices=sorted(VARIANTS))
    p.add_argument('--cases', nargs='*', help='only these case names (default: every converged case with an EDA)')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for run in map(Path, a.runs):
        s = json.loads((run / 'summary.json').read_text())
        base = Settings(**s['eda_settings'])
        for name, c in s['cases'].items():
            if c.get('eda') == 'success' and (not a.cases or name in a.cases):
                jobs.append((name, run, c['seed'], base))
    jobs.sort(key=lambda j: j[0])
    summary_path = out / 'summary.json'
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {'variants': {}, 'cases': {}}
    summary['variants'].update({v: VARIANTS[v] for v in a.variants})
    for v in a.variants:  # variant-major: a complete variant is usable before the next starts
        for name, run, seed, base in jobs:
            case = summary['cases'].setdefault(name, {'run': str(run), 'seed': seed})
            if case.get(v, {}).get('eda') == 'success':
                continue
            settings = replace(base, **VARIANTS[v])
            d = out / v / name
            try:
                symbols, xyz, meta = load(run / name / 'optimization' / 'optimized.xyz', seed + '.json')
                d.mkdir(parents=True, exist_ok=True)
                print(f'{v} {name}: GKS-EDA(sol) ({len(symbols)} atoms)', flush=True)
                r = decompose(symbols, xyz, meta, settings, d)
                write(d / 'result.json', r)
                comps = dict(r['components_hartree']); comps['total'] = r['counterpoise_interaction_energy_hartree']
                case[v] = dict(eda='success', components_kcal_mol={k: x * KCAL for k, x in comps.items()},
                               complex_energy_hartree=r['reference_states']['complex']['energy_hartree'],
                               closure_residual_hartree=r['total_closure_residual_hartree'], seconds=r['seconds'],
                               solvent=r['solvent'])
                print(f'{v} {name}: ' + '  '.join(f'{k}={x:+.2f}' for k, x in case[v]['components_kcal_mol'].items()), flush=True)
            except Exception as exc:
                case[v] = dict(eda='failed', error=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc())
                print(f'{v} {name}: FAILED {case[v]["error"]}', flush=True)
            write(summary_path, summary)


if __name__ == '__main__':
    main()
