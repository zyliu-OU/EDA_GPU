#!/usr/bin/env python
"""Three-body GKS-EDA(sol) analysis of an optimized residue (R) / TEAH+ (T) / proxy (S) complex.

At the common ternary geometry: (RT)|S and R|(TS) partitions of the full complex, and the pairs
R|T, T|S, R|S extracted from it. Since E(RTS) - E(RT) - E(S) = dE((RT)|S),
the three-body term is dE3 = dE((RT)|S) - dE(T|S) - dE(R|S), evaluated component by component.
Counterpoise caveat: (RT)|S uses the trimer basis, each pair its own dimer basis.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
from eda_gpu.engine import Settings, hardware
from eda_gpu.gks import decompose
from eda_gpu.model import load

KCAL = 627.509474


def subsystem(symbols, coords, meta, names):
    """Complex made of the named fragments, re-indexed, with merged fragments where names are tuples."""
    frags = {f['name']: f for f in meta['fragments']}
    groups = [n if isinstance(n, tuple) else (n,) for n in names]
    atoms = [a for g in groups for n in g for a in frags[n]['atoms']]
    index = {a: i for i, a in enumerate(atoms)}
    fragments = []
    for g in groups:
        fragments.append({'name': '+'.join(g), 'atoms': [index[a] for n in g for a in frags[n]['atoms']],
                          'charge': sum(frags[n]['charge'] for n in g), 'spin': sum(frags[n]['spin'] for n in g)})
    sub = {'identity': f"{meta['identity']} [subsystem {' | '.join(f['name'] for f in fragments)}]", 'protonation': meta['protonation'],
           'charge': sum(f['charge'] for f in fragments), 'spin': sum(f['spin'] for f in fragments), 'fragments': fragments,
           'bonds': [[index[i], index[j], o] for i, j, o in meta['bonds'] if i in index and j in index],
           'stereochemistry': meta['stereochemistry']}
    return [symbols[a] for a in atoms], [coords[a] for a in atoms], sub


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('xyz'); p.add_argument('--metadata', required=True); p.add_argument('--output', required=True)
    p.add_argument('--names', nargs=3, required=True, metavar=('R', 'T', 'S'), help='fragment names for residue, reagent, surface')
    p.add_argument('--basis', default='def2-svp'); p.add_argument('--xc', default='m06-2x')
    p.add_argument('--solvent', default='cpcm'); p.add_argument('--solvent-eps', type=float, default=78.3553)
    a = p.parse_args()
    symbols, coords, meta = load(a.xyz, a.metadata)
    r, t, s = a.names
    settings = Settings(backend='gpu', xc=a.xc, basis=a.basis, density_fitting=False, grid_level=3, conv_tol=1e-9, conv_tol_grad=1e-6,
                        max_cycle=200, threads=16, compute_duty_cycle=1.0, gpu_temperature_limit=80.0, cpu_temperature_limit=85.0, solvent=a.solvent, solvent_eps=a.solvent_eps)
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    jobs = {'RT|S': [(r, t), s], 'R|TS': [r, (t, s)], 'R|T': [r, t], 'T|S': [t, s], 'R|S': [r, s]}
    report_path = out / 'three_body.json'
    report = json.loads(report_path.read_text()) if report_path.exists() else {'jobs': {}}
    report.update(hardware=hardware(), settings=asdict(settings), source_geometry=a.xyz,
                  definition='dE3 = dE((RT)|S) - dE(T|S) - dE(R|S) at the ternary geometry; pairs in their own dimer basis')
    for label, names in jobs.items():
        if report['jobs'].get(label, {}).get('status') == 'success':
            continue
        ss, sc, sm = subsystem(symbols, coords, meta, names)
        d = out / label.replace('|', '_'); d.mkdir(exist_ok=True)
        print(f'{label}: {len(ss)} atoms', flush=True)
        try:
            res = decompose(ss, sc, sm, settings, d)
            (d / 'result.json').write_text(json.dumps(res, indent=2, allow_nan=False) + '\n')
            comps = dict(res['components_hartree']); comps['total'] = res['counterpoise_interaction_energy_hartree']
            report['jobs'][label] = {'status': 'success', 'components_kcal_mol': {k: v * KCAL for k, v in comps.items()}}
        except Exception as exc:
            report['jobs'][label] = {'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    j = report['jobs']
    if all(j.get(k, {}).get('status') == 'success' for k in ('RT|S', 'T|S', 'R|S')):
        report['three_body_kcal_mol'] = {k: j['RT|S']['components_kcal_mol'][k] - j['T|S']['components_kcal_mol'][k] - j['R|S']['components_kcal_mol'][k]
                                         for k in j['RT|S']['components_kcal_mol']}
        print('three-body: ' + '  '.join(f'{k}={v:+.2f}' for k, v in report['three_body_kcal_mol'].items()), flush=True)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
