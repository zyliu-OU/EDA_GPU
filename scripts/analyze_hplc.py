#!/usr/bin/env python
"""Tabulate an HPLC campaign: per complex group, the lowest optimized PCM free energy placement
and the spread of GKS-EDA(sol) components across all successful placements (kcal/mol)."""
import argparse
import json
import re
from pathlib import Path

KCAL = 627.509474
ORDER = ['electrostatic', 'exchange', 'repulsion', 'polarization', 'generalized_correlation', 'desolvation', 'total']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign'); p.add_argument('--output', required=True)
    a = p.parse_args()
    summary = json.loads((Path(a.campaign) / 'summary.json').read_text())
    groups = {}
    for name, case in summary['cases'].items():
        groups.setdefault(re.sub(r'_\d+_\d+$', '', name), []).append((name, case))
    table = {}
    for group, cases in sorted(groups.items()):
        done = [(n, c) for n, c in cases if c.get('eda') == 'success']
        row = {'placements': len(cases), 'successful': len(done), 'failed': {n: c.get('error') for n, c in cases if c.get('error')}}
        if done:
            best_name, best = min(done, key=lambda t: t[1]['optimized_energy_hartree'])
            spread = {k: [c['components_kcal_mol'][k] for _, c in done] for k in ORDER}
            row.update(selected=best_name, selected_relative_free_energy_kcal_mol={n: (c['optimized_energy_hartree'] - best['optimized_energy_hartree']) * KCAL for n, c in done},
                       components_kcal_mol={k: best['components_kcal_mol'][k] for k in ORDER},
                       component_range_kcal_mol={k: [min(v), max(v)] for k, v in spread.items()})
        table[group] = row
    out = {'campaign': a.campaign, 'optimization_settings': summary.get('optimization_settings'), 'eda_settings': summary.get('eda_settings'),
           'selection': 'lowest optimized PCM free energy among successful placements; ranges over all successful placements',
           'groups': table}
    Path(a.output).parent.mkdir(parents=True, exist_ok=True); Path(a.output).write_text(json.dumps(out, indent=2) + '\n')
    print('| Group | n | ' + ' | '.join(ORDER) + ' |')
    print('| --- | ---: | ' + ' | '.join('---:' for _ in ORDER) + ' |')
    for g, r in table.items():
        if 'components_kcal_mol' in r:
            print(f"| {g} | {r['successful']}/{r['placements']} | " + ' | '.join(f"{r['components_kcal_mol'][k]:.2f}" for k in ORDER) + ' |')
        else:
            print(f"| {g} | 0/{r['placements']} | " + ' | '.join('—' for _ in ORDER) + ' |')


if __name__ == '__main__':
    main()
