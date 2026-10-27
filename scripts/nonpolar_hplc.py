#!/usr/bin/env python
"""SMD CDS (non-electrostatic) association terms of optimized HPLC complexes, beside GKS-EDA(sol).

C-PCM GKS-EDA(sol) has no cavitation/dispersion solvation term. For every converged, topology-preserving
case, this script computes dG_CDS = G_CDS(AB) - G_CDS(A) - G_CDS(B) at the complex geometry (SMD,
water and acetonitrile; eda_gpu.cds) and the buried nonpolar SASA (scripts/buried_sasa.py), for each
fragment pair and, for three-fragment complexes, the R|T, T|S, R|S, RT|S, R|TS splits (fragment order
R, T, S). EDA totals come from the campaign summary (pairs) or from --three-body directories holding
<case>/three_body.json. EDA + dG_CDS(water) is a mixed-model estimate (SMD CDS was fitted with SMD
electrostatics and radii, not C-PCM), reported only as a sensitivity check of the hydrophobic term.
"""
import argparse
import json
import re
from pathlib import Path
import numpy as np
from eda_gpu.cds import SOLVENTS, association_cds
from eda_gpu.model import load
from buried_sasa import buried


def splits(names):
    if len(names) == 2:
        return {'|'.join(names): [names[:1], names[1:]]}
    r, t, s = names
    return {'R|T': [[r], [t]], 'T|S': [[t], [s]], 'R|S': [[r], [s]], 'RT|S': [[r, t], [s]], 'R|TS': [[r], [t, s]]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaigns', nargs='+', help='run directories containing summary.json')
    p.add_argument('--three-body', action='append', default=[], help='directory of <case>/three_body.json (repeatable)')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    result = {'method': 'SMD CDS (Marenich, Cramer, Truhlar 2009) via PySCF/MNSOL at the complex geometry; '
                        'buried SASA as in scripts/buried_sasa.py; kcal/mol and A^2', 'cases': {}}
    for campaign in a.campaigns:
        summary = json.loads(Path(campaign, 'summary.json').read_text())
        for name, case in summary['cases'].items():
            if case.get('optimization') != 'converged' or not case.get('bonding_topology_preserved'):
                continue
            symbols, coords, meta = load(Path(campaign, name, 'optimization', 'optimized.xyz'), case['seed'] + '.json')
            coords = np.asarray(coords, float)
            frags = {f['name']: f['atoms'] for f in meta['fragments']}
            eda = {}
            if case.get('eda') == 'success':
                eda['|'.join(frags)] = case['components_kcal_mol']['total']
            for d in a.three_body:
                tb = Path(d, name, 'three_body.json')
                if tb.exists():
                    eda.update({k: j['components_kcal_mol']['total'] for k, j in json.loads(tb.read_text())['jobs'].items()
                                if j.get('status') == 'success'})
            rows = {}
            for label, groups in splits(list(frags)).items():
                idx = [[i for n in g for i in frags[n]] for g in groups]
                row = {f'dG_cds_{s}': association_cds(symbols, coords, idx, s) for s in SOLVENTS}
                row.update(buried(symbols, coords, meta, groups))
                if label in eda:
                    row['eda_total'] = eda[label]
                    row['eda_plus_cds_water'] = eda[label] + row['dG_cds_water']
                rows[label] = row
            if len(frags) == 3:
                rows['three_body'] = {k: rows['RT|S'][k] - rows['T|S'][k] - rows['R|S'][k]
                                      for k in rows['RT|S'] if all(k in rows[x] for x in ('T|S', 'R|S'))}
            key = f'{Path(campaign).name}/{name}'
            result['cases'][key] = dict(group=re.sub(r'(_contact)?_\d+_\d+$', '', name),
                                        optimized_energy_hartree=case['optimized_energy_hartree'], splits=rows)
            print(key, '  '.join(f"{k}: cds_w {v['dG_cds_water']:+.2f} cds_acn {v['dG_cds_acetonitrile']:+.2f} "
                                 f"npSASA {v['buried_nonpolar_A2']:.0f}" + (f" eda {v['eda_total']:+.2f}" if 'eda_total' in v else '')
                                 for k, v in rows.items()), flush=True)
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
