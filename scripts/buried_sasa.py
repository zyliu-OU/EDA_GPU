#!/usr/bin/env python
"""Buried solvent-accessible surface area of optimized HPLC complexes (empirical hydrophobic descriptor).

C-PCM GKS-EDA(sol) has no cavitation term, so the hydrophobic effect that drives alkyl-chain adsorption
is absent from every EDA component. This script reports, separately from the EDA, the SASA buried on
association at the optimized complex geometry: dSASA(A|B) = SASA(A) + SASA(B) - SASA(AB), total and
nonpolar (C, and H bonded to C or Si). Shrake-Rupley, 1.4 A probe, Bondi radii (Si 2.10, H 1.20),
960 golden-spiral points per atom. Areas only; a free energy needs a surface-tension coefficient,
whose literature values span ~0.005-0.03 kcal/mol/A^2, so none is applied here.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from eda_gpu.model import load

RADII = {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80, 'Si': 2.10}
PROBE, POINTS = 1.4, 960


def sphere(n):
    k = np.arange(n) + 0.5
    phi, theta = np.arccos(1 - 2 * k / n), np.pi * (1 + 5 ** 0.5) * k
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


def sasa(symbols, coords):
    """Per-atom SASA (A^2)."""
    r = np.array([RADII[s] for s in symbols]) + PROBE
    unit, out = sphere(POINTS), np.zeros(len(symbols))
    for i in range(len(symbols)):
        near = np.where(np.linalg.norm(coords - coords[i], axis=1) < r + r[i])[0]
        near = near[near != i]
        pts = coords[i] + r[i] * unit
        free = np.all(np.linalg.norm(pts[:, None] - coords[near][None], axis=2) >= r[near], axis=1)
        out[i] = 4 * np.pi * r[i] ** 2 * free.mean()
    return out


def nonpolar(symbols, bonds):
    partner = {}
    for i, j, _ in bonds:
        partner.setdefault(i, []).append(j); partner.setdefault(j, []).append(i)
    return np.array([s == 'C' or (s == 'H' and any(symbols[k] in ('C', 'Si') for k in partner.get(i, []))) for i, s in enumerate(symbols)])


def buried(symbols, coords, meta, groups):
    """dSASA between two groups of fragment names, both at the complex geometry."""
    frags = {f['name']: f['atoms'] for f in meta['fragments']}
    a = [i for n in groups[0] for i in frags[n]]; b = [i for n in groups[1] for i in frags[n]]
    np_mask = nonpolar(symbols, meta['bonds'])
    whole = sasa([symbols[i] for i in a + b], coords[a + b])
    parts = np.r_[sasa([symbols[i] for i in a], coords[a]), sasa([symbols[i] for i in b], coords[b])]
    m = np_mask[a + b]
    return {'buried_total_A2': float(parts.sum() - whole.sum()), 'buried_nonpolar_A2': float(parts[m].sum() - whole[m].sum())}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaigns', nargs='+', help='run directories containing summary.json')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    result = {'method': 'Shrake-Rupley, probe 1.4 A, Bondi radii, 960 points; nonpolar = C and H on C/Si', 'cases': {}}
    for campaign in a.campaigns:
        for name, case in json.loads(Path(campaign, 'summary.json').read_text())['cases'].items():
            xyz = Path(campaign, name, 'optimization', 'optimized.xyz')
            if case.get('optimization') != 'converged' or not case.get('bonding_topology_preserved'): continue
            symbols, coords, meta = load(xyz, case['seed'] + '.json')
            coords, names = np.asarray(coords, float), [f['name'] for f in meta['fragments']]
            splits = {f'{x}|{y}': [[x], [y]] for i, x in enumerate(names) for y in names[i + 1:]}
            if len(names) == 3:
                splits[f'{names[0]}+{names[1]}|{names[2]}'] = [names[:2], names[2:]]
            result['cases'][f'{Path(campaign).name}/{name}'] = {k: buried(symbols, coords, meta, g) for k, g in splits.items()}
            print(name, '  '.join(f"{k}: {v['buried_nonpolar_A2']:.0f}/{v['buried_total_A2']:.0f}" for k, v in result['cases'][f'{Path(campaign).name}/{name}'].items()), flush=True)
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
