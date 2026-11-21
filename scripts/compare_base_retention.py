#!/usr/bin/env python
"""Compare per-base GKS-EDA(sol) descriptors from a patch campaign (capped RNA residues on the C8 patch) with base retention data.

Per base, each descriptor is aggregated over its successful placements as: best (value at the lowest optimized C-PCM
energy placement), Boltzmann (298 K, optimized C-PCM electronic energy, as in fit_patch_retention.py), min (most bound)
and mean. Descriptors: GKS-EDA(sol) components and total, EDA + SMD-CDS(water), buried nonpolar SASA, and how many
patch-carbon contacts (< 4.5 A, heavy atoms) come from the base ring vs the sugar-phosphate backbone.
Retention: a CSV with base, model, retention_time_min (more retained = larger). With 4 bases, Spearman rho = +/-1 is the
only result with exact one-sided p < 0.05 (1/24). Sign convention: stronger binding (more negative energy) should mean
more retained, so the expected rho of an energy descriptor with retention time is negative.
"""
import argparse, csv, json, itertools
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

KT = 0.0019872041 * 298.15; H = 627.509474
COMP = ['electrostatic', 'exchange', 'repulsion', 'polarization', 'generalized_correlation', 'desolvation', 'total']


def read_xyz(p):
    L = Path(p).read_text().split('\n'); n = int(L[0])
    s = [l.split()[0] for l in L[2:2 + n]]; X = np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + n]])
    return s, X


def contacts(seed_json, xyz):
    """Patch-carbon contacts (< 4.5 A) of guest base-ring heavy atoms vs backbone (sugar + phosphate + caps) heavy atoms."""
    m = json.loads(Path(seed_json).read_text())
    frags = m['fragments'] if isinstance(m['fragments'], list) else eval(m['fragments'])
    bonds = m['bonds'] if isinstance(m['bonds'], list) else eval(m['bonds'])
    s, X = read_xyz(xyz)
    host = set(frags[0]['atoms']); guest = [i for i in range(len(s)) if i not in host]
    nb = {i: set() for i in range(len(s))}
    for i, j, *_ in bonds: nb[i].add(j); nb[j].add(i)
    gly = None  # glycosidic N-C1': N bonded to a C whose neighbour is the ring ether O4' (O with two C neighbours)
    for n in guest:
        if s[n] != 'N': continue
        for c in nb[n]:
            if s[c] == 'C' and any(s[o] == 'O' and sum(s[k] == 'C' for k in nb[o]) == 2 for o in nb[c]): gly = (n, c)
    base, stack = set(), [gly[0]]
    while stack:
        a = stack.pop()
        if a in base: continue
        base.add(a); stack += [b for b in nb[a] if b != gly[1]]
    pc = [i for i in host if s[i] == 'C']
    r = np.linalg.norm(X[guest][:, None] - X[pc][None], axis=-1)
    heavy = [k for k, i in enumerate(guest) if s[i] != 'H']
    cnt = {i: int((r[k] < 4.5).sum()) for k, i in enumerate(guest) if k in heavy}
    return {'base_contacts': sum(v for i, v in cnt.items() if i in base), 'backbone_contacts': sum(v for i, v in cnt.items() if i not in base)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign'); p.add_argument('--data', required=True); p.add_argument('--seeds', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    camp = Path(a.campaign)
    cases = json.loads((camp / 'summary.json').read_text())['cases']
    npj = json.loads((camp / 'nonpolar.json').read_text())['cases']
    rows = [r for r in csv.DictReader(l for l in open(a.data) if not l.startswith('#'))]
    per = {}
    for r in rows:
        model = r['model']
        names = sorted(k for k, v in cases.items() if k.startswith(f'patch_{model}_') and k.split('_')[-1].isdigit() and v.get('eda') == 'success')
        pl = {}
        for n in names:
            sp = next(iter(next(v for k, v in npj.items() if k.endswith('/' + n))['splits'].values()))
            d = {k: cases[n]['components_kcal_mol'][k] for k in COMP}
            d.update(pauli=d['exchange'] + d['repulsion'], eda_plus_cds=sp['eda_plus_cds_water'], cds_water=sp['dG_cds_water'], buried_np_A2=sp['buried_nonpolar_A2'])
            d.update(contacts(Path(a.seeds) / f'{n}.json', camp / n / 'optimization' / 'optimized.xyz'))
            d['E'] = cases[n]['optimized_energy_hartree']; pl[n] = d
        E0 = min(v['E'] for v in pl.values())
        w = np.array([np.exp(-(v['E'] - E0) * H / KT) for v in pl.values()]); w /= w.sum()
        keys = [k for k in next(iter(pl.values())) if k != 'E']
        best = min(pl, key=lambda n: pl[n]['E'])
        agg = {}
        for k in keys:
            x = np.array([v[k] for v in pl.values()])
            agg[f'{k}/best'] = pl[best][k]; agg[f'{k}/boltzmann'] = float(w @ x); agg[f'{k}/min'] = float(x.min()); agg[f'{k}/mean'] = float(x.mean())
        per[r['base']] = dict(model=model, retention_time_min=float(r['retention_time_min']), n_placements=len(pl), selected=best,
                              weights={n: float(x) for n, x in zip(pl, w)}, relative_energy_kcal_mol={n: (v['E'] - E0) * H for n, v in pl.items()},
                              placements={n: {k: v for k, v in d.items() if k != 'E'} for n, d in pl.items()}, descriptors=agg)
    bases = list(per); t = np.array([per[b]['retention_time_min'] for b in bases])
    # exact permutation null for Spearman with n = 4
    nperm = list(itertools.permutations(range(len(bases))))
    corr = {}
    for k in per[bases[0]]['descriptors']:
        x = np.array([per[b]['descriptors'][k] for b in bases])
        rho = spearmanr(x, t)[0]; rr = pearsonr(x, t)[0] if np.ptp(x) > 0 else float('nan')
        p_neg = sum(spearmanr(x[list(q)], t)[0] <= rho + 1e-12 for q in nperm) / len(nperm)
        corr[k] = dict(spearman_rho=float(rho), pearson_r=float(rr), p_one_sided_negative=p_neg,
                       order_weak_to_strong=''.join(b for _, b in sorted(zip(x, bases), reverse=True)))
    out = dict(campaign=str(camp), data=a.data, retention_order_low_to_high=''.join(b for _, b in sorted(zip(t, bases))),
               note='expected sign for energies: negative rho (stronger binding = more retained); for buried SASA / contacts: positive',
               bases=per, correlations=corr)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True); Path(a.output).write_text(json.dumps(out, indent=1) + '\n')
    print('retention order (least -> most retained):', out['retention_order_low_to_high'])
    print('| descriptor | ' + ' | '.join(bases) + ' | rho | r | weak->strong |')
    print('| --- | ' + ' | '.join('---:' for _ in bases) + ' | ---: | ---: | --- |')
    for k, c in corr.items():
        if k.split('/')[1] in ('best', 'boltzmann'):
            print(f'| {k} | ' + ' | '.join(f"{per[b]['descriptors'][k]:.2f}" for b in bases) + f" | {c['spearman_rho']:+.2f} | {c['pearson_r']:+.2f} | {c['order_weak_to_strong']} |")


if __name__ == '__main__':
    main()
