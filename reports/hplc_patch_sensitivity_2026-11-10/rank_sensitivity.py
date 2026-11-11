"""Per-variant reagent ranking for the solvent sensitivity single points.
Usage (repo root): python rank_sensitivity.py <sensitivity summary.json> <sampling cases.json> <out ranking.json>
Variant rows: tot = variant GKS-EDA(sol) total; E = variant complex energy (Boltzmann weights, 298 K, within reagent);
epc = tot + SMD CDS of the matching solvent (water for base/vdw14, acetonitrile for acn; CDS is geometry-only).
'base' = the original 10-27/10-28 values from cases.json."""
import json, sys, numpy as np
from scipy.stats import spearmanr
sens = json.load(open(sys.argv[1]))['cases']; base = json.load(open(sys.argv[2]))
pub = ['teah', 'dipeah', 'dibutylammonium', 'hexylammonium', 'tributylammonium']  # least -> most retained
kT = 0.0019872041 * 298.15; H = 627.509474
tables = {'base': {k: dict(E=v['E'], tot=v['tot'], epc=v['epc']) for k, v in base.items()}}
for var, cds in [('vdw14', 'cdsw'), ('acn', 'cdsm')]:
    tables[var] = {k: dict(E=c[var]['complex_energy_hartree'], tot=c[var]['components_kcal_mol']['total'],
                           epc=c[var]['components_kcal_mol']['total'] + base[k][cds]) for k, c in sens.items()}
out = dict(published_order_least_to_most_retained=pub, variants={})
for var, rows in tables.items():
    R = {}
    for g in pub:
        c = {k: v for k, v in rows.items() if k.startswith(f'patch_{g}_')}
        E0 = min(v['E'] for v in c.values())
        w = np.array([np.exp(-(v['E'] - E0) * H / kT) for v in c.values()]); w /= w.sum()
        low = min(c, key=lambda k: c[k]['E'])
        d = dict(lowest_energy_placement=low, weights={k: round(float(x), 4) for k, x in zip(c, w)})
        for q in ['tot', 'epc']:
            a = np.array([c[k][q] for k in c])
            d[q] = dict(lowest_energy=c[low][q], best=float(a.min()), boltzmann=float(w @ a), mean=float(a.mean()))
        R[g] = d
    sp = {}
    for q in ['tot', 'epc']:
        for m in ['lowest_energy', 'best', 'boltzmann', 'mean']:
            x = [-R[g][q][m] for g in pub]; r, p = spearmanr(range(5), x)
            sp[f'{q}/{m}'] = dict(rho=round(float(r), 2), p=round(float(p), 2), order=[pub[i] for i in np.argsort(x)])
    out['variants'][var] = dict(reagents=R, spearman=sp)
# per-placement shifts vs base
out['shifts_vs_base'] = {var: {k: dict(dtot=round(tables[var][k]['tot'] - tables['base'][k]['tot'], 3)) for k in sorted(sens)}
                         for var in ['vdw14', 'acn']}
for var in ['vdw14', 'acn']:
    s = np.array([v['dtot'] for v in out['shifts_vs_base'][var].values()])
    xb = [tables['base'][k]['tot'] for k in sorted(sens)]; xv = [tables[var][k]['tot'] for k in sorted(sens)]
    out['shifts_vs_base'][var + '_stats'] = dict(mean=round(float(s.mean()), 3), sd=round(float(s.std()), 3),
        min=round(float(s.min()), 3), max=round(float(s.max()), 3), rho_placements=round(float(spearmanr(xb, xv)[0]), 3))
json.dump(out, open(sys.argv[3], 'w'), indent=1)
for var, o in out['variants'].items():
    print(f'== {var}')
    for g in pub:
        d = o['reagents'][g]; print(f'  {g:18s} low={d["lowest_energy_placement"][-1]} ' + '  '.join(
            f'{q}:{m[:4]}={d[q][m]:+.2f}' for q in ['tot', 'epc'] for m in ['lowest_energy', 'best', 'boltzmann', 'mean']))
    for k, v in o['spearman'].items(): print(f'  {k:20s} rho={v["rho"]:+.2f} p={v["p"]:.2f} {" < ".join(x[:5] for x in v["order"])}')
for var in ['vdw14', 'acn']: print(var, out['shifts_vs_base'][var + '_stats'])
