#!/usr/bin/env python
"""Fit published IP reagent retention (%ACN at elution) against per-reagent GKS-EDA(sol) descriptors on the C8 patch.

Inputs: a per-placement cases.json (reports/hplc_patch_extended_2026-11-11/check_patch_cases.py) and a retention CSV
(examples/hplc/donegan2022_table4.csv; '#' lines are comments, empty retention = rank-only, most retained).
Per reagent, each descriptor is aggregated over its placements as Boltzmann (298 K, optimized C-PCM electronic energy)
and best (tot, epc, np) or value at the lowest-energy placement (single terms). Placement sets: 'all', and 'off_centre' which drops placements inserted in the patch's central hollow
(lateral r < R_IN and N height h < H_IN), a check of whether the open three-chain pocket drives the result.
Per descriptor: Pearson r, LOO RMSE of a 1-parameter linear fit vs the group-mean RMSE, Spearman over the numeric
reagents and over all (rank-only reagents last). Two-descriptor fits add each EDA descriptor to logD (hydrophobicity
baseline) to test whether EDA explains anything beyond it. Selection check: the retention values are permuted
(N_PERM times) and the whole search repeated; p_select = fraction of permutations whose best LOO RMSE (single, or
logD + X) is at least as good as the observed best. n is small (7) and many descriptors are tested; a good fit
is a hypothesis, not a validation.
"""
import argparse, csv, json, itertools
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

N_PERM = 2000
KT = 0.0019872041 * 298.15; H = 627.509474; R_IN, H_IN = 2.5, 8.0
EDA = {'tot': 'EDA total', 'epc': 'EDA + CDS(water)', 'np': 'buried npSASA (A^2, sign flipped)', 'elec': 'electrostatic',
       'pauli': 'Pauli (exch+rep)', 'pol': 'polarization', 'gc': 'generalized correlation', 'des': 'desolvation'}


def aggregate(cases, model, subset):
    c = {k: v for k, v in cases.items() if k.startswith(f'patch_{model}_') and k.split('_')[-1].isdigit()}
    if subset == 'off_centre': c = {k: v for k, v in c.items() if not (v['r'] < R_IN and v['h'] < H_IN)}
    E0 = min(v['E'] for v in c.values()); w = np.array([np.exp(-(v['E'] - E0) * H / KT) for v in c.values()]); w /= w.sum()
    out = {'n_placements': len(c), 'placements': sorted(c)}
    for q in EDA:
        a = np.array([c[k][q] for k in c]) * (-1 if q == 'np' else 1)
        out[f'{q}/boltzmann'] = float(w @ a)
        if q in ('tot', 'epc', 'np'): out[f'{q}/best'] = float(a.min())
        else: out[f'{q}/lowest_energy'] = float(a[int(np.argmax(w))])  # 'best' of a single term is not a binding measure
    return out


def loo(X, y):
    X = np.column_stack([np.ones(len(y)), X]); err = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i; b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]; err.append(X[i] @ b - y[i])
    return float(np.sqrt(np.mean(np.square(err))))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--cases', required=True); p.add_argument('--data', required=True); p.add_argument('--output', required=True)
    a = p.parse_args()
    cases = json.loads(Path(a.cases).read_text())
    rows = list(csv.DictReader(l for l in open(a.data) if not l.startswith('#')))
    num = [r for r in rows if r['pct_acn_20T']]; rank_only = [r for r in rows if not r['pct_acn_20T']]
    res = {'data': a.data, 'cases': a.cases, 'pocket_cut': dict(r_lt=R_IN, h_lt=H_IN), 'reagents': {}, 'fits': {}}
    for r in rows:
        res['reagents'][r['reagent']] = {s: aggregate(cases, r['model'], s) for s in ('all', 'off_centre')}
    for target in ('pct_acn_20T', 'pct_acn_60T'):
        y = np.array([float(r[target]) for r in num]); base = loo(np.zeros((len(y), 0)), y); fits = {}
        desc = {f'baseline/{k}': ([float(r[k]) for r in num], [float(r[k]) for r in rank_only]) for k in ('logP', 'logD', 'n_carbon')}
        for s in ('all', 'off_centre'):
            for key in res['reagents'][num[0]['reagent']][s]:
                if '/' not in key: continue
                f = lambda rr: [res['reagents'][x['reagent']][s][key] for x in rr]
                desc[f'{s}/{key}'] = ([-v for v in f(num)], [-v for v in f(rank_only)])  # sign: larger = more retained
        for k, (x, xr) in desc.items():
            x = np.array(x); ranks_all = list(range(len(num))) + [len(num)] * len(rank_only)
            order = np.argsort(y); yr = np.empty(len(y)); yr[order] = np.arange(len(y))
            fits[k] = dict(pearson_r=float(np.corrcoef(x, y)[0, 1]), loo_rmse=loo(x[:, None], y), spearman_numeric=float(spearmanr(x, y)[0]),
                           spearman_with_rank_only=float(spearmanr(np.r_[x, xr], np.r_[yr, [len(y)] * len(xr)])[0]))
        logd = np.array(desc['baseline/logD'][0])
        two = {k: loo(np.column_stack([logd, np.array(desc[k][0])]), y) for k in desc if not k.startswith('baseline')}
        X1 = {k: np.array(v[0]) for k, v in desc.items()}; rng = np.random.default_rng(0); null1, null2 = [], []
        for _ in range(N_PERM):
            yp = rng.permutation(y)
            null1.append(min(loo(x[:, None], yp) for x in X1.values()))
            null2.append(min(loo(np.column_stack([logd, x]), yp) for k, x in X1.items() if not k.startswith('baseline')))
        b1, b2 = min(v['loo_rmse'] for v in fits.values()), min(two.values())
        sel = dict(n_perm=N_PERM, best_single=b1, p_select_single=float(np.mean(np.array(null1) <= b1)),
                   best_logD_plus=b2, p_select_logD_plus=float(np.mean(np.array(null2) <= b2)),
                   p_logD_alone=float(np.mean([loo(logd[:, None], rng.permutation(y)) <= fits['baseline/logD']['loo_rmse'] for _ in range(N_PERM)])))
        res['fits'][target] = dict(selection_check=sel, n=len(y), loo_rmse_group_mean=base, single=dict(sorted(fits.items(), key=lambda t: t[1]['loo_rmse'])),
                                   logD_plus=dict(sorted(two.items(), key=lambda t: t[1])))
    Path(a.output).write_text(json.dumps(res, indent=1) + '\n')
    for t, f in res['fits'].items():
        print(f"== {t}: n {f['n']}, LOO RMSE of group mean {f['loo_rmse_group_mean']:.2f} %ACN")
        for k, v in f['single'].items():
            print(f"  {k:32s} r {v['pearson_r']:+.2f}  LOO {v['loo_rmse']:5.2f}  rho {v['spearman_numeric']:+.2f}  rho(+TBA) {v['spearman_with_rank_only']:+.2f}")
        print('  selection check:', {k: round(v, 3) for k, v in f['selection_check'].items()})
        print('  logD + X (LOO):', ', '.join(f"{k} {v:.2f}" for k, v in list(f['logD_plus'].items())[:6]))
    for g, d in res['reagents'].items(): print(g, 'off-centre placements', d['off_centre']['n_placements'], 'of', d['all']['n_placements'])


if __name__ == '__main__':
    main()
