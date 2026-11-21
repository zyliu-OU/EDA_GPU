#!/usr/bin/env python
"""Per-base retention contributions from published oligonucleotide data (Sturm et al. 2007 Table S1).

Fits, at each temperature, the Gilar 2002 form RT = sum_i f_i (a_i ln N + e_i) (f_i = base fraction) and the additive form
RT = b0 + b1 ln N + sum_{i != C} d_i f_i (C reference), with leave-one-out RMSE. Also averages the RT shift of all
same-length pairs that differ by exactly one base exchange (canonical direction X -> Y with X before Y in 'CGAT')."""
import argparse, csv, json, itertools
import numpy as np


def loo(X, y):
    e = [y[i] - X[i] @ np.linalg.lstsq(X[np.arange(len(y)) != i], y[np.arange(len(y)) != i], rcond=None)[0] for i in range(len(y))]
    return float(np.sqrt(np.mean(np.square(e))))


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--data', required=True); p.add_argument('--output', required=True)
    a = p.parse_args()
    R = [r for r in csv.DictReader(l for l in open(a.data) if not l.startswith('#'))]
    N = np.array([int(r['length']) for r in R]); n = {b: np.array([int(r['n' + b]) for r in R]) for b in 'ATCG'}
    temps = [k[3:] for k in R[0] if k.startswith('rt_')]
    out = {'data': a.data, 'n_oligos': len(R), 'fits': {}, 'single_exchange_pairs': {}}
    for T in temps:
        y = np.array([float(r['rt_' + T]) for r in R])
        X = np.column_stack([n[b] / N * np.log(N) for b in 'ATCG'] + [n[b] / N for b in 'ATCG'])
        c = np.linalg.lstsq(X, y, rcond=None)[0]; r2 = 1 - np.sum((y - X @ c) ** 2) / np.sum((y - y.mean()) ** 2)
        X2 = np.column_stack([np.ones(len(y)), np.log(N)] + [n[b] / N for b in 'ATG']); c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
        homo20 = {b: float(c[k] * np.log(20) + c[4 + k]) for k, b in enumerate('ATCG')}
        out['fits'][T] = dict(gilar_form=dict(r2=float(r2), loo_rmse_min=loo(X, y), homo_20mer_rt=homo20),
                              additive_C_reference=dict(per_full_fraction={'C': 0.0, 'A': float(c2[2]), 'T': float(c2[3]), 'G': float(c2[4])},
                                                        per_base_in_24mer={'C': 0.0, 'A': float(c2[2] / 24), 'T': float(c2[3] / 24), 'G': float(c2[4] / 24)},
                                                        loo_rmse_min=loo(X2, y)))
        pairs = {}
        order = 'CGAT'
        for i, j in itertools.combinations(range(len(R)), 2):
            if N[i] != N[j]: continue
            d = np.array([n[b][j] - n[b][i] for b in order])
            if np.abs(d).sum() != 2: continue
            frm, to = order[list(d).index(-1)], order[list(d).index(1)]; dy = y[j] - y[i]
            if order.index(frm) > order.index(to): frm, to, dy = to, frm, -dy
            pairs.setdefault(f'{frm}->{to}', []).append(float(dy))
        out['single_exchange_pairs'][T] = {k: dict(n=len(v), mean=float(np.mean(v)), sd=float(np.std(v, ddof=1)) if len(v) > 1 else None) for k, v in sorted(pairs.items())}
    open(a.output, 'w').write(json.dumps(out, indent=1) + '\n')
    for T in temps:
        f = out['fits'][T]; ad = f['additive_C_reference']['per_base_in_24mer']
        print(f"{T}: Gilar-form R2 {f['gilar_form']['r2']:.3f} LOO {f['gilar_form']['loo_rmse_min']:.3f} min | additive per base in a 24mer (C = 0): "
              + ', '.join(f'{b} {v:+.3f}' for b, v in sorted(ad.items(), key=lambda t: t[1])) + f" (LOO {f['additive_C_reference']['loo_rmse_min']:.3f})")
        print('   single exchanges (mean dRT, n): ' + ', '.join(f"{k} {v['mean']:+.2f} ({v['n']})" for k, v in out['single_exchange_pairs'][T].items()))


if __name__ == '__main__':
    main()
