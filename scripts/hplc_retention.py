#!/usr/bin/env python
"""Compare measured IP-RP-HPLC retention with computed GKS-EDA(sol) descriptors (correlation, not prediction).

Input: a CSV in the format of examples/hplc/retention_template.csv (one row per injection; rows that
should be compared share a condition_group, i.e. same column, gradient, temperature and buffer except
for the variable studied) and one or more analysis.json files from scripts/analyze_hplc.py.

Descriptors per sample, from the lowest-free-energy placement of each group (kcal/mol):
- composition-weighted residue terms: sum_b n_b * E(residue_b...reagent) and sum_b n_b * E(residue_b...C8),
  total and per component, also divided by length (requires all four bases computed with that reagent);
- reagent terms: E(C8...reagent) and E(residue_A...reagent), compared across reagents.
Per condition group and descriptor: Pearson r, least-squares slope, and leave-one-out RMSE of a linear fit
against the RMSE of predicting the group mean. Groups with fewer than 4 samples are reported, not fitted.
Single local placements without sampling: treat agreement as a hypothesis, see docs/HPLC_MODELS.md.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

REAGENT_MODELS = {'teaa': 'teah', 'tea': 'teah', 'teah': 'teah', 'tea/hfip': 'teah',
                  'haa': 'hexylammonium', 'hexylammonium': 'hexylammonium',
                  'dbaa': 'dibutylammonium', 'dibutylammonium': 'dibutylammonium',
                  'tbuaa': 'tributylammonium', 'tbaa': 'tributylammonium', 'tributylammonium': 'tributylammonium',
                  'dipea': 'dipeah', 'dipea/hfip': 'dipeah', 'dipeah': 'dipeah'}
BASES = 'ACGU'
TERMS = ['electrostatic', 'exchange', 'repulsion', 'polarization', 'generalized_correlation', 'desolvation', 'total']
MIN_FIT = 4


def composition(row):
    seq = (row.get('sequence') or '').strip().upper().replace('T', 'U')
    if seq:
        if set(seq) - set(BASES): raise ValueError(f"{row['sample_id']}: sequence has non-ACGU letters")
        return {b: seq.count(b) for b in BASES}
    try:
        return {b: int(row[f'n_{b}']) for b in BASES}
    except (KeyError, ValueError):
        raise ValueError(f"{row['sample_id']}: give a sequence or all of n_A, n_C, n_G, n_U") from None


def descriptors(row, groups):
    reagent = REAGENT_MODELS.get((row.get('reagent') or '').strip().lower())
    if reagent is None: return {}, f"reagent {row.get('reagent')!r} has no model"
    counts, out, missing = composition(row), {}, []
    length = sum(counts.values())
    for label, pattern in ((f'residue...{reagent}', 'residue_{}_' + reagent), ('residue...C8', 'residue_{}_c8')):
        rows = {b: groups.get(pattern.format(b)) for b in BASES}
        if all(rows.values()):
            for t in TERMS:
                s = sum(counts[b] * rows[b][t] for b in BASES)
                out[f'sum {label} {t}'] = s
                out[f'per-nt {label} {t}'] = s / length
        else:
            missing.append(f"{label} for bases {''.join(b for b in BASES if not rows[b])}")
    for label, key in ((f'C8...reagent', f'c8_{reagent}'), ('residue_A...reagent', f'residue_A_{reagent}')):
        if key in groups:
            for t in TERMS: out[f'{label} {t}'] = groups[key][t]
        else: missing.append(key)
    return out, '; '.join(missing)


def fit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if np.ptp(x) == 0: return None
    r = float(np.corrcoef(x, y)[0, 1]); slope = float(np.polyfit(x, y, 1)[0])
    loo, base = [], []
    for i in range(len(x)):
        m = np.arange(len(x)) != i
        if np.ptp(x[m]) == 0: return None
        loo.append(np.polyval(np.polyfit(x[m], y[m], 1), x[i]) - y[i]); base.append(y[m].mean() - y[i])
    return dict(n=len(x), pearson_r=r, slope_min_per_kcal_mol=slope,
                loo_rmse_min=float(np.sqrt(np.mean(np.square(loo)))), loo_rmse_group_mean_min=float(np.sqrt(np.mean(np.square(base)))))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data', required=True, help='retention CSV')
    p.add_argument('--analysis', nargs='+', required=True, help='analysis.json files from analyze_hplc.py')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    groups = {}
    for f in a.analysis:
        for name, g in json.loads(Path(f).read_text())['groups'].items():
            if 'components_kcal_mol' in g: groups[name] = g['components_kcal_mol']
    samples, by_group = [], {}
    with open(a.data, newline='') as fh:
        for row in csv.DictReader(fh):
            if not (row.get('retention_time_min') or '').strip(): continue
            d, missing = descriptors(row, groups)
            samples.append(dict(sample_id=row['sample_id'], condition_group=row['condition_group'], reagent=row.get('reagent'),
                                retention_time_min=float(row['retention_time_min']), descriptors_kcal_mol=d, missing_models=missing))
            by_group.setdefault(row['condition_group'], []).append(samples[-1])
    results = {}
    for name, members in sorted(by_group.items()):
        keys = sorted(set.intersection(*(set(s['descriptors_kcal_mol']) for s in members))) if members else []
        entry = dict(samples=len(members), reagents=sorted({s['reagent'] for s in members}))
        if len(members) < MIN_FIT:
            entry['status'] = f'fewer than {MIN_FIT} samples; not fitted'
        elif not keys:
            entry['status'] = 'no descriptor available for every sample'
        else:
            y = [s['retention_time_min'] for s in members]
            fits = {k: f for k in keys if (f := fit([s['descriptors_kcal_mol'][k] for s in members], y))}
            entry.update(status='fitted', fits=dict(sorted(fits.items(), key=lambda t: t[1]['loo_rmse_min'])))
        results[name] = entry
        print(f"{name}: {len(members)} samples, {entry['status']}")
        for k, f in list(entry.get('fits', {}).items())[:5]:
            print(f"  {k:55s} r {f['pearson_r']:+.2f}  LOO RMSE {f['loo_rmse_min']:.2f} min (group mean {f['loo_rmse_group_mean_min']:.2f})")
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(data=a.data, analysis=a.analysis, interpretation='correlation of retention with computed '
                                   'descriptors of single local model placements; not a retention prediction. Many descriptors are '
                                   'tested at once, and when composition varies along one direction (e.g. A for U swaps) all '
                                   'composition descriptors are collinear and fit equally well, so a good fit does not identify the term.',
                                   condition_groups=results, samples=samples), indent=2) + '\n')


if __name__ == '__main__':
    main()
