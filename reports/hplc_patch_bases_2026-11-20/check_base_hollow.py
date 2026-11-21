"""Hollow-insertion check for the base placements on the 3-chain C8 patch (same frame and thresholds as
reports/hplc_patch_extended_2026-11-11/check_patch_cases.py, but the guest position is the base-ring heavy-atom
centroid instead of the reagent N). Inserted = lateral distance from the patch axis r < 2.5 A and height h < 8 A."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, 'scripts')
from compare_base_retention import read_xyz
R_IN, H_IN = 2.5, 8.0
camp = Path('runs/hplc_patch_bases_cpcm_def2svp'); seeds = Path('examples/hplc/patch_bases')
cases = json.loads((camp / 'summary.json').read_text())['cases']
out = {}
for n in sorted(cases):
    m = json.loads((seeds / f'{n}.json').read_text())
    frags = m['fragments'] if isinstance(m['fragments'], list) else eval(m['fragments'])
    bonds = m['bonds'] if isinstance(m['bonds'], list) else eval(m['bonds'])
    s, X = read_xyz(camp / n / 'optimization' / 'optimized.xyz')
    host = set(frags[0]['atoms']); guest = [i for i in range(len(s)) if i not in host]
    nb = {i: set() for i in range(len(s))}
    for i, j, *_ in bonds: nb[i].add(j); nb[j].add(i)
    gly = next((nn, c) for nn in guest if s[nn] == 'N' for c in nb[nn]
               if s[c] == 'C' and any(s[o] == 'O' and sum(s[k] == 'C' for k in nb[o]) == 2 for o in nb[c]))
    base, st = set(), [gly[0]]
    while st:
        a = st.pop()
        if a not in base: base.add(a); st += [b for b in nb[a] if b != gly[1]]
    ring = [i for i in base if s[i] != 'H']
    C1 = X[[0, 26, 52]]; cen = C1.mean(0); nrm = np.cross(C1[1] - C1[0], C1[2] - C1[0]); nrm /= np.linalg.norm(nrm)
    if np.dot(X[[7, 33, 59]].mean(0) - cen, nrm) < 0: nrm = -nrm
    def pos(p): v = p - cen; h = v @ nrm; return float(h), float(np.linalg.norm(v - h * nrm))
    hb, rb = pos(X[ring].mean(0))
    P = next(i for i in guest if s[i] == 'P'); hp, rp = pos(X[P])
    gmin = float(min((X[i] - cen) @ nrm for i in guest if s[i] != 'H'))
    out[n] = dict(base_h=hb, base_r=rb, P_h=hp, P_r=rp, guest_min_h=gmin, inserted=bool(rb < R_IN and hb < H_IN))
json.dump(out, open('reports/hplc_patch_bases_2026-11-20/hollow.json', 'w'), indent=1)
print(f'| placement | base h | base r | P h | P r | min guest h | inserted |')
for n, v in out.items():
    print(f"| {n} | {v['base_h']:.2f} | {v['base_r']:.2f} | {v['P_h']:.2f} | {v['P_r']:.2f} | {v['guest_min_h']:.2f} | {'yes' if v['inserted'] else ''} |")
