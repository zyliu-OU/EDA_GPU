#!/usr/bin/env python
"""Bridged ternary seeds: C8 proxy docked on the alkyl chains of an optimized residue...reagent ion pair.

The usual IP-RP picture: reagent alkyl chains adsorb on the bonded phase, the phosphate pairs with the
reagent N-H, and the residue does not touch the phase. The C8 proxy is placed as a rigid body around the
lowest-free-energy optimized ion pair per reagent (uniform random orientations, chain centroid 4-8 A from
the reagent carbon centroid; fixed RNG seed). Placements with heavy atoms closer than 3.0 A, H...H closer
than 2.0 A, or any proxy heavy atom within 4.0 A of a residue heavy atom are rejected; the placement with
the most proxy-carbon...reagent-carbon pairs within 4.5 A wins (ties: larger clearance). Seeds only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from eda_gpu.model import load, xyz_text
from build_hplc_complexes import complex_metadata

CONTACT, CLASH, HH, NO_RS = 4.5, 3.0, 2.0, 4.0


def random_rotation(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    return np.array([[1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)],
                     [2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)],
                     [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)]])


def dock(hs, hc, hm, gs, gc, reagent, trials, rng):
    hc, gc = np.asarray(hc, float), np.asarray(gc, float)
    frags = {f['name']: f['atoms'] for f in hm['fragments']}
    residue = next(a for n, a in frags.items() if n.startswith('residue'))
    r_heavy = [i for i in residue if hs[i] != 'H']
    t_carbon = [i for i in frags[reagent] if hs[i] == 'C']
    g_carbon = [i for i, s in enumerate(gs) if s == 'C']
    heavy_h, heavy_g = np.array([s != 'H' for s in hs]), np.array([s != 'H' for s in gs])
    centre, chain = hc[t_carbon].mean(axis=0), gc[g_carbon].mean(axis=0)
    best = None
    for _ in range(trials):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        g = (gc - chain) @ random_rotation(rng).T + centre + rng.uniform(4, 8) * u
        d = np.linalg.norm(hc[:, None] - g[None], axis=2)
        if d[np.outer(heavy_h, heavy_g)].min() < CLASH or d[np.outer(~heavy_h, ~heavy_g)].min() < HH: continue
        if d[np.ix_(r_heavy, np.where(heavy_g)[0])].min() < NO_RS: continue
        score = (int((d[np.ix_(t_carbon, g_carbon)] < CONTACT).sum()), float(d[np.outer(heavy_h, heavy_g)].min()))
        if best is None or score > best[0]: best = (score, g)
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pair', nargs=3, action='append', required=True, metavar=('XYZ', 'METADATA', 'REAGENT'),
                   help='optimized ion-pair geometry, its metadata, and the reagent fragment name (repeatable)')
    p.add_argument('--proxy', default='examples/hplc/c8_silane_proxy')
    p.add_argument('--trials', type=int, default=20000)
    p.add_argument('--output', default='examples/hplc/bridged')
    a = p.parse_args()
    gs, gc, gm = load(a.proxy + '.xyz', a.proxy + '.json')
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    for xyz, metadata, reagent in a.pair:
        hs, hc, hm = load(xyz, metadata)
        best = dock(hs, hc, hm, gs, gc, reagent, a.trials, np.random.default_rng(20261021))
        if best is None: raise RuntimeError(f'{xyz}: no acceptable bridged placement')
        (contacts, clearance), g = best
        residue = next(f['name'] for f in hm['fragments'] if f['name'].startswith('residue'))
        stem = out / f'{residue}_{reagent}_c8_bridged'
        meta = complex_metadata(hm, gm, len(hs), f"{hm['identity']} + {gm['identity']} (bridged)")
        meta['geometry_provenance'] = f'C8 proxy docked on optimized ion pair {xyz}; not optimized'
        meta['placement'] = {'ion_pair_geometry': xyz, 'ion_pair_metadata': metadata, 'proxy': a.proxy, 'reagent_fragment': reagent,
                             'trials': a.trials, 'rng_seed': 20261021, 'selection': 'most proxy-C...reagent-C pairs within 4.5 A',
                             'reagent_carbon_contacts': contacts, 'min_heavy_atom_distance_angstrom': clearance,
                             'min_residue_proxy_heavy_distance_required_angstrom': NO_RS}
        stem.with_suffix('.xyz').write_text(xyz_text(hs + gs, np.vstack([np.asarray(hc, float), g]), f'{stem.name}; bridged seed; angstrom'))
        stem.with_suffix('.json').write_text(json.dumps(meta, indent=2) + '\n')
        print(f'{stem}: {contacts} reagent-C contacts, clearance {clearance:.2f} A')


if __name__ == '__main__':
    main()
