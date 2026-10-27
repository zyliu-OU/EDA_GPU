#!/usr/bin/env python
"""Silanol-free three-chain alkyl patch and reagent-on-patch seeds for the HPLC model study (optional rdkit==2026.9.1).

Patch: three all-anti n-octane chains standing parallel to +z on an equilateral triangle of side 7.6 A (one
ligand per 50 A^2, i.e. ~3.3 umol/m^2 as on a C18 phase), each chain turned by 120 deg. C1 stands in for the
Si-CH2 tether; C1 and C2 of every chain are frozen, so the chains are tethered but can bend and collapse. No
silanol, siloxane or other polar site: the reagent can only contact alkyl chains (an endcapped C18 or a
polymeric alkyl phase). The central channel is 4.4 A from each chain axis, room for one inserted alkyl chain.

Seeds: each reagent monomer is placed as a rigid body (uniform random orientations, reagent carbon centroid
within 6 A of the patch axis and 3-12 A above the C1 plane; fixed RNG seed). Rejected: heavy atoms closer than
3.0 A, H...H closer than 2.0 A, any reagent atom below the C2 plane (the support surface). Ranked by the number
of reagent-C...patch-C pairs within 4.5 A (ties: clearance); kept: the best and the best whose carbon centroid
lies >= 3 A from every placement already kept. Seeds only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
from eda_gpu.model import load, xyz_text
from build_hplc_bridged import random_rotation
from build_hplc_complexes import complex_metadata
from build_hplc_models import OUT, bonds, embed

SIDE, CONTACT, CLASH, HH, SEED = 7.6, 4.5, 3.0, 2.0, 20261026
REAGENTS = ['teah', 'hexylammonium', 'dipeah', 'dibutylammonium', 'tributylammonium']


def octane():
    """All-anti n-octane, C1 at the origin, C1->C8 axis along +z."""
    mol = Chem.AddHs(Chem.MolFromSmiles('CCCCCCCC'))
    ids, energies = embed(mol, 50, 2026)
    cid = ids[int(np.argmin(energies))]
    conf = mol.GetConformer(cid)
    for k in range(5):
        T.SetDihedralDeg(conf, k, k + 1, k + 2, k + 3, 180.0)
    AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid).Minimize(maxIts=5000)
    if not all(abs(T.GetDihedralDeg(conf, k, k + 1, k + 2, k + 3)) >= 170 for k in range(5)):
        raise RuntimeError('octane did not stay all-anti')
    x = conf.GetPositions(); x -= x[0]
    axis = x[7] / np.linalg.norm(x[7])
    v, c = np.cross(axis, [0, 0, 1]), axis[2]
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return mol, x @ (np.eye(3) + k + k @ k / (1 + c)).T


def patch():
    mol, x = octane()
    n, symbols, coords, bond_list, frozen = mol.GetNumAtoms(), [], [], [], []
    for i in range(3):
        a = 2 * np.pi * i / 3
        rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        coords.append(x @ rz.T + SIDE / np.sqrt(3) * np.array([np.cos(a), np.sin(a), 0]))
        symbols += [atom.GetSymbol() for atom in mol.GetAtoms()]
        bond_list += [[p + i * n, q + i * n, o] for p, q, o in bonds(mol)]
        frozen += [i * n, i * n + 1]
    meta = {'identity': f'silanol-free C8 patch: three tethered all-anti n-octane chains, triangle side {SIDE} A',
            'protonation': 'neutral alkane; no polar site', 'stereochemistry': 'no stereocenters', 'charge': 0, 'spin': 0,
            'bonds': bond_list, 'frozen_atoms': frozen,
            'frozen_atoms_note': 'C1 and C2 of each chain: tethers on the support surface (C1 stands in for Si-CH2)',
            'fragments': [{'name': 'c8x3_patch', 'atoms': list(range(len(symbols))), 'charge': 0, 'spin': 0}],
            'geometry_provenance': 'RDKit ETKDG + MMFF94 all-anti octane, placed by construction; not quantum-optimized'}
    return symbols, np.vstack(coords), meta


def dock(ps, pc, gs, gc, trials, keep, rng):
    p_carbon = [i for i, s in enumerate(ps) if s == 'C']
    g_carbon = [i for i, s in enumerate(gs) if s == 'C']
    heavy_p, heavy_g = np.array([s != 'H' for s in ps]), np.array([s != 'H' for s in gs])
    floor = pc[[i for i in range(len(ps)) if i % (len(ps) // 3) == 1], 2].max()  # C2 plane
    centre = gc[g_carbon].mean(axis=0)
    found = []
    for _ in range(trials):
        r, z = 6 * np.sqrt(rng.uniform()), rng.uniform(3, 12)
        phi = rng.uniform(0, 2 * np.pi)
        g = (gc - centre) @ random_rotation(rng).T + [r * np.cos(phi), r * np.sin(phi), z]
        if g[:, 2].min() < floor: continue
        d = np.linalg.norm(pc[:, None] - g[None], axis=2)
        if d[np.outer(heavy_p, heavy_g)].min() < CLASH or d[np.outer(~heavy_p, ~heavy_g)].min() < HH: continue
        found.append(((int((d[np.ix_(p_carbon, g_carbon)] < CONTACT).sum()), float(d[np.outer(heavy_p, heavy_g)].min())), g))
    kept = []
    for score, g in sorted(found, key=lambda t: t[0], reverse=True):
        c = g[g_carbon].mean(axis=0)
        if all(np.linalg.norm(c - k[1][g_carbon].mean(axis=0)) >= 3 for k in kept): kept.append((score, g))
        if len(kept) == keep: break
    return kept


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reagents', nargs='+', default=REAGENTS)
    p.add_argument('--trials', type=int, default=50000)
    p.add_argument('--keep', type=int, default=2)
    p.add_argument('--output', default='examples/hplc/patch')
    a = p.parse_args()
    ps, pc, pm = patch()
    (OUT / 'c8x3_patch.xyz').write_text(xyz_text(ps, pc, 'c8x3_patch; constructed seed; angstrom'))
    (OUT / 'c8x3_patch.json').write_text(json.dumps(pm, indent=2) + '\n')
    print(f'c8x3_patch: {len(ps)} atoms')
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    for reagent in a.reagents:
        gs, gc, gm = load(OUT / f'{reagent}.xyz', OUT / f'{reagent}.json')
        kept = dock(ps, pc, gs, np.asarray(gc, float), a.trials, a.keep, np.random.default_rng(SEED))
        if not kept: raise RuntimeError(f'{reagent}: no acceptable placement')
        for n, ((contacts, clearance), g) in enumerate(kept):
            stem = out / f'patch_{reagent}_{n}'
            meta = complex_metadata(pm, gm, len(ps), f"{pm['identity']} + {gm['identity']}")
            meta['geometry_provenance'] = 'reagent seed docked on the constructed patch; not optimized'
            meta['placement'] = {'reagent': reagent, 'trials': a.trials, 'rng_seed': SEED, 'rank': n,
                                 'selection': f'most reagent-C...patch-C pairs within {CONTACT} A; kept placements >= 3 A apart',
                                 'reagent_carbon_contacts': contacts, 'min_heavy_atom_distance_angstrom': clearance,
                                 'reagent_N_height_angstrom': float(g[gs.index('N'), 2])}
            stem.with_suffix('.xyz').write_text(xyz_text(ps + gs, np.vstack([pc, g]), f'{stem.name}; patch seed; angstrom'))
            stem.with_suffix('.json').write_text(json.dumps(meta, indent=2) + '\n')
            print(f'{stem}: {contacts} C-C contacts, clearance {clearance:.2f} A, N at z {g[gs.index("N"), 2]:.1f} A, {len(ps) + len(gs)} atoms')


if __name__ == '__main__':
    main()
