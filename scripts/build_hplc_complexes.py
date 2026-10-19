#!/usr/bin/env python
"""Seeded donor-acceptor placements of two labelled monomers for the HPLC model study.

The guest donor X-H is aligned so H points at a host acceptor atom at a set H...A distance.
Approach directions are sampled on cones (0/30/50/70/85 deg, 8 azimuths) around the acceptor's
outward direction (away from its bonded neighbours), and the guest is rotated about each
approach axis in fixed steps. Outside the H-bond motif (acceptor and its neighbours versus donor X and
its neighbours), heavy atoms closer than 2.6 A or H...H closer than 1.8 A reject a placement;
the largest remaining heavy-atom clearance wins. The host keeps its coordinates; fragments, charges, spins,
bonds and frozen atoms are carried into complex metadata. These are optimization seeds only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from eda_gpu.model import load, xyz_text


def rotation(axis, angle):
    axis = axis / np.linalg.norm(axis); x, y, z = axis; c, s = np.cos(angle), np.sin(angle)
    return np.array([[c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
                     [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
                     [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])


def align(a, b):
    """Rotation taking unit vector a onto unit vector b."""
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v, c = np.cross(a, b), float(a @ b)
    if np.linalg.norm(v) < 1e-12:
        if c > 0: return np.eye(3)
        perp = np.cross(a, [1, 0, 0] if abs(a[0]) < .9 else [0, 1, 0])
        return rotation(perp, np.pi)
    return rotation(v, np.arctan2(np.linalg.norm(v), c))


def neighbours(meta, atom):
    return [j if i == atom else i for i, j, _ in meta['bonds'] if atom in (i, j)]


def place(hs, hc, hm, gs, gc, gm, acceptor, donor_h, distance, steps):
    hc, gc = np.asarray(hc, float), np.asarray(gc, float)
    (heavy,) = neighbours(gm, donor_h)
    outward = hc[acceptor] - hc[neighbours(hm, acceptor)].mean(axis=0)
    outward /= np.linalg.norm(outward)
    perp = np.cross(outward, [1, 0, 0] if abs(outward[0]) < .9 else [0, 1, 0]); perp /= np.linalg.norm(perp)
    heavy_h = np.array([s != 'H' for s in hs]); heavy_g = np.array([s != 'H' for s in gs])
    accepted = []
    # Contacts inside the H-bond motif (acceptor + neighbours vs donor X + neighbours) are fixed by the placement itself.
    pair = np.ones((len(hs), len(gs)), bool)
    pair[np.ix_([acceptor] + neighbours(hm, acceptor), [heavy] + neighbours(gm, heavy))] = False
    for tilt in (0, 30, 50, 70, 85):
        for azimuth in range(8 if tilt else 1):
            axis = rotation(outward, 2 * np.pi * azimuth / 8) @ (rotation(perp, np.radians(tilt)) @ outward)
            target_h = hc[acceptor] + distance * axis
            r = align(gc[donor_h] - gc[heavy], -axis)  # X-H points along -axis, H toward the acceptor
            base = (gc - gc[donor_h]) @ r.T + target_h
            for k in range(steps):
                g = (base - target_h) @ rotation(axis, 2 * np.pi * k / steps).T + target_h
                d = np.linalg.norm(hc[:, None, :] - g[None, :, :], axis=2)
                clearance = float(d[np.outer(heavy_h, heavy_g) & pair].min())
                if clearance < 2.6 or d[np.outer(~heavy_h, ~heavy_g) & pair].min() < 1.8:
                    continue
                accepted.append(((tilt, azimuth, k), g, clearance))
    return accepted


def complex_metadata(hm, gm, nh, identity):
    shift = lambda atoms: [a + nh for a in atoms]
    meta = {'identity': identity, 'protonation': f"host: {hm['protonation']}; guest: {gm['protonation']}",
            'charge': hm['charge'] + gm['charge'], 'spin': hm['spin'] + gm['spin'],
            'fragments': [dict(f) for f in hm['fragments']] + [dict(f, atoms=shift(f['atoms'])) for f in gm['fragments']],
            'bonds': hm['bonds'] + [[i + nh, j + nh, o] for i, j, o in gm['bonds']],
            'stereochemistry': f"host: {hm['stereochemistry']}; guest: {gm['stereochemistry']}",
            'geometry_provenance': 'seeded donor-acceptor placement of RDKit/MMFF94 monomer seeds; not optimized'}
    frozen = hm.get('frozen_atoms', []) + shift(gm.get('frozen_atoms', []))
    if frozen: meta['frozen_atoms'] = frozen
    return meta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host', required=True, help='monomer stem, e.g. examples/hplc/residue_A')
    p.add_argument('--host-xyz', help='host coordinates overriding <host>.xyz (e.g. an optimized complex)')
    p.add_argument('--guest', required=True)
    p.add_argument('--acceptor', type=int, nargs='+', required=True, help='host acceptor atom indices (0-based)')
    p.add_argument('--donor-h', type=int, required=True, help='guest donor hydrogen index (0-based)')
    p.add_argument('--distance', type=float, default=1.8)
    p.add_argument('--steps', type=int, default=12)
    p.add_argument('--keep', type=int, default=1, help='placements kept per acceptor (largest heavy-atom clearance)')
    p.add_argument('--output', required=True, help='output stem; _<acceptor>_<k>.xyz/.json appended')
    a = p.parse_args()
    hs, hc, hm = load(a.host_xyz or a.host + '.xyz', a.host + '.json')
    gs, gc, gm = load(a.guest + '.xyz', a.guest + '.json')
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    for acc in a.acceptor:
        found = sorted(place(hs, hc, hm, gs, gc, gm, acc, a.donor_h, a.distance, a.steps), key=lambda t: -t[2])[:a.keep]
        if not found: raise RuntimeError(f'no clash-free placement at acceptor {acc}')
        for n, (k, g, clearance) in enumerate(found):
            stem = Path(f'{out}_{acc}_{n}')
            meta = complex_metadata(hm, gm, len(hs), f"{hm['identity']} + {gm['identity']}")
            if a.host_xyz: meta['geometry_provenance'] = f'seeded placement on optimized host {a.host_xyz}; not optimized'
            meta['placement'] = {'host': a.host, 'host_geometry': a.host_xyz or a.host + '.xyz', 'guest': a.guest, 'host_acceptor': acc, 'guest_donor_h': a.donor_h + len(hs),
                                 'h_acceptor_distance_angstrom': a.distance, 'tilt_deg_azimuth_step_rotation_step': list(k), 'rotation_steps': a.steps,
                                 'min_heavy_atom_distance_angstrom': clearance}
            stem.with_suffix('.xyz').write_text(xyz_text(hs + gs, np.vstack([hc, g]), f'{stem.name}; placement seed; angstrom'))
            stem.with_suffix('.json').write_text(json.dumps(meta, indent=2) + '\n')
            print(f'{stem}: clearance {clearance:.2f} A')


if __name__ == '__main__':
    main()
