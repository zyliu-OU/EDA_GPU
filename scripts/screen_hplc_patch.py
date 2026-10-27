#!/usr/bin/env python
"""GFN2-xTB/ALPB screen of extra reagent-on-patch placements; picks seeds for DFT.

Each seed is pre-optimized exactly as in run_hplc_stage.py (same constraints and topology checks), then
dE_xTB = E(complex) - E(patch) - E(reagent) is taken from ALPB single points at the optimized complex geometry
(no fragment relaxation; a ranking descriptor only). Placements whose reagent carbon centroid and N both lie
within --same A of an already-kept placement (including the --known xTB geometries already run with DFT) are
the same basin and dropped. Per reagent the --pick lowest dE_xTB distinct placements are written to selection.json.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import numpy as np
from eda_gpu.model import load, xyz_text
from run_hplc_stage import ALPB, XTB, read_xyz, write, xtb_preoptimize

HARTREE_KCAL = 627.509474


def single_point(symbols, coords, chrg, eps, folder):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'in.xyz').write_text(xyz_text(symbols, coords))
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='1G')
    out = subprocess.run([str(XTB), 'in.xyz', '--gfn', '2', '--alpb', ALPB[eps], '--chrg', str(chrg)],
                         cwd=folder, capture_output=True, text=True, env=env, check=True).stdout
    (folder / 'xtb.out').write_text(out)
    return float(re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)', out)[-1])


def descriptor(symbols, coords, reagent_atoms):
    g = [i for i in reagent_atoms if symbols[i] == 'C']
    return coords[g].mean(axis=0), coords[[i for i in reagent_atoms if symbols[i] == 'N'][0]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('seeds', nargs='+', help='seed stems patch_<reagent>_<rank> (without .xyz/.json)')
    p.add_argument('--known', nargs='*', default=[], help='run directories whose xtb/xtbopt.xyz are already-computed placements')
    p.add_argument('--output', required=True)
    p.add_argument('--pick', type=int, default=3)
    p.add_argument('--same', type=float, default=1.5)
    p.add_argument('--solvent-eps', type=float, default=78.3553)
    a = p.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    path = out / 'screen.json'
    screen = json.loads(path.read_text()) if path.exists() else {'cases': {}}
    screen['method'] = (f'GFN2-xTB/ALPB({ALPB[a.solvent_eps]}) pre-optimization as run_hplc_stage.py; dE = complex - fragments at the '
                        f'complex geometry (kcal/mol); duplicates: reagent C centroid and N both within {a.same} A')
    for stem in a.seeds:
        name = Path(stem).name
        case = screen['cases'].setdefault(name, {'seed': stem})
        if 'dE_xtb_kcal_mol' in case or 'error' in case: continue
        symbols, coords, meta = load(stem + '.xyz', stem + '.json')
        try:
            xyz = xtb_preoptimize(symbols, coords, meta, a.solvent_eps, out / name / 'xtb')
            e = {}
            for label, atoms, chrg in [('complex', range(len(symbols)), meta['charge'])] + \
                    [(f['name'], f['atoms'], f['charge']) for f in meta['fragments']]:
                e[label] = single_point([symbols[i] for i in atoms], xyz[list(atoms)], chrg, a.solvent_eps, out / name / f'sp_{label}')
            case['dE_xtb_kcal_mol'] = (e['complex'] - sum(v for k, v in e.items() if k != 'complex')) * HARTREE_KCAL
        except Exception as exc:
            case['error'] = f'{type(exc).__name__}: {exc}'
        print(f"{name}: {case.get('dE_xtb_kcal_mol', case.get('error'))}", flush=True)
        write(path, screen)
    selection = {}
    for name, case in screen['cases'].items():
        reagent = name.rsplit('_', 1)[0]
        selection.setdefault(reagent, {'known': [], 'candidates': []})
        if 'dE_xtb_kcal_mol' in case: selection[reagent]['candidates'].append(name)
    for reagent, s in selection.items():
        meta = json.loads(Path(screen['cases'][s['candidates'][0]]['seed'] + '.json').read_text())
        g_atoms = meta['fragments'][1]['atoms']
        symbols = load(screen['cases'][s['candidates'][0]]['seed'] + '.xyz', screen['cases'][s['candidates'][0]]['seed'] + '.json')[0]
        kept = []
        for run in a.known:
            for d in sorted(Path(run).glob(f'{reagent}_[0-9]*')):
                kept.append((d.name, descriptor(symbols, read_xyz(d / 'xtb' / 'xtbopt.xyz'), g_atoms)))
                s['known'].append(d.name)
        picked, dropped = [], {}
        for name in sorted(s['candidates'], key=lambda n: screen['cases'][n]['dE_xtb_kcal_mol']):
            c, n = descriptor(symbols, read_xyz(out / name / 'xtb' / 'xtbopt.xyz'), g_atoms)
            twin = next((k for k, (c0, n0) in kept if np.linalg.norm(c - c0) < a.same and np.linalg.norm(n - n0) < a.same), None)
            if twin: dropped[name] = twin; continue
            kept.append((name, (c, n)))
            if len(picked) < a.pick: picked.append(name)
        s.update(picked=picked, picked_seeds=[screen['cases'][n]['seed'] for n in picked], duplicate_of=dropped)
        print(f"{reagent}: picked {picked}; duplicates {dropped}", flush=True)
    write(out / 'selection.json', selection)


if __name__ == '__main__':
    main()
