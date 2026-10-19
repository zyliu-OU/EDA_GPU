#!/usr/bin/env python
"""Resumable optimize -> GKS-EDA(sol) campaign over labelled complex seeds.

For each <seed>.xyz/.json: optional GFN2-xTB/ALPB pre-optimization of the seed (anchor internal
geometry constrained; rejected if any bond breaks or a hydrogen changes partner, e.g. proton transfer), GPU
optimization (density fitting, PCM, frozen atoms from metadata), then direct-integral counterpoise GKS-EDA(sol) at the optimized geometry with the
same functional, basis and solvent. Completed stages are skipped on restart; a failure is
recorded and the campaign continues. Optimized structures are stationary points, not
established minima, and EDA values are electronic descriptors, not retention free energies.
"""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import subprocess
import traceback
import numpy as np
from pyscf.data import radii
from pyscf.data.elements import charge
from eda_gpu.engine import Settings, hardware, optimize
from eda_gpu.gks import decompose
from eda_gpu.model import load, xyz_text

KCAL = 627.509474
XTB = Path(__file__).resolve().parents[1] / '.xtb-env' / 'bin' / 'xtb'
ALPB = {78.3553: 'water', 35.688: 'acetonitrile'}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    return np.array([[float(v) for v in line.split()[1:4]] for line in lines[2:2 + int(lines[0])]])


def topology_preserved(symbols, coords, meta):
    """Every metadata bond within 1.25x covalent radii and every H nearest its bonded partner."""
    c = np.asarray(coords); cov = lambda x: radii.COVALENT[charge(x)] * radii.BOHR
    partner = {}
    for i, j, _ in meta['bonds']:
        if np.linalg.norm(c[i] - c[j]) > 1.25 * (cov(symbols[i]) + cov(symbols[j])): return False
        if symbols[i] == 'H': partner[i] = j
        if symbols[j] == 'H': partner[j] = i
    heavy = [k for k, x in enumerate(symbols) if x != 'H']
    return all(min(heavy, key=lambda k: np.linalg.norm(c[h] - c[k])) == p for h, p in partner.items())


def xtb_preoptimize(symbols, coords, meta, eps, folder):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'seed.xyz').write_text(xyz_text(symbols, coords))
    control = f"$chrg {meta['charge']}\n$spin {meta['spin']}\n$opt\n maxcycle=3000\n$end\n"
    if meta.get('frozen_atoms'):
        # Constrain the anchor's internal geometry (rigid-body motion allowed); exact $fix fights ANCopt here.
        control += '$constrain\n force constant=1.0\n atoms: ' + ','.join(str(i + 1) for i in meta['frozen_atoms']) + '\n$end\n'
    (folder / 'xcontrol').write_text(control)
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_STACKSIZE='1G')
    with open(folder / 'xtb.out', 'w') as log:
        subprocess.run([str(XTB), 'seed.xyz', '--gfn', '2', '--alpb', ALPB[eps], '--opt', 'normal', '--input', 'xcontrol'],  # seed refinement only; DFT converges tightly
                       cwd=folder, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    if 'GEOMETRY OPTIMIZATION CONVERGED' not in (folder / 'xtb.out').read_text():
        raise RuntimeError('xTB pre-optimization did not converge')
    xyz = read_xyz(folder / 'xtbopt.xyz')
    if not topology_preserved(symbols, xyz, meta):
        raise RuntimeError('xTB pre-optimization changed bonding topology (bond broken or proton moved)')
    frozen = meta.get('frozen_atoms', [])
    if frozen:
        internal = lambda x: np.linalg.norm(x[frozen][:, None] - x[frozen][None], axis=2)
        if np.abs(internal(xyz) - internal(np.asarray(coords))).max() > 0.02:
            raise RuntimeError('xTB distorted the frozen anchor geometry')
    return xyz


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('seeds', nargs='+', help='seed stems (without .xyz/.json)')
    p.add_argument('--output', required=True)
    p.add_argument('--xc', default='m06-2x'); p.add_argument('--basis', default='def2-svp')
    p.add_argument('--auxbasis', default='def2-universal-jkfit')
    p.add_argument('--solvent', default='cpcm', choices=['none', 'cpcm', 'iefpcm']); p.add_argument('--solvent-eps', type=float, default=78.3553)
    p.add_argument('--maxsteps', type=int, default=200); p.add_argument('--threads', type=int, default=16)
    p.add_argument('--skip-eda', action='store_true')
    p.add_argument('--continue-unconverged', action='store_true', help='restart a not-converged optimization from its last geometry')
    p.add_argument('--preopt', choices=['none', 'xtb'], default='none', help='GFN2-xTB/ALPB seed refinement before DFT')
    p.add_argument('--gpu-temperature-limit', type=float, default=80.0, help='pause until 65 C at this GPU temperature; 0 disables')
    p.add_argument('--cpu-temperature-limit', type=float, default=85.0, help='pause until 75 C at this CPU (k10temp Tctl) temperature; 0 disables')
    a = p.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    common = dict(backend='gpu', xc=a.xc, basis=a.basis, auxbasis=a.auxbasis, grid_level=3, max_cycle=200, threads=a.threads,
                  compute_duty_cycle=1.0, gpu_temperature_limit=a.gpu_temperature_limit,
                  cpu_temperature_limit=a.cpu_temperature_limit, solvent=a.solvent, solvent_eps=a.solvent_eps)
    opt_settings = Settings(density_fitting=True, conv_tol=1e-10, **common)
    eda_settings = Settings(density_fitting=False, conv_tol=1e-9, conv_tol_grad=1e-6, **common)
    summary_path = out / 'summary.json'
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {'cases': {}}
    summary.update(hardware=hardware(), optimization_settings=asdict(opt_settings), eda_settings=asdict(eda_settings))
    for stem in a.seeds:
        name = Path(stem).name
        case = summary['cases'].setdefault(name, {'seed': stem})
        symbols, coords, meta = load(stem + '.xyz', stem + '.json')
        opt_dir, eda_dir = out / name / 'optimization', out / name / 'eda'
        try:
            if a.preopt == 'xtb' and case.get('preoptimization') != 'GFN2-xTB/ALPB converged':
                print(f'{name}: GFN2-xTB pre-optimization', flush=True)
                coords = xtb_preoptimize(symbols, coords, meta, a.solvent_eps, out / name / 'xtb')
                case['preoptimization'] = 'GFN2-xTB/ALPB converged'; write(summary_path, summary)
            elif a.preopt == 'xtb':
                coords = read_xyz(out / name / 'xtb' / 'xtbopt.xyz')
            if case.get('optimization') != 'converged':
                if a.continue_unconverged and case.get('optimization') == 'not converged' and (opt_dir / 'optimized.xyz').exists():
                    coords = read_xyz(opt_dir / 'optimized.xyz')  # last geometry of the capped run
                    case['restarted_from_unconverged'] = case.get('restarted_from_unconverged', 0) + 1
                if opt_dir.exists(): opt_dir.rename(opt_dir.with_name(f'optimization_failed_{len(list(opt_dir.parent.glob("optimization_failed_*")))}'))
                opt_dir.mkdir(parents=True)
                print(f'{name}: optimizing ({len(symbols)} atoms)', flush=True)
                r = optimize(symbols, coords, meta, opt_settings, opt_dir, maxsteps=a.maxsteps)
                write(opt_dir / 'result.json', r)
                case.update(optimization='converged' if r['optimization_converged'] else 'not converged',
                            optimized_energy_hartree=r['energy_hartree'], evaluations=len(r['evaluations']),
                            frozen_max_displacement_angstrom=r['frozen_max_displacement_angstrom'])
                write(summary_path, summary)
                if not r['optimization_converged']: raise RuntimeError('optimization not converged')
                _, final, _ = load(opt_dir / 'optimized.xyz', stem + '.json')
                case['bonding_topology_preserved'] = bool(topology_preserved(symbols, final, meta)); write(summary_path, summary)
                if not case['bonding_topology_preserved']:
                    raise RuntimeError('DFT optimization changed bonding topology; fragment definitions invalid, EDA skipped')
            if not a.skip_eda and case.get('eda') != 'success':
                _, xyz, _ = load(opt_dir / 'optimized.xyz', stem + '.json')
                if eda_dir.exists(): eda_dir.rename(eda_dir.with_name(f'eda_failed_{len(list(eda_dir.parent.glob("eda_failed_*")))}'))
                eda_dir.mkdir(parents=True)
                print(f'{name}: GKS-EDA(sol)', flush=True)
                r = decompose(symbols, xyz, meta, eda_settings, eda_dir)
                write(eda_dir / 'result.json', r)
                comps = dict(r['components_hartree']); comps['total'] = r['counterpoise_interaction_energy_hartree']
                case.update(eda='success', components_kcal_mol={k: v * KCAL for k, v in comps.items()},
                            closure_residual_hartree=r['total_closure_residual_hartree'])
                print(f'{name}: ' + '  '.join(f'{k}={v:+.2f}' for k, v in case['components_kcal_mol'].items()), flush=True)
        except Exception as exc:
            case['error'] = f'{type(exc).__name__}: {exc}'; case['traceback'] = traceback.format_exc()
            print(f'{name}: FAILED {case["error"]}', flush=True)
        else:
            case.pop('error', None); case.pop('traceback', None)
        write(summary_path, summary)


if __name__ == '__main__':
    main()
