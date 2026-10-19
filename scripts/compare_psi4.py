#!/usr/bin/env python
"""Cross-check EDA_GPU CPU/GPU energies and analytic gradients against an independent Psi4 run.

Both programs use exact (non-density-fitted) integrals, the same Libxc functional IDs,
the same basis representation, unpruned large grids, and identical bohr coordinates.
Tolerances are fixed before running; HF has no quadrature and is the tightest check.
This tests the DFT/gradient layer only; it is not an EDA component validation.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import numpy as np
from pyscf import gto, scf
from pyscf.dft import libxc
from pyscf.data.nist import BOHR
from eda_gpu.engine import Settings, hardware, singlepoint
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
PSI4 = ROOT / '.psi4-env' / 'bin' / 'python'
GRID_LEVEL, RADIAL, SPHERICAL = 8, 150, 1202
TOL = {'hf': (1e-8, 1e-7), 'dft': (2e-6, 2e-5), 'mgga': (2e-5, 5e-5)}  # (energy Eh, max gradient Eh/bohr)
# PySCF (conv_tol, conv_tol_grad); 1e-12/1e-7 is below the GPU meta-GGA noise floor for 425 functions,
# while 1e-6 orbital-gradient norm still bounds the variational energy error near 1e-12 Eh.
SCF_TOL = {'hf': (1e-12, 1e-7), 'dft': (1e-12, 1e-7), 'mgga': (1e-10, 1e-6)}
TRYP = ('tryptamine_water_neutral', ROOT / 'reports/tryptamine_neutral_gpu_2026-10-11/initial.xyz',
        ROOT / 'examples/tryptamine_water/neutral.json')
CASES = [  # name, xyz, metadata, basis, cart, method, tolerance class, psi4 scf_type (or 'df:<auxbasis>')
    ('water', ROOT / 'examples/water.xyz', ROOT / 'examples/water.json', 'def2-svp', False, 'pbe0', 'dft', 'pk'),
    ('water', ROOT / 'examples/water.xyz', ROOT / 'examples/water.json', 'def2-svp', False, 'b3lyp', 'dft', 'pk'),
    ('water_dimer', ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json', 'def2-svp', False, 'hf', 'hf', 'pk'),
    ('water_dimer', ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json', 'def2-svp', False, 'pbe0', 'dft', 'pk'),
    ('water_dimer', ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json', 'def2-svp', False, 'b3lyp', 'dft', 'pk'),
    ('co2_dimer_xeda_ii', ROOT / 'examples/co2_dimer_xeda_ii.xyz', ROOT / 'examples/co2_dimer_xeda_ii.json', '6-31g*', True, 'b3lyp', 'dft', 'pk'),
    (*TRYP, 'aug-cc-pvdz', False, 'm06-2x', 'mgga', 'direct'),
]
WATER_CATION = ('water_cation', ROOT / 'examples/water.xyz', ROOT / 'examples/water_cation.json')
EXTENDED = [  # density-fitted paths used by validate/optimize, and unrestricted open-shell gradients
    ('water_dimer', ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json', 'def2-svp', False, 'b3lyp', 'dft', 'df:def2-universal-jkfit'),
    (*TRYP, 'aug-cc-pvdz', False, 'm06-2x', 'mgga', 'df:aug-cc-pvdz-jkfit'),
    (*WATER_CATION, 'def2-svp', False, 'pbe0', 'dft', 'pk'),
    (*WATER_CATION, 'def2-svp', False, 'm06-2x', 'mgga', 'pk'),
]


def pyscf_hf(symbols, coords, meta, basis, cart):
    mol = gto.M(atom=list(zip(symbols, coords)), unit='Angstrom', charge=meta['charge'], spin=meta['spin'], basis=basis, cart=cart, verbose=0)
    mf = scf.RHF(mol) if not meta['spin'] else scf.UHF(mol)
    mf.conv_tol = 1e-12
    energy = mf.kernel()
    if not mf.converged: raise RuntimeError('PySCF HF did not converge')
    return {'energy_hartree': float(energy), 'gradient_hartree_per_bohr': mf.nuc_grad_method().kernel().tolist(), 'backend_used': 'cpu', 'note': 'direct PySCF HF, not the eda_gpu DFT CLI path'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--gpu', action='store_true', help='also run the actual GPU backend; never substituted by CPU')
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--skip-large', action='store_true', help='omit tryptamine-water M06-2X/aug-cc-pVDZ')
    parser.add_argument('--only', nargs='+', help='run only these system names')
    parser.add_argument('--set', choices=('core', 'extended', 'all'), default='core',
                        help='core = archived 2026-10-18 cases; extended = density fitting and UKS')
    args = parser.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=False)
    report = {'hardware': hardware(), 'gpu_requested': args.gpu, 'gpu_executed': False, 'tolerances': TOL, 'pyscf_scf_tolerances': SCF_TOL,
              'grid': {'pyscf_level': GRID_LEVEL, 'pyscf_prune': None, 'psi4_radial': RADIAL, 'psi4_spherical': SPHERICAL, 'psi4_pruning': 'none'},
              'bohr_angstrom_pyscf': BOHR, 'pyscf_libxc_version': libxc.__version__, 'scope': 'DFT/HF energies and analytic gradients only; not EDA validation', 'cases': [], 'passed': False}
    passed = True
    cases = {'core': CASES, 'extended': EXTENDED, 'all': CASES + EXTENDED}[args.set]
    for name, xyz, meta_path, basis, cart, method, tol_class, scf_type in cases:
        auxbasis = scf_type.split(':', 1)[1] if scf_type.startswith('df:') else None
        if args.skip_large and name == TRYP[0]: continue
        if args.only and name not in args.only: continue
        symbols, coords, meta = load(xyz, meta_path)
        label = f'{name}_{method}_{basis.replace("*", "s")}' + (f'_df_{auxbasis}' if auxbasis else '')
        job = {'symbols': symbols, 'coordinates_bohr': (np.asarray(coords) / BOHR).tolist(), 'charge': meta['charge'],
               'multiplicity': meta['spin'] + 1, 'method': method, 'basis': basis, 'puream': not cart, 'scf_type': 'df' if auxbasis else scf_type, 'df_basis_scf': auxbasis,
               'radial': RADIAL, 'spherical': SPHERICAL, 'threads': args.threads, 'memory': '16 GB'}
        (out / f'{label}_psi4_job.json').write_text(json.dumps(job, indent=2) + '\n')
        print(f'{label}: psi4', flush=True)
        subprocess.run([str(PSI4), '-I', str(ROOT / 'scripts/psi4_reference.py'), str(out / f'{label}_psi4_job.json'), str(out / f'{label}_psi4.json')], check=True)
        psi = json.loads((out / f'{label}_psi4.json').read_text())
        runs = {}
        for backend in (['cpu', 'gpu'] if args.gpu else ['cpu']):
            if method == 'hf':
                if backend == 'gpu': continue
                runs[backend] = pyscf_hf(symbols, coords, meta, basis, cart)
                continue
            print(f'{label}: eda_gpu {backend}', flush=True)
            conv_tol, conv_tol_grad = SCF_TOL[tol_class]
            settings = Settings(backend=backend, xc=method, basis=basis, cart=cart, density_fitting=bool(auxbasis),
                                auxbasis=auxbasis or Settings.auxbasis, grid_level=GRID_LEVEL,
                                conv_tol=conv_tol, conv_tol_grad=conv_tol_grad, max_cycle=200, threads=args.threads)
            try:
                runs[backend] = singlepoint(symbols, coords, meta, settings)
            except RuntimeError as e:
                runs[backend] = {'error': str(e), 'settings': asdict(settings)}
                print(f'  {backend}: FAILED {e}', flush=True)
            (out / f'{label}_{backend}.json').write_text(json.dumps(runs[backend], indent=2, allow_nan=False) + '\n')
            if 'error' in runs[backend]: continue
            if runs[backend]['backend_used'] != backend: raise RuntimeError('backend mismatch')
            if backend == 'gpu': report['gpu_executed'] = True
        e_tol, g_tol = TOL[tol_class]
        case = {'system': name, 'method': method, 'basis': basis, 'cart': cart, 'auxbasis': auxbasis, 'spin': meta['spin'], 'natoms': len(symbols), 'nbf_psi4': psi['nbf'],
                'psi4': psi, 'eda_gpu': {}, 'energy_tolerance_hartree': e_tol, 'gradient_tolerance_hartree_per_bohr': g_tol}
        for backend, r in runs.items():
            if 'error' in r:
                passed = False
                case['eda_gpu'][backend] = {'passed': False, 'error': r['error'], 'settings': r['settings']}
                continue
            de = r['energy_hartree'] - psi['energy_hartree']
            dg = float(np.abs(np.asarray(r['gradient_hartree_per_bohr']) - np.asarray(psi['gradient_hartree_per_bohr'])).max())
            ok = abs(de) <= e_tol and dg <= g_tol
            passed &= ok
            case['eda_gpu'][backend] = {'energy_hartree': r['energy_hartree'], 'energy_minus_psi4_hartree': de,
                                        'max_abs_gradient_difference_hartree_per_bohr': dg, 'passed': ok,
                                        'settings': r.get('settings'), 'seconds': r.get('seconds'), 'note': r.get('note')}
            print(f'  {backend}: dE={de:+.3e} Eh  max|dG|={dg:.3e} Eh/bohr  {"PASS" if ok else "FAIL"}', flush=True)
        report['cases'].append(case)
        (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    report['passed'] = passed
    (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print('PASSED' if passed else 'FAILED')


if __name__ == '__main__':
    main()
