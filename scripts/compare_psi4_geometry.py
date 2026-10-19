#!/usr/bin/env python
"""Compare EDA_GPU (geomeTRIC, default production thresholds) and Psi4 (optking,
gau_tight) optimized minima from identical starting geometries and settings.

Exact integrals, PBE0/def2-SVP, unpruned large grids, bohr coordinates shared via
PySCF's Bohr constant. Tolerances are fixed before running and reflect geomeTRIC's
default convergence (gmax 4.5e-4 Eh/bohr) on stiff (water) and soft (dimer) surfaces.
Both programs start from the same seeded symmetry-broken geometry (the example dimer has
an exact mirror plane that traps symmetry-preserving optimizers on a Cs saddle), and each
final structure is classified by an analytic PySCF Hessian at the same settings.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
from pyscf import dft, gto, lib
from pyscf.data.nist import BOHR
from pyscf.hessian import thermo
from eda_gpu.engine import Settings, hardware, optimize
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
PSI4 = ROOT / '.psi4-env' / 'bin' / 'python'
CASES = [  # name, rmsd tolerance angstrom, extra distance checks (atom pairs, tolerance angstrom)
    ('water', 2e-3, [((0, 1), 1e-3), ((0, 2), 1e-3)]),
    ('water_dimer', 2e-2, [((0, 3), 5e-3)]),
]
ENERGY_TOL = 1e-6
PERTURB_ANGSTROM, SEED = 0.05, 20261018


def aligned_rmsd(a, b):
    a = a - a.mean(axis=0); b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3); correction[-1, -1] = np.linalg.det(u @ vt)
    return float(np.sqrt(np.mean(np.sum((a @ u @ correction @ vt - b) ** 2, axis=1))))


def imaginary_modes(symbols, coords, threads):
    lib.num_threads(threads)
    mol = gto.M(atom=list(zip(symbols, coords)), unit='Angstrom', basis='def2-svp', verbose=0)
    mf = dft.RKS(mol); mf.xc = 'pbe0'; mf.grids.level = 8; mf.grids.prune = None; mf.conv_tol = 1e-12
    mf.kernel()
    freq = thermo.harmonic_analysis(mol, mf.Hessian().kernel())['freq_wavenumber']
    return [float(-abs(f)) for f in freq if np.iscomplex(f)], float(np.min(np.real(freq[~np.iscomplex(freq)])))


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    return np.array([[float(v) for v in line.split()[1:4]] for line in lines[2:2 + int(lines[0])]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=False)
    report = {'hardware': hardware(), 'method': 'pbe0/def2-svp exact integrals', 'energy_tolerance_hartree': ENERGY_TOL,
              'eda_gpu_optimizer': 'geomeTRIC, default EDA_GPU thresholds', 'psi4_optimizer': 'optking gau_tight',
              'cases': [], 'passed': True}
    for name, rmsd_tol, distances in CASES:
        symbols, coords, meta = load(ROOT / f'examples/{name}.xyz', ROOT / f'examples/{name}.json')
        coords = np.asarray(coords) + np.random.default_rng(SEED).uniform(-PERTURB_ANGSTROM, PERTURB_ANGSTROM, np.shape(coords))
        job = {'task': 'optimize', 'symbols': symbols, 'coordinates_bohr': (np.asarray(coords) / BOHR).tolist(),
               'charge': meta['charge'], 'multiplicity': meta['spin'] + 1, 'method': 'pbe0', 'basis': 'def2-svp',
               'puream': True, 'scf_type': 'pk', 'radial': 150, 'spherical': 1202, 'threads': args.threads, 'memory': '8 GB'}
        (out / f'{name}_psi4_job.json').write_text(json.dumps(job, indent=2) + '\n')
        print(f'{name}: psi4 optimize', flush=True)
        subprocess.run([str(PSI4), '-I', str(ROOT / 'scripts/psi4_reference.py'), str(out / f'{name}_psi4_job.json'), str(out / f'{name}_psi4.json')], check=True)
        psi = json.loads((out / f'{name}_psi4.json').read_text())
        print(f'{name}: eda_gpu optimize', flush=True)
        run = out / f'{name}_eda_gpu'; run.mkdir()
        settings = Settings(backend='cpu', xc='pbe0', basis='def2-svp', density_fitting=False, grid_level=8,
                            conv_tol=1e-11, conv_tol_grad=1e-6, threads=args.threads)
        ours = optimize(symbols, coords, meta, settings, run)
        (run / 'result.json').write_text(json.dumps(ours, indent=2, allow_nan=False) + '\n')
        a = read_xyz(run / 'optimized.xyz'); b = np.asarray(psi['coordinates_bohr']) * BOHR
        rmsd = aligned_rmsd(a, b)
        dist = [{'atoms': list(p), 'eda_gpu': float(np.linalg.norm(a[p[0]] - a[p[1]])), 'psi4': float(np.linalg.norm(b[p[0]] - b[p[1]])), 'tolerance': t}
                for p, t in distances]
        de = ours['energy_hartree'] - psi['energy_hartree']
        hessian = {}
        for label, geom in (('eda_gpu', a), ('psi4', b)):
            imag, lowest = imaginary_modes(symbols, geom, args.threads)
            hessian[label] = {'imaginary_cm-1': imag, 'lowest_real_cm-1': lowest, 'minimum': not imag}
        ok = ours['optimization_converged'] and hessian['eda_gpu']['minimum'] and hessian['psi4']['minimum'] and rmsd <= rmsd_tol and abs(de) <= ENERGY_TOL and all(abs(d['eda_gpu'] - d['psi4']) <= d['tolerance'] for d in dist)
        report['passed'] &= ok
        report['cases'].append({'system': name, 'aligned_rmsd_angstrom': rmsd, 'rmsd_tolerance_angstrom': rmsd_tol, 'distances_angstrom': dist,
                                'energy_eda_gpu_hartree': ours['energy_hartree'], 'energy_psi4_hartree': psi['energy_hartree'],
                                'energy_difference_hartree': de, 'hessian_pbe0_def2svp': hessian, 'start_perturbation_angstrom': PERTURB_ANGSTROM, 'seed': SEED, 'eda_gpu_converged': ours['optimization_converged'],
                                'eda_gpu_final_max_gradient': float(np.abs(ours['gradient_hartree_per_bohr']).max()),
                                'psi4_final_max_gradient': float(np.abs(psi['gradient_hartree_per_bohr']).max()), 'passed': ok})
        print(f'  RMSD={rmsd:.2e} A  dE={de:+.2e} Eh  imag(ours/psi4)={hessian["eda_gpu"]["imaginary_cm-1"]}/{hessian["psi4"]["imaginary_cm-1"]}  ' + '  '.join(f"r{d['atoms']}: {d['eda_gpu'] - d['psi4']:+.1e}" for d in dist) + ('  PASS' if ok else '  FAIL'), flush=True)
        (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print('PASSED' if report['passed'] else 'FAILED')


if __name__ == '__main__':
    main()
