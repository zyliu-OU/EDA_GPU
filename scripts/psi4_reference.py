#!/usr/bin/env python
"""Independent Psi4 energy/analytic gradient for one job; run with the separate Psi4 environment.

Input JSON: symbols, coordinates_bohr, charge, multiplicity, method, basis, puream,
scf_type, radial, spherical, threads, memory. Coordinates are bohr so both programs
share one Bohr constant. Output JSON: energy, gradient and the Psi4/Libxc versions.
"""
import json
from pathlib import Path
import sys
import time
import psi4


def libxc_version():
    meta = sorted(Path(sys.prefix, 'conda-meta').glob('libxc-c-*.json'))
    return json.loads(meta[0].read_text())['version'] if meta else None


def main(job_path, out_path):
    job = json.loads(open(job_path).read())
    psi4.core.set_output_file(out_path.removesuffix('.json') + '.out', False)
    psi4.set_memory(job['memory'])
    psi4.set_num_threads(job['threads'])
    lines = [f"{job['charge']} {job['multiplicity']}", 'units bohr', 'symmetry c1', 'no_com', 'no_reorient']
    lines += [f'{s} {x:.12f} {y:.12f} {z:.12f}' for s, (x, y, z) in zip(job['symbols'], job['coordinates_bohr'])]
    psi4.geometry('\n'.join(lines))
    psi4.set_options({'basis': job['basis'], 'puream': job['puream'], 'scf_type': job['scf_type'],
                      'reference': 'rks' if job['multiplicity'] == 1 else 'uks',
                      'e_convergence': 1e-10, 'd_convergence': 1e-9, 'ints_tolerance': 1e-14,
                      'dft_radial_points': job['radial'], 'dft_spherical_points': job['spherical'],
                      'dft_pruning_scheme': 'none', 'guess': 'sad', 'maxiter': 200})
    if job['method'] == 'hf':
        psi4.set_options({'reference': 'rhf' if job['multiplicity'] == 1 else 'uhf'})
    start = time.perf_counter()
    grad, wfn = psi4.gradient(job['method'], return_wfn=True)
    result = {'energy_hartree': wfn.energy(), 'gradient_hartree_per_bohr': grad.np.tolist(),
              'nbf': wfn.basisset().nbf(), 'seconds': time.perf_counter() - start,
              'psi4_version': psi4.__version__, 'libxc_version': libxc_version(),
              'options': {k: job[k] for k in ('method', 'basis', 'puream', 'scf_type', 'radial', 'spherical')} |
                         {'e_convergence': 1e-10, 'd_convergence': 1e-9, 'ints_tolerance': 1e-14, 'pruning': 'none'}}
    open(out_path, 'w').write(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:3])
