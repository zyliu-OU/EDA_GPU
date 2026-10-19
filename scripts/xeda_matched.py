#!/usr/bin/env python
"""Matched-settings GKS-EDA comparison with GAMESS-XEDA.

prepare: write GAMESS-XEDA inputs (keywords from the official XEDA_gms manual) and
         run EDA_GPU at the same geometry, basis representation, functional variant,
         counterpoise policy and a large grid.
compare: parse component lines (ES=, EX=, REP=, POL=, CORR=, DISP=, E=) from GAMESS
         logs named <case>.log and compare every component with preset tolerances.

XEDA prints components to 1e-6 Eh, so agreement below ~5e-7 Eh cannot be resolved.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import re
import numpy as np
from pyscf.data.nist import BOHR
from eda_gpu.engine import Settings, hardware
from eda_gpu.gks import decompose
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
GRID_LEVEL, NRAD, NLEB = 8, 150, 1202
TOL = {'gga': 2e-6, 'mgga': 1e-5}  # per component, Eh; fixed before any XEDA output exists
# PySCF (conv_tol, conv_tol_grad); 1e-11 is below the GPU M06-2X level-8 noise floor for the
# 425-function ghost-basis reference SCFs; 1e-9/1e-6 matches the archived tryptamine EDA.
SCF_TOL = {'gga': (1e-11, 1e-6), 'mgga': (1e-9, 1e-6)}
# case: (example stem, metadata, pyscf xc, GAMESS DFTTYP, pyscf basis, GAMESS $BASIS, cart, tolerance class, backends)
WATER = (ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json')
CO2 = (ROOT / 'examples/co2_dimer_xeda_ii.xyz', ROOT / 'examples/co2_dimer_xeda_ii.json')
TRYP = (ROOT / 'reports/tryptamine_neutral_gpu_2026-10-11/initial.xyz', ROOT / 'examples/tryptamine_water/neutral.json')
CASES = {
    'water_dimer_pbe0_ccpvdz': (*WATER, 'pbe0', 'PBE0', 'cc-pvdz', 'GBASIS=CCD', False, 'gga', ('cpu', 'gpu')),
    'water_dimer_b3lyp5_ccpvdz': (*WATER, 'b3lyp5', 'B3LYP', 'cc-pvdz', 'GBASIS=CCD', False, 'gga', ('cpu', 'gpu')),
    'co2_dimer_pbe0_631gs_cart': (*CO2, 'pbe0', 'PBE0', '6-31g*', 'GBASIS=N31 NGAUSS=6 NDFUNC=1', True, 'gga', ('cpu', 'gpu')),
    'co2_dimer_b3lyp5_631gs_cart': (*CO2, 'b3lyp5', 'B3LYP', '6-31g*', 'GBASIS=N31 NGAUSS=6 NDFUNC=1', True, 'gga', ('cpu', 'gpu')),
    'tryptamine_water_m062x_augccpvdz': (*TRYP, 'm06-2x', 'M06-2X', 'aug-cc-pvdz', 'GBASIS=ACCD', False, 'mgga', ('gpu',)),
}
LABELS = {'ES': 'electrostatic', 'EX': 'exchange', 'REP': 'repulsion', 'POL': 'polarization',
          'CORR': 'generalized_correlation', 'DISP': 'explicit_dispersion', 'E': 'total'}
NUCLEAR = {'H': 1, 'C': 6, 'N': 7, 'O': 8}


def gamess_input(name, symbols, coords, meta, dfttyp, basis, cart):
    frags = meta['fragments']
    order = [i for f in frags for i in f['atoms']]
    if order != list(range(len(symbols))):
        raise ValueError(f'{name}: XEDA requires fragment atoms contiguous and in order')
    bohr = np.asarray(coords) / BOHR
    atoms = '\n'.join(f' {s:<2} {float(NUCLEAR[s]):4.1f} {x:18.12f} {y:18.12f} {z:18.12f}' for s, (x, y, z) in zip(symbols, bohr))
    mult = ' '.join(str(f['spin'] + 1) for f in frags)
    return f""" $CONTRL SCFTYP={'RHF' if meta['spin'] == 0 else 'UHF'} RUNTYP=EDA DFTTYP={dfttyp}
  ICHARG={meta['charge']} MULT={meta['spin'] + 1} MAXIT=200 ISPHER={-1 if cart else 1}
  UNITS=BOHR ICUT=11 ITOL=30 $END
 $LMOEDA MATOM(1)={' '.join(str(len(f['atoms'])) for f in frags)}
  MCHARG(1)={' '.join(str(f['charge']) for f in frags)} MMULT(1)={mult}
  EDATYP=GKS SUPBAS=.T. $END
 $BASIS {basis} $END
 $DFT NRAD={NRAD} NLEB={NLEB} $END
 $SCF DIRSCF=.T. FDIFF=.F. DIIS=.T. SOSCF=.F. CONV=1.0D-08 $END
 $SYSTEM MWORDS=500 MEMDDI=1000 $END
 $DATA
EDA_GPU matched XEDA case {name}; coordinates in bohr with PySCF BOHR={BOHR}
C1
{atoms}
 $END
"""


def prepare(args):
    out = Path(args.output); out.mkdir(parents=True, exist_ok=False)
    (out / 'gamess_inputs').mkdir()
    manifest = {'hardware': hardware(), 'grid': {'pyscf_level': GRID_LEVEL, 'gamess_nrad': NRAD, 'gamess_nleb': NLEB},
                'tolerances_hartree': TOL, 'pyscf_scf_tolerances': SCF_TOL, 'bohr_angstrom': BOHR, 'cases': {}}
    for name, (xyz, meta_path, xc, dfttyp, basis, gbasis, cart, tol, backends) in CASES.items():
        if args.only and name not in args.only: continue
        symbols, coords, meta = load(xyz, meta_path)
        (out / 'gamess_inputs' / f'{name}.inp').write_text(gamess_input(name, symbols, coords, meta, dfttyp, gbasis, cart))
        case = {'xyz': str(Path(xyz).relative_to(ROOT)), 'pyscf_xc': xc, 'gamess_dfttyp': dfttyp, 'basis': basis,
                'gamess_basis': gbasis, 'cart': cart, 'tolerance_class': tol, 'eda_gpu': {}}
        for backend in backends:
            if args.inputs_only: break
            if backend == 'gpu' and not args.gpu: continue
            run = out / f'{name}_{backend}'; run.mkdir()
            settings = Settings(backend=backend, xc=xc, basis=basis, cart=cart, density_fitting=False, grid_level=GRID_LEVEL,
                                conv_tol=SCF_TOL[tol][0], conv_tol_grad=SCF_TOL[tol][1], max_cycle=200, threads=args.threads)
            print(f'{name}: eda_gpu {backend}', flush=True)
            r = decompose(symbols, coords, meta, settings, run)
            (run / 'result.json').write_text(json.dumps(r, indent=2, allow_nan=False) + '\n')
            comps = dict(r['components_hartree'])
            comps['total'] = r['counterpoise_interaction_energy_hartree']
            case['eda_gpu'][backend] = {'components_hartree': comps, 'settings': asdict(settings),
                                        'closure_residual_hartree': r['total_closure_residual_hartree'], 'seconds': r['seconds']}
            print('  ' + '  '.join(f'{k}={v:+.6f}' for k, v in comps.items()), flush=True)
        manifest['cases'][name] = case
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')


def parse_log(path):
    found = {}
    for line in Path(path).read_text(errors='replace').splitlines():
        m = re.search(r'(?:ELECTROSTATIC|EXCHANGE|REPULSION|POLARIZATION|DISPERSION|TOTAL INTERACTION) ENERGY'
                      r'\s+(ES|EX|REP|POL|DISP|E)=\s*(-?\d+\.\d+)|ELECTRON CORRELATION\s+(CORR)=\s*(-?\d+\.\d+)', line)
        if m:  # last decomposition table in the log wins
            key, value = (m.group(1), m.group(2)) if m.group(1) else (m.group(3), m.group(4))
            found[LABELS[key]] = float(value)
    return found


def compare(args):
    manifest = json.loads(Path(args.manifest).read_text())
    report = {'manifest': args.manifest, 'cases': {}, 'passed': True}
    for name, case in manifest['cases'].items():
        log = Path(args.logs) / f'{name}.log'
        if not log.exists():
            print(f'{name}: no log'); continue
        xeda = parse_log(log)
        if 'total' not in xeda: raise ValueError(f'{log}: no XEDA decomposition table found')
        tol = manifest['tolerances_hartree'][case['tolerance_class']]
        rows = {}
        for backend, ours in case['eda_gpu'].items():
            diff = {k: ours['components_hartree'].get(k, 0.0) - v for k, v in xeda.items()}
            ok = all(abs(d) <= tol for d in diff.values())
            report['passed'] &= ok
            rows[backend] = {'difference_hartree': diff, 'passed': ok}
            print(f'{name} {backend}: ' + '  '.join(f'{k}={d:+.2e}' for k, d in diff.items()) + ('  PASS' if ok else '  FAIL'))
        report['cases'][name] = {'xeda_hartree': xeda, 'tolerance_hartree': tol, 'backends': rows}
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print('PASSED' if report['passed'] else 'FAILED')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--output', required=True)
    p.add_argument('--gpu', action='store_true', help='run actual GPU cases; never substituted by CPU')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--only', nargs='+')
    p.add_argument('--inputs-only', action='store_true')
    p.set_defaults(func=prepare)
    c = sub.add_parser('compare')
    c.add_argument('--manifest', required=True)
    c.add_argument('--logs', required=True, help='directory of GAMESS-XEDA logs named <case>.log')
    c.add_argument('--output', required=True)
    c.set_defaults(func=compare)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
