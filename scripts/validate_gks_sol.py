#!/usr/bin/env python
"""Validation of experimental GKS-EDA(sol) beyond closure (which holds by construction).

1. Vacuum limit: eps -> 1 reproduces every gas-phase GKS-EDA component; desolvation -> 0.
2. Separation limit. Preregistered criterion (desolvation -> 0) was physically wrong and is kept
   as a recorded, non-counting failure: the gas-phase-form electrostatic term is the unscreened
   vacuum interaction, so desolvation must approach -(1-1/eps) x electrostatic while the total
   decays to zero. The corrected screening criterion was fixed before this rerun.
3. CPU/GPU component agreement for C-PCM and IEF-PCM.
4. PCM analytic gradient against central finite differences (solvated optimizations).
5. Exploratory comparison with the XEDA II manual BH3-NH3 C-PCM example (rounded values,
   unspecified radii/discretization/RI settings; not a matched validation).
Tolerances are fixed here before running. Run after sourcing scripts/gpu-env.sh.
"""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import numpy as np
from eda_gpu.engine import Settings, hardware, make_mf
from eda_gpu.gks import decompose
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]
TOL = {'vacuum_limit': 1e-6, 'cpu_gpu': 1e-7, 'fd_gradient': 2e-5, 'separation_fraction': 0.05,
       'screening_total_hartree': 1e-6, 'screening_ratio': 0.02}
XEDA_BH3NH3 = {'source': 'XEDA II manual Section 1.2.7, PDF sha256 269131a8383e117af77fd32a806d22ec0c8cd5a08213477c237ba59a76d89952',
               'settings': 'PBE0/6-31G, INT=LR (RI), SOLVENT=Water (C-PCM); radii and discretization unspecified',
               'components_hartree': {'electrostatic': -0.149777, 'exchange': -0.219316, 'repulsion': 0.410406, 'polarization': -0.098240,
                                      'generalized_correlation': -0.023502, 'desolvation': -0.004018, 'total': -0.084447}}
BH3NH3 = [('B', (-1.4930097404, -0.5082255173, 0.1193861744)), ('H', (-0.2909569897, -0.6053374650, 0.0759344076)),
          ('H', (-1.9156069781, 0.5145748522, -0.3618103786)), ('H', (-2.0743035096, -1.4976304007, -0.2543430610)),
          ('N', (-1.8483952471, -0.3946836719, 1.7286231427)), ('H', (-2.8481727853, -0.3078087733, 1.8739751783)),
          ('H', (-1.4009475344, 0.4155452940, 2.1428980703)), ('H', (-1.5302644453, -1.2172120485, 2.2291780163))]


def components(r):
    c = dict(r['components_hartree']); c['total'] = r['counterpoise_interaction_energy_hartree']
    return c


def run(out, tag, symbols, coords, meta, settings):
    d = out / tag; d.mkdir()
    r = decompose(symbols, coords, meta, settings, d)
    (d / 'result.json').write_text(json.dumps(r, indent=2, allow_nan=False) + '\n')
    print(f'{tag}: ' + '  '.join(f'{k}={v:+.7f}' for k, v in components(r).items()), flush=True)
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True); p.add_argument('--gpu', action='store_true'); p.add_argument('--threads', type=int, default=16)
    a = p.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=False)
    report = {'hardware': hardware(), 'tolerances': TOL, 'checks': {}, 'passed': True}
    def record(name, ok, counts=True, **data):
        report['checks'][name] = dict(data, passed=bool(ok), counts_toward_result=counts)
        if counts: report['passed'] &= bool(ok)
        print(f'== {name}: {"PASS" if ok else "FAIL"}', flush=True)
        (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    s, c, m = load(ROOT / 'examples/water_dimer.xyz', ROOT / 'examples/water_dimer.json')
    base = Settings(backend='cpu', xc='pbe0', basis='def2-svp', density_fitting=False, conv_tol=1e-11, conv_tol_grad=1e-6,
                    threads=a.threads, compute_duty_cycle=1.0)
    # 1. vacuum limit
    gas = components(run(out, 'water_dimer_gas', s, c, m, base))
    near = components(run(out, 'water_dimer_eps1', s, c, m, replace(base, solvent='cpcm', solvent_eps=1.000001)))
    diff = {k: near[k] - gas.get(k, 0.0) for k in near}
    record('vacuum_limit', max(abs(v) for v in diff.values()) <= TOL['vacuum_limit'], eps=1.000001, difference_hartree=diff)
    # 2. separation limit along the O...O axis
    axis = np.asarray(c[3]) - np.asarray(c[0]); axis /= np.linalg.norm(axis)
    scan = []
    for shift in (0.0, 4.0, 8.0, 12.0):
        moved = np.asarray(c, dtype=float).copy(); moved[3:] += shift * axis
        r = components(run(out, f'water_dimer_cpcm_shift{shift:g}', s, moved, m, replace(base, solvent='cpcm')))
        scan.append({'shift_angstrom': shift, 'desolvation_hartree': r['desolvation'], 'electrostatic_hartree': r['electrostatic'], 'total_hartree': r['total']})
    contact = abs(scan[0]['desolvation_hartree'])
    decays = all(abs(scan[i + 1]['desolvation_hartree']) <= abs(scan[i]['desolvation_hartree']) + 1e-9 for i in range(len(scan) - 1))
    record('separation_desolvation_vanishes_preregistered', decays and abs(scan[-1]['desolvation_hartree']) <= TOL['separation_fraction'] * contact
           and abs(scan[-1]['total_hartree']) <= abs(scan[0]['total_hartree']) * TOL['separation_fraction'], counts=False, scan=scan,
           note='criterion physically wrong (unscreened gas-phase electrostatics); superseded by separation_screening')
    target = 1 - 1 / replace(base, solvent='cpcm').solvent_eps  # dielectric screening of a distant dipole pair
    ratios = [-x['desolvation_hartree'] / x['electrostatic_hartree'] for x in scan[2:]]
    totals_decay = all(abs(scan[i + 1]['total_hartree']) <= abs(scan[i]['total_hartree']) for i in range(len(scan) - 1))
    record('separation_screening', totals_decay and abs(scan[-1]['total_hartree']) <= TOL['screening_total_hartree']
           and all(abs(r - target) <= TOL['screening_ratio'] for r in ratios), scan=scan, ratio_minus_desolvation_over_electrostatic=ratios,
           expected_ratio=target)
    # 3. CPU/GPU agreement
    if a.gpu:
        for model in ('cpcm', 'iefpcm'):
            cpu = components(run(out, f'water_dimer_{model}_cpu', s, c, m, replace(base, solvent=model)))
            gpu = components(run(out, f'water_dimer_{model}_gpu', s, c, m, replace(base, solvent=model, backend='gpu')))
            d = {k: gpu[k] - cpu[k] for k in cpu}
            record(f'cpu_gpu_{model}', max(abs(v) for v in d.values()) <= TOL['cpu_gpu'], gpu_minus_cpu_hartree=d)
    # 4. PCM gradient finite differences (water, C-PCM)
    ws, wc, wm = load(ROOT / 'examples/water.xyz', ROOT / 'examples/water.json')
    sp = replace(base, solvent='cpcm', density_fitting=True, conv_tol=1e-12, conv_tol_grad=None)
    for backend in (['cpu', 'gpu'] if a.gpu else ['cpu']):
        st = replace(sp, backend=backend)
        mf = make_mf(ws, wc, wm, st); mf.kernel()
        g = mf.nuc_grad_method().kernel(); g = np.asarray(g.get() if hasattr(g, 'get') else g)
        h, fd = 1e-3, np.zeros_like(g)
        from pyscf.data.nist import BOHR
        for atom in range(len(ws)):
            for x in range(3):
                e = []
                for sign in (1, -1):
                    cc = np.asarray(wc, dtype=float).copy(); cc[atom, x] += sign * h * BOHR
                    mm = make_mf(ws, cc, wm, st); e.append(float(mm.kernel()))
                fd[atom, x] = (e[0] - e[1]) / (2 * h)
        err = float(np.abs(fd - g).max())
        record(f'pcm_gradient_fd_{backend}', err <= TOL['fd_gradient'], max_error_hartree_per_bohr=err, step_bohr=h)
    # 5. exploratory XEDA comparison
    bs = [a_ for a_, _ in BH3NH3]; bc = [x for _, x in BH3NH3]
    bm = {'identity': 'BH3-NH3 (XEDA II manual)', 'protonation': 'neutral', 'charge': 0, 'spin': 0,
          'fragments': [{'name': 'BH3', 'atoms': [0, 1, 2, 3], 'charge': 0, 'spin': 0}, {'name': 'NH3', 'atoms': [4, 5, 6, 7], 'charge': 0, 'spin': 0}]}
    ours = components(run(out, 'bh3nh3_cpcm_6-31g', bs, np.asarray(bc), bm, replace(base, basis='6-31g', solvent='cpcm')))
    ref = XEDA_BH3NH3['components_hartree']
    report['exploratory_xeda_bh3nh3'] = dict(XEDA_BH3NH3, eda_gpu_hartree=ours, eda_gpu_minus_xeda_hartree={k: ours[k] - v for k, v in ref.items()},
                                             note='not a validation: radii, surface discretization and RI settings differ or are unknown')
    print('XEDA BH3NH3 diff: ' + '  '.join(f'{k}={ours[k]-v:+.6f}' for k, v in ref.items()), flush=True)
    report['settings'] = asdict(base)
    (out / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print('PASSED' if report['passed'] else 'FAILED')


if __name__ == '__main__':
    main()
