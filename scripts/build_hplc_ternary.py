#!/usr/bin/env python
"""Ternary residue/TEAH+/C8-proxy seeds (stage 7.3) on the lowest optimized ion pair per base.

The silanol O-H of the C8 proxy is placed on the phosphate O not bound to TEAH+, using
build_hplc_complexes.py with the optimized ion pair as host geometry. Seeds only.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

PHOSPHATE_O = {'A': (17, 18), 'C': (16, 17), 'G': (13, 14), 'U': (16, 17)}
SILANOL_H = 39


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', default='runs/hplc_ionpair_cpcm_def2svp')
    p.add_argument('--output', default='examples/hplc/ternary')
    a = p.parse_args()
    cases = json.loads((Path(a.campaign) / 'summary.json').read_text())['cases']
    for base, oxygens in PHOSPHATE_O.items():
        done = {o: cases[f'residue_{base}_teah_{o}_0'] for o in oxygens
                if cases.get(f'residue_{base}_teah_{o}_0', {}).get('bonding_topology_preserved')}
        if not done:
            print(f'{base}: no optimized ion pair', file=sys.stderr); continue
        bound = min(done, key=lambda o: done[o]['optimized_energy_hartree'])
        free = next(o for o in oxygens if o != bound)
        name = f'residue_{base}_teah_{bound}_0'
        subprocess.run([sys.executable, 'scripts/build_hplc_complexes.py', '--host', f'examples/hplc/ionpair/{name}',
                        '--host-xyz', f'{a.campaign}/{name}/optimization/optimized.xyz', '--guest', 'examples/hplc/c8_silane_proxy',
                        '--acceptor', str(free), '--donor-h', str(SILANOL_H), '--output', f'{a.output}/residue_{base}_teah_c8'], check=True)


if __name__ == '__main__':
    main()
