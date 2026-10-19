#!/usr/bin/env python
"""Compare tryptamine-water H-bond distances/angles with Liu et al. (2017) Table 1."""
import argparse
import json
from pathlib import Path
import numpy as np
from eda_gpu.model import load

ROOT = Path(__file__).resolve().parents[1]


def measure(coords, metadata):
    c = np.asarray(coords); rows = []
    for name, (i, j) in metadata['comparison_atoms'].items():
        h, x = (i, j) if metadata['symbols_hint'][i] == 'H' else (j, i)  # H and its acceptor/donor partner
        partner = min((np.linalg.norm(c[h] - c[k]), k) for k in range(len(c)) if k not in (h, x))[1]  # covalent neighbour of H
        v1, v2 = c[x] - c[h], c[partner] - c[h]
        rows.append({'hbond': name, 'atoms': [x, h, partner], 'distance_angstrom': float(np.linalg.norm(c[i] - c[j])),
                     'angle_degree': float(np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2))))})
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('xyz'); p.add_argument('--metadata', required=True)
    p.add_argument('--case', choices=['A_H2O', 'A_cation_H2O'], required=True); p.add_argument('--output')
    a = p.parse_args()
    symbols, coords, meta = load(a.xyz, a.metadata); meta['symbols_hint'] = symbols
    ref = json.loads((ROOT / 'references/tryptamine_water_2017.json').read_text())
    rows = measure(coords, meta)
    for row, d, t in zip(rows, ref['published_geometry_distances_angstrom'][a.case], ref['published_geometry_angles_degree'][a.case]):
        row.update(paper_distance_angstrom=d, paper_angle_degree=t, distance_difference_angstrom=row['distance_angstrom'] - d,
                   angle_difference_degree=row['angle_degree'] - t)
        print(f"{row['hbond']}: {row['distance_angstrom']:.3f} A / {row['angle_degree']:.1f} deg  paper {d:.3f} / {t:.1f}  "
              f"diff {row['distance_difference_angstrom']:+.3f} A / {row['angle_difference_degree']:+.1f} deg")
    result = {'xyz': a.xyz, 'case': a.case, 'reference': ref['doi'], 'notes': ref['published_geometry_notes'], 'hbonds': rows}
    if a.output: Path(a.output).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
