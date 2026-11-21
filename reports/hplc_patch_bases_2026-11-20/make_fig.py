"""Figure: per-base GKS-EDA(sol) total on the C8 patch vs NTP retention time (Maurer 2025 Fig. 9). Reads base_retention.json."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).parent
SURF, INK, INK2, MUTED, GRID, BLUE, BLUE_L = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8984', '#e4e3df', '#2a78d6', '#9ec5f4'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'axes.facecolor': SURF, 'figure.facecolor': SURF, 'axes.grid': True,
    'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titleweight': 'bold', 'axes.titlesize': 10, 'axes.titlecolor': INK, 'savefig.dpi': 220})
d = json.loads((HERE / 'base_retention.json').read_text())
fig, ax = plt.subplots(figsize=(5.2, 3.4))
for b, v in d['bases'].items():
    t = v['retention_time_min']
    pts = [p['total'] for p in v['placements'].values()]
    ax.scatter([t] * len(pts), pts, s=22, color=BLUE_L, edgecolor=SURF, linewidth=1, zorder=2)
    y = v['descriptors']['total/boltzmann']
    ax.scatter([t], [y], s=64, color=BLUE, edgecolor=SURF, linewidth=2, zorder=3)
    ax.annotate(f"{b}  {y:.1f}".replace('-', '\u2212'), (t, y), xytext=(9, 0), textcoords='offset points', va='center', color=INK, fontsize=9)
ax.set_xlabel('NTP retention time, BEH C18 + BA/TPA (min)')
ax.set_ylabel('GKS-EDA(sol) total (kcal/mol)')
ax.set_xlim(1.3, 3.15)
ax.set_title('Base binding on the C8 patch vs nucleotide retention', loc='left', pad=16)
ax.text(0.0, 1.015, 'Large dots: Boltzmann average (298 K); small dots: the 5 placements.  Spearman ρ = −1 (n = 4)',
        transform=ax.transAxes, fontsize=7.5, color=MUTED, va='bottom', ha='left')
fig.savefig(HERE / 'base_vs_retention.png', bbox_inches='tight', facecolor=SURF)
