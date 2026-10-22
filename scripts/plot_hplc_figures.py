"""Figures for the HPLC model study: optimized structures, pair-energy diagrams, EDA vs published retention.

Reads archived geometries/metadata only (reports/, examples/hplc/). Needs numpy + matplotlib.
Usage: python3 scripts/plot_hplc_figures.py [--output docs/figures]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
FRAG = {'residue': '#2a78d6', 'reagent': '#eb6834', 'surface': '#1baf7a'}   # validated categorical slots 1-3
RADIUS = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'Si': 1.11}
SIZE = {'H': 14, 'C': 46, 'N': 70, 'O': 62, 'P': 110, 'Si': 110}
HETERO = {'N', 'O', 'P', 'Si'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.unicode_minus': True, 'font.size': 9, 'text.color': INK, 'axes.labelcolor': INK2,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.edgecolor': GRID, 'figure.facecolor': SURFACE,
                     'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE})


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    el = [l.split()[0] for l in lines[2:2 + n]]
    return el, np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])


def fragment_roles(meta_path):
    """Map atom index -> 'residue' | 'reagent' | 'surface' from the fragment metadata."""
    role = {}
    for f in json.loads(Path(meta_path).read_text())['fragments']:
        r = 'residue' if f['name'].startswith('residue') else 'surface' if 'silane' in f['name'] else 'reagent'
        for a in f['atoms']:
            role[a] = r
    return role


def hbonds(el, X, role):
    """Inter-fragment X-H...O contacts (H on N or O, acceptor O) shorter than 2.2 A."""
    out = []
    for h in (i for i, e in enumerate(el) if e == 'H'):
        donors = [k for k, e in enumerate(el) if e in ('N', 'O') and np.linalg.norm(X[h] - X[k]) < 1.15]
        if not donors:
            continue
        for o in (k for k, e in enumerate(el) if e == 'O' and role[k] != role[h]):
            d = np.linalg.norm(X[h] - X[o])
            if d < 2.2:
                out.append((h, o, d))
    return out


def draw_structure(ax, xyz, meta, title, notes):
    el, X = read_xyz(ROOT / xyz)
    role = fragment_roles(ROOT / meta)
    X = X - X.mean(0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)       # principal axes: widest view in the page plane
    P = X @ vt.T
    n = len(el)
    for i in range(n):
        for j in range(i + 1, n):
            if role[i] == role[j] and np.linalg.norm(X[i] - X[j]) < 1.2 * (RADIUS[el[i]] + RADIUS[el[j]]):
                ax.plot(P[[i, j], 0], P[[i, j], 1], color=FRAG[role[i]], lw=1.3, alpha=0.75, zorder=2 + P[i, 2] * 0.01)
    for h, o, d in hbonds(el, X, role):
        ax.plot(P[[h, o], 0], P[[h, o], 1], color=INK, lw=1.0, ls=(0, (2, 2)), zorder=50)
        mid = P[[h, o], :2].mean(0)
        ax.text(mid[0], mid[1] + 0.35, f'{d:.2f} Å', fontsize=6.5, color=INK, ha='center', zorder=60,
                bbox=dict(boxstyle='round,pad=0.12', fc=SURFACE, ec='none', alpha=0.85))
    for i in np.argsort(P[:, 2]):
        hetero = el[i] in HETERO
        ax.scatter(P[i, 0], P[i, 1], s=SIZE[el[i]], color=FRAG[role[i]] if el[i] != 'H' else '#ffffff',
                   edgecolors=INK if hetero else FRAG[role[i]], linewidths=1.1 if hetero else 0.6,
                   zorder=10 + P[i, 2] * 0.01)
        if el[i] in ('N', 'P', 'Si'):
            ax.text(P[i, 0], P[i, 1], el[i], fontsize=5.5, color='#ffffff', ha='center', va='center',
                    fontweight='bold', zorder=40)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=9.5, loc='left', color=INK, fontweight='bold')
    ax.text(0.0, -0.02, notes, transform=ax.transAxes, fontsize=7.8, color=INK2, va='top', linespacing=1.4)


def fragment_legend(fig):
    labels = {'residue': 'RNA residue (capped nucleotide, phosphate −1)', 'reagent': 'ion-pair reagent (alkylammonium +1)',
              'surface': 'C8 stationary-phase proxy (octyl-Si–O–Si(OH)₃)'}
    handles = [plt.Line2D([], [], marker='o', ls='', color=FRAG[k], markersize=7, label=v) for k, v in labels.items()]
    handles.append(plt.Line2D([], [], color=INK, ls=(0, (2, 2)), label='inter-fragment H-bond (H···O)'))
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=8)


def structures(out):
    panels = [
        ('reports/hplc_ionpair_2026-10-19/residue_A_teah_17_0/optimized.xyz', 'examples/hplc/ionpair/residue_A_teah_17_0.json',
         'a  Ion pair: residue A···TEAH⁺', 'N–H···O–P salt bridge\nR···T  −18.8 kcal/mol'),
        ('reports/hplc_reagents_2026-10-21/surface/c8_teah_contact_14_0/optimized.xyz',
         'examples/hplc/reagents/surface/c8_teah_contact_14_0.json',
         'b  Reagent on surface: TEAH⁺···C8', 'N–H binds a silanol O (all 17 reagent···C8 pairs)\nT···S  −10.6 kcal/mol'),
        ('reports/hplc_ternary_2026-10-20/residue_A_teah_c8_18_0/optimized.xyz', 'examples/hplc/ternary/residue_A_teah_c8_18_0.json',
         'c  Silanol-bound ternary (direct mode)', 'residue H-bonds the surface directly\nR···S −13.8 · T···S −2.0 · ΔE₃ +1.00'),
        ('reports/hplc_bridged_2026-10-22/residue_A_teah_c8_bridged/optimized.xyz', 'examples/hplc/bridged/residue_A_teah_c8_bridged.json',
         'd  Bridged ternary, TEAH⁺ (ionic mode)', 'residue never touches C8\nR···T −19.2 · T···S −2.1 · R···S 0.0 · ΔE₃ +0.02'),
        ('reports/hplc_bridged_2026-10-22/residue_A_hexylammonium_c8_bridged/optimized.xyz',
         'examples/hplc/bridged/residue_A_hexylammonium_c8_bridged.json',
         'e  Bridged ternary, hexylammonium', 'linear chain on C8\nR···T −18.3 · T···S −3.7 · R···S −0.1 · ΔE₃ +0.01'),
        ('reports/hplc_bridged_2026-10-22/residue_A_tributylammonium_c8_bridged/optimized.xyz',
         'examples/hplc/bridged/residue_A_tributylammonium_c8_bridged.json',
         'f  Bridged ternary, tributylammonium', 'bulky tertiary amine\nR···T −24.7 · T···S −3.2 · R···S 0.0 · ΔE₃ +0.03'),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 9.2))
    for ax, p in zip(axes.flat, panels):
        draw_structure(ax, *p)
    fragment_legend(fig)
    fig.suptitle('Optimized model systems — M06-2X/def2-SVP, C-PCM water (2-D projection on principal axes)',
                 x=0.01, ha='left', fontsize=11, fontweight='bold')
    fig.tight_layout(rect=(0, 0.04, 1, 0.97), h_pad=3.5)
    fig.savefig(out / 'hplc_structures.png', dpi=150)
    plt.close(fig)


def pair_diagram(out):
    """Triangle diagrams: edge width proportional to the pair interaction (kcal/mol)."""
    cases = [('Silanol-bound (direct mode)', {'RT': -18.28, 'TS': -2.02, 'RS': -13.78}, '+1.00'),
             ('Bridged (ionic mode)', {'RT': -19.22, 'TS': -2.149, 'RS': 0.02}, '+0.02')]
    pos = {'R': (0.0, 0.0), 'T': (1.0, 0.0), 'S': (0.5, 0.86)}
    name = {'R': ('residue A', 'residue'), 'T': ('TEAH⁺', 'reagent'), 'S': ('C8 surface', 'surface')}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for ax, (title, e, e3) in zip(axes, cases):
        for k, v in e.items():
            (x0, y0), (x1, y1) = pos[k[0]], pos[k[1]]
            w = max(abs(v) * 0.75, 0.6)
            ax.plot([x0, x1], [y0, y1], color=INK2 if abs(v) > 0.5 else GRID, lw=w, solid_capstyle='round', zorder=1,
                    ls='-' if abs(v) > 0.5 else (0, (2, 2)))
            off = np.array([-(y1 - y0), x1 - x0]) * (0.13 if k != 'RT' else -0.13)
            ax.text((x0 + x1) / 2 + off[0], (y0 + y1) / 2 + off[1], f'{v:.1f}'.replace('-', '−') if abs(v) >= 0.05 else '0.0', ha='center', va='center',
                    fontsize=10, color=INK, fontweight='bold')
        for k, (x, y) in pos.items():
            ax.scatter(x, y, s=1500, color=FRAG[name[k][1]], zorder=3, edgecolors=SURFACE, linewidths=2)
            ax.text(x, y - 0.2 if k != 'S' else y + 0.2, name[k][0], ha='center', va='center', fontsize=9.5, color=INK)
        ax.text(0.5, -0.42, f'three-body ΔE₃ = {e3} kcal/mol', ha='center', fontsize=9, color=INK2)
        ax.set_title(title, fontsize=10.5, fontweight='bold', color=INK)
        ax.set_xlim(-0.35, 1.35); ax.set_ylim(-0.55, 1.15); ax.set_aspect('equal'); ax.axis('off')
    fig.suptitle('Pair interactions at the ternary geometries (GKS-EDA(sol) totals, kcal/mol; line width ∝ |ΔE|)',
                 fontsize=10, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out / 'hplc_pair_diagram.png', dpi=150)
    plt.close(fig)


def bar_panel(ax, labels, values, title, unit, note=None, missing=(), decimals=1):
    x = np.arange(len(labels))
    v = np.array([np.nan if l in missing else val for l, val in zip(labels, values)], float)
    ax.bar(x, np.nan_to_num(v), width=0.6, color=FRAG['residue'], zorder=2)
    for xi, vi, l in zip(x, v, labels):
        if np.isnan(vi):
            ax.text(xi, 0, 'n/a', ha='center', va='bottom', fontsize=8, color=INK2)
        else:
            ax.text(xi, vi, f'{vi:.{decimals}f}'.replace('-', '−'), ha='center', va='bottom' if vi >= 0 else 'top', fontsize=8, color=INK)
    ax.set_xticks(x, labels, fontsize=8.5)
    ax.axhline(0, color=INK2, lw=0.8, zorder=3)
    ax.grid(axis='y', color=GRID, lw=0.6, zorder=0)
    ax.set_title(title, fontsize=9.5, loc='left', fontweight='bold', color=INK)
    ax.set_ylabel(unit)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    lo, hi = np.nanmin(np.append(v, 0)), np.nanmax(np.append(v, 0))
    pad = 0.18 * (hi - lo)
    ax.set_ylim(lo - (pad if lo < 0 else 0), hi + (pad if hi > 0 else pad * 0.3))
    if note:
        ax.text(0, -0.2, note, transform=ax.transAxes, fontsize=7.5, color=INK2, va='top')


def comparison(out):
    reagents = ['TEA', 'DIPEA', 'DBA', 'HA', 'TBA']          # order of increasing retention in both papers
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6))
    bar_panel(axes[0, 0], reagents, [11.92, 14.96, 32.89, 35.26, 0], 'Experiment: % ACN at 20T elution',
              '% acetonitrile', 'Donegan et al. 2022, Table 4 (100 mM, C18)\nTBA: no value reported', missing={'TBA'})
    bar_panel(axes[0, 1], reagents, [-18.77, -19.08, -23.42, -17.74, -23.52], 'EDA: ion pair residue A···reagent',
              'kcal/mol', 'rank correlation with experiment ρ = 0.4\nHA pairs weakest but is among most retained')
    bar_panel(axes[0, 2], reagents, [-10.55, -12.08, -10.82, -12.52, -12.67], 'EDA: isolated reagent···C8 pair',
              'kcal/mol', 'dominated by N–H···silanol H-bond;\nspread < placement range (5.3 kcal/mol)')
    bar_panel(axes[1, 0], reagents, [-2.149, -2.588, -2.373, -3.746, -3.172], 'EDA: chain-only reagent···C8 (bridged)',
              'kcal/mol', 'ρ = 0.8 (n = 5, single seeds — not established)\nC-PCM has no cavitation / hydrophobic term')
    bar_panel(axes[1, 1], reagents, [125, 170, 181, 128, 223], 'Buried nonpolar SASA (reagent···C8)',
              'Å²', 'ρ = 0.7; tracks logD/bp (r ≈ 0.9)\nbut ranks linear HA with TEA', decimals=0)
    bar_panel(axes[1, 2], ['A', 'C', 'G', 'U'], [-14.92, -18.87, -15.43, -14.86], 'EDA: residue···C8 direct contact',
              'kcal/mol', 'experiment: C/G-rich elute first (T > A > G > C)\nEDA puts C strongest — opposite order')
    fig.suptitle('EDA descriptors vs published IP-RP retention (reagents ordered by measured retention, weakest → strongest)',
                 fontsize=10.5, x=0.01, ha='left', fontweight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=4.5, w_pad=2.5)
    fig.savefig(out / 'hplc_eda_vs_experiment.png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', default='docs/figures')
    out = ROOT / ap.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    structures(out)
    pair_diagram(out)
    comparison(out)
    print(f'wrote figures to {out}')
