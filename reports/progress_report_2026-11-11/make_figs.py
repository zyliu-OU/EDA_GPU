"""Figures for the EDA_GPU progress report (2026-11-11). Reads archived data only."""
import json, glob
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Polygon

ROOT = Path('/home/zyliu/Documents/EDA_GPU')
OUT = Path(__file__).parent / 'fig'
OUT.mkdir(exist_ok=True)

SURF, INK, INK2, MUTED, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
CAT = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
BLUE_SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
CRIT, GOOD = '#d03b3b', '#0ca30c'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2,
    'xtick.color': INK2, 'ytick.color': INK2, 'axes.facecolor': SURF, 'figure.facecolor': SURF,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.axisbelow': True,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.titleweight': 'bold',
    'axes.titlesize': 10, 'axes.titlecolor': INK, 'legend.frameon': False, 'savefig.dpi': 220,
})


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches='tight', facecolor=SURF)
    plt.close(fig)


# reagent order = published retention (least -> most retained, Donegan 2022 Table 4)
ORDER = ['diethylammonium', 'teah', 'dipeah', 'butylammonium', 'tripropylammonium', 'dibutylammonium',
         'hexylammonium', 'tributylammonium']
LAB = {'diethylammonium': 'DEA', 'teah': 'TEA', 'dipeah': 'DIPEA', 'butylammonium': 'BA',
       'tripropylammonium': 'TPA', 'dibutylammonium': 'DBA', 'hexylammonium': 'HA', 'tributylammonium': 'TBA'}
RCOL = {k: CAT[i] for i, k in enumerate(ORDER)}

cases = json.load(open(ROOT / 'reports/hplc_patch_extended_2026-11-11/cases.json'))
rank = json.load(open(ROOT / 'reports/hplc_patch_extended_2026-11-11/ranking8.json'))
fit = json.load(open(ROOT / 'reports/hplc_retention_fit_2026-11-11/fit.json'))


def reagent_of(name):
    return name[len('patch_'):].rsplit('_', 1)[0]


# ---------------------------------------------------------------- 1 workflow
def fig_workflow():
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    ax.set_xlim(0, 102); ax.set_ylim(0, 46); ax.axis('off')
    steps = [
        ('1  Build', 'capped RNA\nresidue, reagent\ncation, 3-chain\nC8 patch', CAT[0]),
        ('2  Seeds', 'RDKit seeds;\nGFN2-xTB/ALPB\nscreen, best 3\nof 8 kept', CAT[0]),
        ('3  Optimize', 'M06-2X/\ndef2-SVP, DF,\nC-PCM water,\nfrozen anchors', CAT[1]),
        ('4  GKS-EDA(sol)', 'elec, exch,\nrep, pol,\ngen. corr.,\ndesolvation', CAT[1]),
        ('5  Descriptors', 'SMD CDS,\nburied nonpolar\nSASA, Boltzmann\nweights', CAT[2]),
        ('6  Compare', 'Spearman vs\npublished order;\nLOO fit to\n% ACN', CAT[2]),
    ]
    w, gap = 14.6, 2.4
    for i, (t, body, c) in enumerate(steps):
        x = 1 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 14), w, 24, boxstyle='round,pad=0.3,rounding_size=1.2',
                                    fc='white', ec=c, lw=1.6))
        ax.add_patch(FancyBboxPatch((x, 33.5), w, 4.5, boxstyle='round,pad=0.3,rounding_size=1.2',
                                    fc=c, ec=c, lw=1.6))
        ax.text(x + w / 2, 35.8, t, ha='center', va='center', color='white', fontsize=7.0, weight='bold')
        ax.text(x + w / 2, 23.5, body, ha='center', va='center', color=INK, fontsize=6.6, linespacing=1.4)
        if i < len(steps) - 1:
            ax.annotate('', xy=(x + w + gap - 0.2, 26), xytext=(x + w + 0.5, 26),
                        arrowprops=dict(arrowstyle='-|>', color=INK2, lw=1.2))
    bands = [(1, 2 * w + gap, 'structure preparation (CPU)', CAT[0]),
             (1 + 2 * (w + gap), 2 * w + gap, 'electronic structure (GPU, RTX 4080 SUPER)', CAT[1]),
             (1 + 4 * (w + gap), 2 * w + gap, 'analysis (CPU)', CAT[2])]
    for x, ww, t, c in bands:
        ax.plot([x, x + ww], [9.5, 9.5], color=c, lw=3, solid_capstyle='round')
        ax.text(x + ww / 2, 5.5, t, ha='center', va='center', color=INK2, fontsize=6.6)
    ax.text(51, 43.5, 'Guards on every run: GPU cooldown 75→65 °C, CPU 85→75 °C, 50 % duty cycle; '
            'convergence, topology and EDA-closure checks', ha='center', fontsize=6.6, color=INK2, style='italic')
    save(fig, 'f1_workflow.png')


# ---------------------------------------------------------------- 2 validation
def fig_validation():
    comps = ['Electrostatic', 'Pauli', 'Polarization', 'Gen. corr.', 'Total']
    neu = dict(ours=[-17.279, 19.535, -7.048, -5.489, -10.282], paper=[-17.34, 19.69, -7.07, -5.54, -10.27])
    cat = dict(ours=[-22.33, 23.68, -9.62, -5.30, -13.56], paper=[-22.24, 23.48, -9.57, -5.22, -13.56])
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axs[0]
    x = np.arange(5); bw = 0.19
    for j, (nm, d, c) in enumerate([('neutral', neu, CAT[0]), ('radical cation', cat, CAT[1])]):
        ax.bar(x + (j - 0.5) * 2 * bw - bw / 2, d['ours'], bw * 0.92, color=c, label=f'EDA_GPU, {nm}')
        ax.bar(x + (j - 0.5) * 2 * bw + bw / 2, d['paper'], bw * 0.92, color='white', ec=c, lw=1.1, hatch='////',
               label=f'Liu et al. 2017, {nm}')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x, ['Elec.', 'Pauli', 'Pol.', 'Gen. corr.', 'Total'], fontsize=7.8); ax.set_ylabel('kcal/mol')
    ax.set_title('a  Tryptamine–water GKS-EDA components', loc='left')
    ax.legend(fontsize=6.6, ncol=2, loc='lower left', bbox_to_anchor=(0, -0.42))
    ax = axs[1]
    dn = np.array(neu['ours']) - neu['paper']; dc = np.array(cat['ours']) - cat['paper']
    ax.axvspan(-0.2, 0.2, color='#f0efec', zorder=0)
    ax.axvline(0, color=INK2, lw=0.8)
    y = np.arange(5)[::-1]
    ax.scatter(dn, y + 0.12, s=36, color=CAT[0], zorder=3, ec='white', lw=1, label='neutral')
    ax.scatter(dc, y - 0.12, s=36, color=CAT[1], zorder=3, ec='white', lw=1, label='cation')
    ax.set_yticks(y, comps, fontsize=7.8); ax.set_xlim(-0.3, 0.3)
    ax.set_xlabel('EDA_GPU − paper (kcal/mol)')
    ax.text(0.19, -0.75, '±0.2 band', ha='right', fontsize=6.5, color=INK2)
    ax.set_ylim(-1, 4.5)
    ax.set_title('b  Deviation from published', loc='left')
    ax.legend(fontsize=6.8, loc='lower right')
    fig.tight_layout()
    save(fig, 'f2_validation.png')


# ---------------------------------------------------------------- 3 scaling + memory
def fig_scaling():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]
    labels = ['water dimer\n48 bf', 'tryptamine–H₂O\ndef2-SVP 252 bf', 'tryptamine–H₂O\naug-cc-pVDZ 425 bf']
    eg = [(1.2, 0.5), (56.9, 8.4), (137.5, 22.0)]
    ed = [(2.4, 14.1), (169.5, 53.9), (548.6, 121.1)]
    x = np.arange(3); bw = 0.19
    for k, (data, tag, off) in enumerate([(eg, 'energy + gradient', -1), (ed, 'GKS-EDA', 1)]):
        cpu = [a for a, b in data]; gpu = [b for a, b in data]
        xs = x + off * bw * 1.1
        ax.bar(xs - bw / 2, cpu, bw * 0.92, color=CAT[0] if k == 0 else CAT[2], label=f'CPU, {tag}')
        ax.bar(xs + bw / 2, gpu, bw * 0.92, color=CAT[1] if k == 0 else CAT[3], label=f'GPU, {tag}')
        for xi, a, b in zip(xs, cpu, gpu):
            ax.text(xi, max(a, b) * 1.25, f'{a / b:.1f}×', ha='center', fontsize=6.8, color=INK)
    ax.set_yscale('log'); ax.set_ylim(0.3, 2500)
    ax.set_xticks(x, labels, fontsize=7.2); ax.set_ylabel('median wall time (s, log)')
    ax.set_title('a  CPU vs GPU timing (speedup = CPU/GPU)', loc='left')
    ax.legend(fontsize=6.4, ncol=2, loc='upper left')
    ax = axs[1]
    chains = ['3 chains\n119 atoms', '5 chains\n171 atoms', '7 chains\n223 atoms']
    buf = [0.9 * 2, 2.7 * 2, 6.0 * 2]
    cols = [CAT[0], CAT[0], CAT[7]]
    ax.bar(range(3), buf, 0.55, color=cols)
    ax.axhline(16, color=CRIT, lw=1.4, ls='--')
    ax.text(-0.35, 16.5, 'GPU memory 16 GB', fontsize=7, color=CRIT)
    ax.scatter([2], [15.3], marker='v', s=50, color=INK, zorder=4)
    ax.text(1.55, 13.4, 'observed peak\n15.3 GB → OOM', fontsize=6.6, color=INK, ha='right')
    for i, b in enumerate(buf):
        ax.text(i, b + 0.4, f'{b:.1f}', ha='center', fontsize=7)
    ax.set_ylim(0, 19); ax.set_xticks(range(3), chains, fontsize=7.2)
    ax.set_ylabel('fixed DF K-gradient buffers (GB)')
    ax.set_title('b  Patch size vs GPU memory (+TBA)', loc='left')
    fig.tight_layout()
    save(fig, 'f3_scaling_memory.png')


# ---------------------------------------------------------------- 4 patch EDA components
def fig_components():
    comp = [('elec', 'Electrostatic', CAT[0]), ('pol', 'Polarization', CAT[2]), ('gc', 'Gen. correlation', CAT[6]),
            ('pauli', 'Pauli (exch + rep)', CAT[1]), ('des', 'Desolvation', CAT[3])]
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    x = np.arange(len(ORDER)); bw = 0.56
    for i, r in enumerate(ORDER):
        d = rank['reagents'][r]
        pos = neg = 0.0
        for k, nm, c in comp:
            v = d[k]['boltzmann']
            if v >= 0:
                ax.bar(i, v, bw, bottom=pos, color=c, ec=SURF, lw=1.2, label=nm if i == 0 else None); pos += v
            else:
                ax.bar(i, v, bw, bottom=neg, color=c, ec=SURF, lw=1.2, label=nm if i == 0 else None); neg += v
        tot = d['tot']['boltzmann']
        ax.scatter(i, tot, marker='D', s=34, color=INK, zorder=5, ec='white', lw=1,
                   label='EDA total' if i == 0 else None)
        ax.text(i + 0.34, tot, f'{tot:.1f}', va='center', fontsize=7, color=INK)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x, [LAB[r] for r in ORDER])
    ax.set_ylabel('kcal/mol (Boltzmann over 5 placements)')
    ax.set_xlabel('reagent, ordered by published retention (least → most retained) →', fontsize=8)
    ax.set_title('Reagent⁺···C8-patch GKS-EDA(sol) decomposition, M06-2X/def2-SVP, C-PCM water', loc='left')
    ax.legend(fontsize=7, ncol=6, loc='upper center', bbox_to_anchor=(0.5, -0.2))
    save(fig, 'f4_patch_components.png')


# ---------------------------------------------------------------- 5 patch geometry / hollow
def read_xyz(p):
    L = Path(p).read_text().split('\n'); n = int(L[0]); s = []; c = []
    for l in L[2:2 + n]:
        t = l.split(); s.append(t[0]); c.append([float(v) for v in t[1:4]])
    return s, np.array(c)


def frame(X):
    C1 = X[[0, 26, 52]]; cen = C1.mean(0)
    nrm = np.cross(C1[1] - C1[0], C1[2] - C1[0]); nrm /= np.linalg.norm(nrm)
    if np.dot(X[[7, 33, 59]].mean(0) - cen, nrm) < 0: nrm = -nrm
    e1 = C1[0] - cen; e1 -= (e1 @ nrm) * nrm; e1 /= np.linalg.norm(e1); e2 = np.cross(nrm, e1)
    return cen, nrm, e1, e2


def xyz_of(name):
    return next(iter(glob.glob(str(ROOT / f'reports/hplc_patch_*/{name}/optimized.xyz'))))


def fig_hollow():
    fig = plt.figure(figsize=(7.2, 3.6))
    ax = fig.add_axes([0.06, 0.12, 0.42, 0.78])
    names = list(cases)
    pts = []
    for nm in names:
        s, X = read_xyz(xyz_of(nm)); cen, nrm, e1, e2 = frame(X)
        N = [i for i in range(78, len(s)) if s[i] == 'N'][0]; v = X[N] - cen
        pts.append((nm, v @ e1, v @ e2, cases[nm]['tot'], cases[nm]['h'], cases[nm]['r']))
        if nm == names[0]:
            anchors = [((X[i] - cen) @ e1, (X[i] - cen) @ e2) for i in (0, 26, 52)]
            tips = [((X[i] - cen) @ e1, (X[i] - cen) @ e2) for i in (7, 33, 59)]
    ax.add_patch(Polygon(anchors, closed=True, fc='#f0efec', ec=MUTED, lw=0.8, zorder=0))
    ax.add_patch(Circle((0, 0), 2.5, fc='none', ec=CRIT, lw=1.1, ls='--', zorder=1))
    ax.text(0, -2.95, 'hollow: r < 2.5 Å', ha='center', fontsize=6.8, color=CRIT)
    for a in anchors:
        ax.scatter(*a, s=70, marker='s', color=INK2, zorder=2)
    ax.scatter([], [], s=40, marker='s', color=INK2, label='chain anchor (C1)')
    vals = np.array([p[3] for p in pts])
    sc = None
    for nm, x, y, tot, h, r in pts:
        inside = (r < 2.5 and h < 8.0)
        sc = ax.scatter(x, y, c=[tot], cmap=matplotlib.colors.ListedColormap(BLUE_SEQ[::-1]),
                        vmin=-9, vmax=-2.5, s=52 if inside else 34, marker='o' if inside else '^',
                        ec=INK if inside else 'white', lw=1.0, zorder=3)
    ax.scatter([], [], marker='o', s=40, color=BLUE_SEQ[4], ec=INK, label='N inserted in hollow (h < 8 Å)')
    ax.scatter([], [], marker='^', s=34, color=BLUE_SEQ[2], ec='white', label='other placement')
    ax.set_aspect('equal'); ax.set_xlim(-8, 8); ax.set_ylim(-8, 8.5)
    ax.set_xlabel('in-plane x (Å)'); ax.set_ylabel('in-plane y (Å)')
    ax.set_title('a  Top view: reagent N of all 40 placements', loc='left', fontsize=9)
    ax.legend(fontsize=6.4, loc='upper left', bbox_to_anchor=(-0.02, 1.02))
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02); cb.set_label('EDA total (kcal/mol)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5)

    ax2 = fig.add_axes([0.62, 0.14, 0.36, 0.74])
    ys = np.arange(len(ORDER))[::-1]
    for i, r in enumerate(ORDER):
        rows = [(nm, cases[nm]) for nm in names if reagent_of(nm) == r]
        for nm, c in rows:
            inside = c['r'] < 2.5 and c['h'] < 8.0
            ax2.scatter(c['tot'], ys[i], marker='o' if inside else '^', s=40 if inside else 30,
                        color=RCOL[r] if inside else 'white', ec=RCOL[r], lw=1.3, zorder=3)
        ax2.scatter(rank['reagents'][r]['tot']['boltzmann'], ys[i], marker='|', s=160, color=INK, zorder=4)
    ax2.set_yticks(ys, [LAB[r] for r in ORDER]); ax2.set_xlabel('EDA total per placement (kcal/mol)')
    ax2.set_title('b  Placements per reagent', loc='left', fontsize=9)
    ax2.text(-9.3, -1.25, '● hollow   △ outside   | Boltzmann', fontsize=6.8, color=INK2)
    ax2.set_ylim(-1.6, 7.6)
    ax2.text(1.02, 0.5, 'published retention ↑', transform=ax2.transAxes, rotation=90, va='center',
             fontsize=7, color=INK2)
    save(fig, 'f5_hollow.png')


def fig_structures():
    """Side views of an inserted (DEA 6) and outside (TEA 2) placement."""
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9))
    sel = [('patch_diethylammonium_6', 'a  DEA⁺, placement 6: in hollow'),
           ('patch_dibutylammonium_0', 'b  DBA⁺, placement 0: inserted'),
           ('patch_teah_2', 'c  TEA⁺, placement 2: outside')]
    rad = {'H': 9, 'C': 30, 'N': 55}
    for ax, (nm, t) in zip(axs, sel):
        s, X = read_xyz(xyz_of(nm)); cen, nrm, e1, e2 = frame(X)
        P = np.c_[(X - cen) @ e1, (X - cen) @ nrm, (X - cen) @ e2]
        order = np.argsort(P[:, 2])
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                d = np.linalg.norm(X[i] - X[j])
                lim = 1.2 if 'H' in (s[i], s[j]) else 1.7
                if d < lim:
                    col = CAT[2] if i < 78 and j < 78 else CAT[1]
                    ax.plot(P[[i, j], 0], P[[i, j], 1], color=col, lw=1.0 if 'H' in (s[i], s[j]) else 1.8,
                            alpha=0.85, zorder=1)
        for i in order:
            col = CAT[2] if i < 78 else (CAT[0] if s[i] == 'N' else CAT[1])
            ax.scatter(P[i, 0], P[i, 1], s=rad[s[i]], color=col if s[i] != 'H' else 'white',
                       ec=col, lw=0.8, zorder=2)
        ax.axhline(0, color=MUTED, lw=0.6, ls=':')
        ax.set_aspect('equal'); ax.axis('off'); ax.set_title(t, fontsize=8, loc='left')
        c = cases[nm]
        ax.text(0.0, -0.06, f"EDA {c['tot']:.2f} kcal/mol · N r {c['r']:.1f} Å, h {c['h']:.1f} Å · "
                f"{c['cc']} C···C", transform=ax.transAxes, fontsize=6.4, color=INK2)
    fig.text(0.01, 0.0, 'Side views (anchor plane dotted). Green: C8 patch; orange: reagent; blue: N.',
             fontsize=6.8, color=INK2)
    fig.tight_layout()
    save(fig, 'f6_structures.png')


# ---------------------------------------------------------------- 7 retention fit
def fig_fit():
    import csv
    rows = [r for r in csv.reader(l for l in open(ROOT / 'examples/hplc/donegan2022_table4.csv') if not l.startswith('#'))]
    hdr, rows = rows[0], rows[1:]
    D = {r[1]: dict(zip(hdr, r)) for r in rows}
    reg = [r for r in ORDER if D[r]['pct_acn_20T']]
    y = np.array([float(D[r]['pct_acn_20T']) for r in reg])
    key = {v['model'] if 'model' in v else None: k for k, v in []}
    tag = {D[r]['model']: D[r]['reagent'] for r in D}
    xs = [
        ('logD (ionized amine)', np.array([float(D[r]['logD']) for r in reg]), 'baseline/logD', 'logD'),
        ('EDA total, all placements', np.array([fit['reagents'][tag[r]]['all']['tot/boltzmann'] for r in reg]),
         'all/tot/boltzmann', 'kcal/mol'),
        ('EDA total, hollow excluded', np.array([fit['reagents'][tag[r]]['off_centre']['tot/boltzmann'] for r in reg]),
         'off_centre/tot/boltzmann', 'kcal/mol'),
    ]
    single = fit['fits']['pct_acn_20T']['single']
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    for k, (ax, (t, x, fk, unit)) in enumerate(zip(axs, xs)):
        st = single.get(fk)
        if st is None:
            r_ = np.corrcoef(x, y)[0, 1]; loo = None
        else:
            r_, loo = st['pearson_r'], st['loo_rmse']
        a, b = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 10)
        ax.plot(xx, a * xx + b, color=MUTED, lw=1.2, ls='--', zorder=1)
        for xi, yi, r in zip(x, y, reg):
            ax.scatter(xi, yi, s=40, color=RCOL[r], ec='white', lw=1, zorder=3)
            ax.annotate(LAB[r], (xi, yi), xytext=(4, 3), textcoords='offset points', fontsize=6.6, color=INK)
        ax.set_title(['a', 'b', 'c'][k] + '  ' + t, loc='left', fontsize=8.2)
        ax.set_xlabel(unit)
        lbl = f'r = {r_:.2f}' + (f'\nLOO RMSE {loo:.1f} %' if loo else '')
        ax.text(0.04 if k == 0 else 0.97, 0.97 if k == 0 else 0.04, lbl, transform=ax.transAxes, fontsize=7,
                color=INK, ha='left' if k == 0 else 'right', va='top' if k == 0 else 'bottom')
        if unit == 'kcal/mol':
            ax.invert_xaxis()
    axs[0].set_ylabel('% ACN at elution, dT20')
    fig.text(0.01, -0.04, 'Mean-only model LOO RMSE 11.2 % ACN. EDA axes inverted (stronger binding to the right). '
             'n = 7 (TBA not eluted).', fontsize=6.8, color=INK2)
    fig.tight_layout()
    save(fig, 'f7_retention_fit.png')


# ---------------------------------------------------------------- 8 basis + rho history
def fig_basis_rho():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[0.85, 1.4]))
    ax = axs[0]
    data = [('teah', -5.33, -7.60), ('dipeah', -5.49, -7.35), ('dibutylammonium', -8.85, -12.19),
            ('hexylammonium', -5.97, -8.30), ('tributylammonium', -7.14, -10.06)]
    for r, a, b in data:
        ax.plot([0, 1], [a, b], color=RCOL[r], lw=2, marker='o', ms=6, mec='white')
        ax.text(1.06, b, f'{LAB[r]} {b:.1f}', va='center', fontsize=7, color=INK)
    ax.set_xticks([0, 1], ['def2-SVP', 'aug-cc-pVDZ']); ax.set_xlim(-0.15, 1.5)
    ax.set_ylabel('EDA total (kcal/mol)'); ax.grid(axis='x', visible=False)
    ax.set_title('a  Basis-set check (ρ = 0.9)', loc='left', fontsize=9)
    ax = axs[1]
    hist = [
        ('10-26', 'single C8, bridged T···S', 5, 0.8),
        ('10-27', '3-chain patch, 2 placements', 5, 0.3),
        ('10-28', '+3 placements (5 each)', 5, 0.7),
        ('11-10', 'Bondi ×1.4 radii', 5, 0.5),
        ('11-11', 'aug-cc-pVDZ', 5, 0.6),
        ('11-11', '+DEA, BA, TPA', 8, 0.45),
        ('11-11', 'n = 8, hollow excluded*', 8, 0.83),
    ]
    yv = np.arange(len(hist))[::-1]
    for yy, (d, t, n, rho) in zip(yv, hist):
        c = CAT[0] if n == 5 else CAT[1]
        ax.plot([0, rho], [yy, yy], color=c, lw=2, solid_capstyle='round')
        ax.scatter(rho, yy, s=40, color=c, ec='white', lw=1, zorder=3)
        ax.text(rho + 0.03, yy, f'{rho:.2f}', va='center', fontsize=7)
    ax.set_yticks(yv, [f'{d}  {t}' for d, t, n, r in hist], fontsize=7)
    ax.set_xlim(0, 1.05); ax.set_xlabel('Spearman ρ, EDA total vs published reagent order')
    ax.axvline(0.738, color=MUTED, lw=0.8, ls=':')
    ax.text(0.745, -0.5, 'p = .05 (n = 8)', fontsize=6, color=INK2, ha='left')
    ax.set_ylim(-0.7, 7.9)
    ax.scatter([], [], color=CAT[0], label='n = 5 reagents'); ax.scatter([], [], color=CAT[1], label='n = 8 reagents')
    ax.legend(fontsize=6.6, loc='upper left', ncol=2, bbox_to_anchor=(0.0, 1.0))
    ax.set_title('b  Agreement with experiment over time', loc='left', fontsize=9)
    fig.tight_layout()
    save(fig, 'f8_basis_rho.png')


# ---------------------------------------------------------------- 9 timeline
def fig_timeline():
    ev = [
        ('2026-10-11', 'CPU/GPU optimizer, GKS-EDA (gas phase), validation suite', 0),
        ('2026-10-18', 'Psi4 cross-check, tryptamine benchmark, GKS-EDA(sol), scaling', 0),
        ('2026-10-19', 'HPLC 7.1/7.2: ion pairs, residue/reagent on C8', 1),
        ('2026-10-20', 'Silanol-bound ternaries, three-body EDA', 1),
        ('2026-10-21', 'Reagent series (5 reagents)', 1),
        ('2026-10-22', 'Bridged ternaries', 1),
        ('2026-10-26', 'SMD CDS nonpolar term', 2),
        ('2026-10-27', 'Silanol-free 3-chain C8 patch', 2),
        ('2026-10-28', 'Placement sampling; 7-chain brush GPU OOM', 2),
        ('2026-11-09', 'GPU cache fix; 7-chain retest OOM → dropped', 3),
        ('2026-11-10', 'Solvent sensitivity (radii, ε)', 3),
        ('2026-11-11', 'aug-cc-pVDZ, n = 8 reagents, retention fit', 3),
    ]
    import datetime as dt
    t0 = dt.date(2026, 10, 9)
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    lanes = ['Method & validation', 'HPLC models (silanol C8)', 'Silanol-free patch', 'Sensitivity & fit']
    for k, (d, t, lane) in enumerate(ev):
        day = (dt.date.fromisoformat(d) - t0).days
        y = 3 - lane
        ax.scatter(day, y, s=46, color=CAT[lane], ec='white', lw=1, zorder=3)
        dy = 0.22 if k % 2 == 0 else -0.3
        ax.annotate(t, (day, y), xytext=(day + 0.2, y + dy), fontsize=6.1, color=INK, va='center')
    ax.set_yticks(range(4), lanes[::-1], fontsize=7.5)
    ticks = [0, 7, 14, 21, 28, 33]
    ax.set_xticks(ticks, [(t0 + dt.timedelta(days=t)).strftime('%b %d') for t in ticks])
    ax.set_xlim(-1, 45); ax.set_ylim(-0.7, 3.6); ax.grid(axis='y', visible=False)
    ax.set_title('Project timeline (from handoffs and commits)', loc='left')
    save(fig, 'f9_timeline.png')


for f in [fig_workflow, fig_validation, fig_scaling, fig_components, fig_hollow, fig_structures, fig_fit,
          fig_basis_rho]:
    f(); print('ok', f.__name__)
