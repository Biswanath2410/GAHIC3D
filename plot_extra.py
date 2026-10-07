#figures and tables for run_extra.py
#  fig6_depth        accuracy vs sequencing depth (1 Mb all chromosomes; 400 kb chr19)
#  figS_sensitivity  GAHIC3D settings (chr19, 400 kb, full and 1 % depth)
#  figS_convergence  best deviance vs number of local relaxations, GAHIC3D vs multi-start
#usage: python plot_extra.py
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from plot_results import COLORS, METHODS, LABELS, MUTED, INK, FIG, ROOT

TAB = os.path.join(ROOT, 'results', 'tables')
DEPTHS = [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]
DEPTH_METHODS = ['GAHIC3D', 'pastis-pm2', 'pastis-pm1', 'pastis-mds', 'lordg', 'shrec3d', 'minimds', 'chromsde']
SENS = [('GAHIC3D', 'Default settings'), ('GAHIC3D-pop5', 'Population size 5'), ('GAHIC3D-pop20', 'Population size 20'),
        ('GAHIC3D-relax20', '20 local optimizations'), ('GAHIC3D-relax80', '80 local optimizations'),
        ('GAHIC3D-relax120', '120 local optimizations'),
        ('GAHIC3D-pivot0', 'No pivot mutation'), ('GAHIC3D-pivot0.6', 'Pivot probability 0.6'),
        ('GAHIC3D-iter100', '100 L-BFGS-B iterations'), ('GAHIC3D-iter400', '400 L-BFGS-B iterations'),
        ('GAHIC3D-a1.5', 'Exponent a = 1.5'), ('GAHIC3D-a2.5', 'Exponent a = 2.5'), ('GAHIC3D-a3', 'Exponent a = 3')]


def load():
    d = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(ROOT, 'results', 'evals', '*.csv'))], ignore_index=True)
    d['CHROMOSOME'] = d['CHROMOSOME'].astype(str)
    d = d[~d.TAG.str.contains('nb')]
    d['COV'] = d['TAG'].str.extract(r'cov([0-9.]+)')[0].astype(float)
    return d


def fig_depth(d):
    panels = [(1000000, None, '1 Mb, chromosomes 1, 9, 19 and X (32 structures)'), (400000, '19', '400 kb, chromosome 19 (8 structures)')]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.4), sharex=True)
    rows = []
    for j, (r, ch, title) in enumerate(panels):
        s = d[(d.RESOLUTION == r) & d.COV.isin(DEPTHS) & d.METHOD.isin(DEPTH_METHODS)]
        if ch is not None:
            s = s[s.CHROMOSOME == ch]
        for i, (col, ylab) in enumerate([('RMSD', 'RMSD to ground-truth structure'), ('PCoeff', 'Pearson r of pairwise distances')]):
            ax = axes[i, j]
            for m in DEPTH_METHODS:
                g = s[s.METHOD == m].groupby('COV')[col]
                if g.ngroups == 0:
                    continue
                med, q1, q3 = g.median(), g.quantile(0.25), g.quantile(0.75)
                x = med.index.values * 100
                c = COLORS[METHODS.index(m)]
                lw = 2.4 if m == 'GAHIC3D' else 1.4
                ax.plot(x, med.values, '-o', color=c, lw=lw, ms=4, label=LABELS[m], zorder=3 if m == 'GAHIC3D' else 2)
                if m == 'GAHIC3D':
                    ax.fill_between(x, q1.values, q3.values, color=c, alpha=0.15, lw=0)
                if col == 'RMSD':
                    for cov in med.index:
                        rows.append([r, ch or 'all', m, cov, med[cov], q1[cov], q3[cov], g.count()[cov]])
            ax.set_xscale('log')
            ax.invert_xaxis()
            ax.set_xticks([100, 30, 10, 3, 1, 0.3])
            ax.set_xticklabels(['100', '30', '10', '3', '1', '0.3'])
            ax.grid(alpha=0.3, lw=0.5)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(title, fontsize=9, color=INK)
            if i == 1:
                ax.set_xlabel('Sequencing depth (% of full depth)')
    axes[0, 0].legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'fig6_depth.' + ext), dpi=200)
    plt.close(fig)
    t = pd.DataFrame(rows, columns=['RESOLUTION', 'CHROMOSOMES', 'METHOD', 'DEPTH', 'RMSD_median', 'RMSD_q1', 'RMSD_q3', 'n'])
    t.round(3).to_csv(os.path.join(TAB, 'tableS4_depth_titration.csv'), index=False)
    # paired tests GAHIC3D vs each method at each depth (1 Mb, all chromosomes)
    out = []
    for r, ch, _ in panels:
        s = d[(d.RESOLUTION == r) & d.COV.isin(DEPTHS)]
        if ch is not None:
            s = s[s.CHROMOSOME == ch]
        p = s.pivot_table(index=['COV', 'CELL', 'CHROMOSOME'], columns='METHOD', values='RMSD')
        for cov in DEPTHS:
            if cov not in p.index.get_level_values(0):
                continue
            q = p.loc[cov]
            for m in DEPTH_METHODS[1:]:
                if m not in q or 'GAHIC3D' not in q:
                    continue
                a = q[['GAHIC3D', m]].dropna()
                if len(a) < 5:
                    continue
                out.append([r, cov, m, len(a), a.GAHIC3D.median(), a[m].median(), int((a.GAHIC3D < a[m]).sum()),
                            stats.wilcoxon(a.GAHIC3D, a[m]).pvalue])
    pd.DataFrame(out, columns=['RESOLUTION', 'DEPTH', 'METHOD', 'n', 'GAHIC3D', 'other', 'GAHIC3D_better', 'p']).to_csv(
        os.path.join(TAB, 'tableS4b_depth_tests.csv'), index=False)


def fig_sensitivity(d):
    s = d[(d.RESOLUTION == 400000) & (d.CHROMOSOME == '19') & d.COV.isin([1.0, 0.01])]
    names = [m for m, _ in SENS if m in set(s.METHOD)]
    labels = dict(SENS)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    rows = []
    for ax, cov, title in zip(axes, [1.0, 0.01], ['full depth', '1% depth']):
        p = s[s.COV == cov].pivot_table(index='CELL', columns='METHOD', values='RMSD')
        for k, m in enumerate(names):
            if m not in p:
                continue
            v = p[m].dropna()
            y = len(names) - 1 - k
            ax.scatter(v.values, np.full(len(v), y) + np.linspace(-0.15, 0.15, len(v)), s=12,
                       color=COLORS[0] if m == 'GAHIC3D' else MUTED, zorder=2)
            ax.plot([v.median()] * 2, [y - 0.3, y + 0.3], color=INK, lw=1.5, zorder=3)
            if m != 'GAHIC3D':
                a = p[['GAHIC3D', m]].dropna()
                pv = stats.wilcoxon(a.GAHIC3D, a[m]).pvalue if len(a) >= 5 else np.nan
                rows.append([cov, m, labels[m], len(a), a.GAHIC3D.median(), a[m].median(), pv])
            else:
                rows.append([cov, m, labels[m], len(v), v.median(), v.median(), np.nan])
        ax.axvline(p['GAHIC3D'].median(), color=COLORS[0], lw=0.8, ls='--', zorder=1)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([labels[m] for m in names][::-1])
        ax.set_xlabel('RMSD to ground-truth structure')
        ax.set_title('Chromosome 19, 400 kb, ' + title, fontsize=10, color=INK)
        ax.grid(axis='x', alpha=0.3, lw=0.5)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'figS_sensitivity.' + ext), dpi=200)
    plt.close(fig)
    pd.DataFrame(rows, columns=['DEPTH', 'METHOD', 'setting', 'n', 'default_median', 'setting_median', 'p']).round(4).to_csv(
        os.path.join(TAB, 'tableS5_sensitivity.csv'), index=False)


def fig_convergence():
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=False)
    rows = []
    for ax, tag, title in zip(axes, ['cov1.0', 'cov0.01'], ['full depth', '1% depth']):
        curves = {'GAHIC3D': [], 'MultiStart': []}
        for ce in range(1, 9):
            dd = os.path.join(ROOT, 'results', 'structures', tag, '{0}_19_400000'.format(ce))
            fg, fm = os.path.join(dd, 'GAHIC3D-relax120_trace.txt'), os.path.join(dd, 'MultiStart-relax120_trace.txt')
            if not (os.path.exists(fg) and os.path.exists(fm)):
                continue
            g, m = np.loadtxt(fg), np.loadtxt(fm)
            P = 10
            xg = g[:, 0] + P                          # relaxations used: P initial + one per offspring
            yg, ym = g[:, 1], m[:, 1]                  # best deviance so far
            ref = min(yg.min(), ym.min())
            curves['GAHIC3D'].append(np.interp(np.arange(P, 121), xg, yg - ref))
            curves['MultiStart'].append((ym - ref)[P - 1:120])
            rows.append([tag, ce, yg[-1], ym[-1], np.interp(40, xg, yg), ym[39]])
        x = np.arange(10, 121)
        for name, c, lab in [('GAHIC3D', COLORS[0], 'GAHIC3D'), ('MultiStart', MUTED, 'Multi-start')]:
            if not curves[name]:
                continue
            a = np.array(curves[name]) + 1e-6
            ax.plot(x, np.median(a, 0), color=c, lw=2, label=lab)
            ax.fill_between(x, np.percentile(a, 25, 0), np.percentile(a, 75, 0), color=c, alpha=0.15, lw=0)
        ax.axvline(40, color=INK, lw=0.6, ls=':')
        ax.text(41, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else 0, ' default (40)', fontsize=7, color=INK, va='top')
        ax.set_yscale('log')
        ax.set_xlabel('Number of local optimizations')
        ax.set_title('Chromosome 19, 400 kb, ' + title, fontsize=10, color=INK)
        ax.grid(alpha=0.3, lw=0.5)
    axes[0].set_ylabel('Excess deviance over best structure found')
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'figS_convergence.' + ext), dpi=200)
    plt.close(fig)
    pd.DataFrame(rows, columns=['TAG', 'CELL', 'GAHIC3D_best120', 'MultiStart_best120', 'GAHIC3D_best40', 'MultiStart_best40']).to_csv(
        os.path.join(TAB, 'tableS6_convergence.csv'), index=False)


def main():
    d = load()
    fig_depth(d)
    fig_sensitivity(d)
    fig_convergence()


if __name__ == '__main__':
    main()
