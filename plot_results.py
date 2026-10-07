#figures and tables for the manuscript from results/benchmark_all.csv
#usage: python plot_results.py
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(ROOT, 'figures')
TAB = os.path.join(ROOT, 'results', 'tables')
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

# fixed categorical order (validated palette, slots 1-6); GAHIC3D always slot 1
METHODS = ['GAHIC3D', 'pastis-pm2', 'pastis-pm1', 'pastis-mds', 'lordg', 'shrec3d',
           'minimds', 'hsa1', 'chromsde', 'chromosome3d', 'gem']
LABELS = {'GAHIC3D': 'GAHIC3D', 'pastis-pm2': 'PASTIS-PM2', 'pastis-pm1': 'PASTIS-PM1',
          'pastis-mds': 'PASTIS-MDS', 'lordg': 'LorDG', 'shrec3d': 'ShRec3D',
          'MultiStart': 'Multi-start (40 local optimizations)', 'GA-poissonNZ': 'Poisson, non-zero pairs only',
          'GA-wstress': 'Relative stress', 'GA-stress': 'Unweighted stress',
          'GA-generational': 'Generational GA, no local optimization', 'GA-2022': 'GAHIC3D 2022',
          'minimds': 'miniMDS', 'hsa1': 'HSA', 'chromsde': 'ChromSDE', 'chromosome3d': 'Chromosome3D', 'gem': 'GEM'}
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948',
          '#6b6b66', '#8a8a85', '#a9a9a3']
INK = '#2b2b2b'
MUTED = '#8a8a85'

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e6e3',
                     'grid.linewidth': 0.6, 'axes.axisbelow': True, 'legend.frameon': False})


def load():
    files = glob.glob(os.path.join(ROOT, 'results', 'evals', '*.csv'))
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d['CHROMOSOME'] = d['CHROMOSOME'].astype(str)
    d['RES'] = d['RESOLUTION'].map({1000000: '1 Mb', 400000: '400 kb'})
    d['NB'] = d['TAG'].str.contains('nb')
    d['COV'] = d['TAG'].str.extract(r'cov([0-9.]+)')[0].astype(float)
    d['COVL'] = d['COV'].map({1.0: 'full depth', 0.01: '1% depth'})
    # depth titration, parameter sensitivity and convergence runs (run_extra.py) are plotted by plot_extra.py
    extra = d.METHOD.str.contains(r'-(?:pop|relax|pivot|iter|a)[0-9]') | d.COVL.isna()
    return d[~extra].copy()


def runtime_table(d):
    rows = []
    for f in glob.glob(os.path.join(ROOT, 'results', 'structures', '*', '*', '*_param.txt')):
        p = pd.read_csv(f, sep='\t')
        cond = f.split(os.sep)
        ce, ch, r = cond[-2].split('_')
        method = os.path.basename(f).replace('_param.txt', '')
        rows.append([method, cond[-3], ch, r, float(p['seconds'].iloc[0])])
    t = pd.DataFrame(rows, columns=['METHOD', 'TAG', 'CHROMOSOME', 'RESOLUTION', 'seconds'])
    # main benchmark only (the depth series and sensitivity runs of run_extra.py are excluded)
    t = t[~t.TAG.isin(['cov0.3', 'cov0.1', 'cov0.03', 'cov0.003'])
          & ~t.METHOD.str.contains(r'-(?:pop|relax|pivot|iter|a)[0-9]')]
    t.to_csv(os.path.join(TAB, 'runtime_raw.csv'), index=False)
    s = t.groupby(['METHOD', 'RESOLUTION'])['seconds'].median().unstack().round(1)
    s.to_csv(os.path.join(TAB, 'runtime_median_seconds.csv'))
    return s


def box_panel(ax, sub, metric):
    """one box + points per method; methods without results in this condition are left out"""
    ms = [m for m in METHODS if len(sub[sub.METHOD == m][metric].dropna())]
    for x, m in enumerate(ms):
        k = METHODS.index(m)
        v = sub[sub.METHOD == m][metric].dropna().values
        bp = ax.boxplot(v, positions=[x], widths=0.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color=INK, linewidth=1.2),
                        whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
        bp['boxes'][0].set(facecolor=COLORS[k] + '55', edgecolor=COLORS[k], linewidth=1.2)
        ax.scatter(x + np.random.default_rng(k).uniform(-0.18, 0.18, len(v)), v, s=9, color=COLORS[k],
                   edgecolor='white', linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(ms)))
    ax.set_xticklabels([LABELS[m] for m in ms], rotation=45, ha='right')
    ax.set_xlim(-0.6, len(ms) - 0.4)
    return ms


def fig_main(d):
    """Fig 2 / Fig S4: RMSD and distance correlation per method, by resolution x depth"""
    main = d[d.METHOD.isin(METHODS)]
    panels = [(r, c) for c in ['full depth', '1% depth'] for r in ['1 Mb', '400 kb']]
    for metric, ylab, fname in [('RMSD', 'RMSD to ground-truth structure', 'fig2_rmsd'),
                                ('PCoeff', 'Pearson r of pairwise distances', 'fig2_pearson')]:
        fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharey=True)
        for ax, (r, c) in zip(axes.ravel(), panels):
            box_panel(ax, main[(main.RES == r) & (main.COVL == c)], metric)
            ax.set_title('{0}, {1}'.format(r, c), color=INK, fontsize=10)
        axes[0, 0].set_ylabel(ylab)
        axes[1, 0].set_ylabel(ylab)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, fname + '.png'), dpi=200)
        fig.savefig(os.path.join(FIG, fname + '.pdf'))
        plt.close(fig)


def summary_tables(d):
    main = d[d.METHOD.isin(METHODS + ['MultiStart'])]
    t = main.groupby(['RES', 'COVL', 'METHOD'])[['RMSD', 'RMSE', 'PCoeff', 'SCoeff']].agg(['mean', 'std']).round(3)
    t.to_csv(os.path.join(TAB, 'table1_summary.csv'))
    # per chromosome median RMSD
    p = main.pivot_table(index=['RES', 'COVL', 'CHROMOSOME'], columns='METHOD', values='RMSD', aggfunc='median').round(2)
    p.to_csv(os.path.join(TAB, 'tableS1_rmsd_by_chromosome.csv'))
    # paired Wilcoxon: GAHIC3D vs each method (pairs = cell x chromosome)
    rows = []
    for (r, c), sub in main.groupby(['RES', 'COVL']):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index(['CELL', 'CHROMOSOME'])
        for m in METHODS[1:] + ['MultiStart']:
            o = sub[sub.METHOD == m].set_index(['CELL', 'CHROMOSOME'])
            k = g.index.intersection(o.index)
            if len(k) < 5:
                continue
            diff = g.loc[k, 'RMSD'] - o.loc[k, 'RMSD']
            w = stats.wilcoxon(g.loc[k, 'RMSD'], o.loc[k, 'RMSD'])
            rows.append([r, c, m, len(k), round(g.loc[k, 'RMSD'].median(), 3), round(o.loc[k, 'RMSD'].median(), 3),
                         round(diff.median(), 3), int((diff < 0).sum()), w.pvalue])
    w = pd.DataFrame(rows, columns=['RES', 'DEPTH', 'vs', 'n_pairs', 'median_RMSD_GAHIC3D', 'median_RMSD_other',
                                    'median_diff', 'GAHIC3D_better', 'wilcoxon_p'])
    w['wilcoxon_p'] = w['wilcoxon_p'].map(lambda x: '{0:.2g}'.format(x))
    w.to_csv(os.path.join(TAB, 'table2_paired_tests.csv'), index=False)
    return t, w


def fig_ablation(d):
    """Fig 3: what matters - objective vs optimiser (chr19, 400 kb)"""
    order = ['GAHIC3D', 'MultiStart', 'GA-generational', 'GA-poissonNZ', 'GA-wstress', 'GA-stress']
    sub = d[(d.CHROMOSOME == '19') & (d.RES == '400 kb') & d.METHOD.isin(order)]
    if sub.empty:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    for ax, c in zip(axes, ['full depth', '1% depth']):
        s = sub[sub.COVL == c]
        for k, m in enumerate(order):
            v = s[s.METHOD == m]['RMSD'].values
            if len(v) == 0:
                continue
            col = COLORS[0] if m == 'GAHIC3D' else MUTED
            ax.scatter(np.full(len(v), k) + np.random.default_rng(k).uniform(-0.15, 0.15, len(v)), v,
                       s=12, color=col, edgecolor='white', linewidth=0.4, zorder=3)
            ax.hlines(np.median(v), k - 0.3, k + 0.3, color=INK, linewidth=1.5, zorder=4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([LABELS[m] for m in order], rotation=40, ha='right')
        ax.set_title('Chromosome 19, 400 kb, ' + c, fontsize=10, color=INK)
    axes[0].set_ylabel('RMSD to ground-truth structure')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig3_ablation.png'), dpi=200)
    fig.savefig(os.path.join(FIG, 'fig3_ablation.pdf'))
    plt.close(fig)
    t = sub.groupby(['COVL', 'METHOD'])[['RMSD', 'PCoeff']].agg(['median', 'mean', 'std']).round(3)
    t.to_csv(os.path.join(TAB, 'table3_ablation.csv'))
    return t


def fig_chirality(d):
    """Fig S2: effect of allowing the mirror image in the superposition, one panel per method"""
    main = d[d.METHOD.isin(METHODS)]
    ms = [m for m in METHODS if (main.METHOD == m).any()]
    lim = [0, np.nanmax(main['RMSD_noreflect']) * 1.03]
    fig, axes = plt.subplots(3, 4, figsize=(10, 7.6), sharex=True, sharey=True)
    for ax, m in zip(axes.ravel(), ms):
        s = main[main.METHOD == m]
        k = METHODS.index(m)
        mir = (s['RMSD_noreflect'] > s['RMSD'] + 1e-6)
        ax.plot(lim, lim, color=MUTED, linewidth=0.8, zorder=1)
        ax.scatter(s['RMSD'], s['RMSD_noreflect'], s=7, color=COLORS[k], alpha=0.8, edgecolor='none', zorder=2)
        ax.set_title('{0}  ({1:.0f}% mirrored)'.format(LABELS[m], 100 * mir.mean()), fontsize=9, color=INK)
        ax.set_xlim(lim); ax.set_ylim(lim)
    for ax in axes.ravel()[len(ms):]:
        ax.axis('off')
    for ax in axes[:, 0]:
        ax.set_ylabel('RMSD, rotation only')
    for ax in axes[-1, :]:
        ax.set_xlabel('RMSD, reflection allowed')
    for ax in axes.ravel()[len(ms) - 4:len(ms)]:
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.set_xlabel('RMSD, reflection allowed')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'figS1_chirality.png'), dpi=200)
    fig.savefig(os.path.join(FIG, 'figS1_chirality.pdf'))
    plt.close(fig)
    frac = (main['RMSD_noreflect'] > main['RMSD'] + 1e-6).groupby(main.METHOD).mean().round(2)
    frac.to_csv(os.path.join(TAB, 'tableS2_fraction_mirrored.csv'))
    return frac


def fig_example(cell=1, ch='X', r='400000', tag='cov0.01'):
    """Fig 1C: ground truth vs reconstructions, superposed"""
    from simulateHiC import load_groundtruth
    from GAHIC3D import getTransformation
    gi, g = load_groundtruth(os.path.join(ROOT, 'data', 'groundtruth', 'groundtruth_{0}_{1}_{2}.txt'.format(cell, ch, r)), ch, r)
    show = ['GAHIC3D', 'pastis-pm1', 'shrec3d']
    fig = plt.figure(figsize=(10, 3.3))
    for k, m in enumerate(show):
        f = os.path.join(ROOT, 'results', 'structures', tag, '{0}_{1}_{2}'.format(cell, ch, r), m + '_struct.txt')
        if not os.path.exists(f):
            return
        e = np.loadtxt(f)
        common, ie, ig = np.intersect1d(e[:, 0].astype(int), gi, return_indices=True)
        X, Y = getTransformation(e[ie, 1:4], g[ig], reflection=True)
        ax = fig.add_subplot(1, 3, k + 1, projection='3d')
        ax.plot(*Y.T, color=MUTED, linewidth=1.0, label='ground truth')
        ax.plot(*X.T, color=COLORS[METHODS.index(m)], linewidth=1.2, label=LABELS[m])
        rmsd = np.sqrt(((X - Y) ** 2).sum(axis=1).mean())
        ax.set_title('{0}  (RMSD {1:.2f})'.format(LABELS[m], rmsd), fontsize=9, color=INK)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.legend(fontsize=7, loc='upper left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'fig1c_example_{0}_{1}_{2}_{3}.{4}'.format(cell, ch, r, tag, ext)), dpi=200)
    plt.close(fig)


def fig_runtime(s):
    if s is None or s.empty:
        return
    s = s.rename(columns=lambda c: int(c))
    s = s.reindex([m for m in METHODS if m in s.index])
    fig, ax = plt.subplots(figsize=(6, 3.6))
    x = np.arange(len(s))
    for j, (col, lab) in enumerate([(1000000, '1 Mb'), (400000, '400 kb')]):
        if col in s.columns:
            ax.bar(x + (j - 0.5) * 0.38, s[col].values, width=0.36, color=[COLORS[0], MUTED][j], label=lab)
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in s.index], rotation=40, ha='right')
    ax.set_ylabel('Median running time per structure (s)')
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis='y', alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'figS2_runtime.png'), dpi=200)
    fig.savefig(os.path.join(FIG, 'figS2_runtime.pdf'))
    plt.close(fig)


def fig_robust(d):
    """Fig 4: misspecified model (simulated alpha=-3 + negative-binomial noise; methods assume alpha=2)"""
    rb = d[d.NB & d.METHOD.isin(METHODS)]
    if rb.empty:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
    for ax, c in zip(axes, ['full depth', '1% depth']):
        box_panel(ax, rb[rb.COVL == c], 'RMSD')
        ax.set_title('1 Mb, ' + c, fontsize=10, color=INK)
    axes[0].set_ylabel('RMSD to ground-truth structure')
    fig.suptitle('Data simulated with exponent \u22123 and negative binomial noise; methods assume \u22122',
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fig4_misspecified.png'), dpi=200)
    fig.savefig(os.path.join(FIG, 'fig4_misspecified.pdf'))
    plt.close(fig)
    t = rb.groupby(['COVL', 'METHOD'])[['RMSD', 'PCoeff']].agg(['median', 'mean']).round(3)
    t.to_csv(os.path.join(TAB, 'table5_misspecified.csv'))
    rows = []
    for c, sub in rb.groupby('COVL'):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index(['CELL', 'CHROMOSOME']).sort_index()
        for m in METHODS[1:]:
            o = sub[sub.METHOD == m].set_index(['CELL', 'CHROMOSOME']).sort_index()
            k = g.index.intersection(o.index)
            rows.append([c, m, len(k), int((g.loc[k, 'RMSD'] < o.loc[k, 'RMSD']).sum()),
                         '{0:.2g}'.format(stats.wilcoxon(g.loc[k, 'RMSD'], o.loc[k, 'RMSD']).pvalue)])
    pd.DataFrame(rows, columns=['DEPTH', 'vs', 'n', 'GAHIC3D_better', 'wilcoxon_p']).to_csv(
        os.path.join(TAB, 'table5b_misspecified_tests.csv'), index=False)
    return t


def main():
    d_all = load()
    d_all.to_csv(os.path.join(ROOT, 'results', 'benchmark_all.csv'), index=False)
    fig_robust(d_all)
    d = d_all[~d_all.NB]
    fig_main(d)
    t, w = summary_tables(d)
    a = fig_ablation(d)
    frac = fig_chirality(d)
    s = runtime_table(d)
    fig_runtime(s)
    fig_example()
    print(w.to_string())
    print(a)
    print(frac)
    print(s)


if __name__ == '__main__':
    main()
