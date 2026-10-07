#figures and tables for the real-data validation (results/real)
import os
import sys
import glob
import numpy as np
import pandas as pd
pd.set_option('display.width', 250)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from plot_results import LABELS, INK, MUTED, FIG, TAB

METHODS = ['GAHIC3D', 'pastis-pm2', 'pastis-pm1', 'pastis-mds', 'lordg', 'shrec3d',
           'minimds', 'hsa1', 'chromsde', 'chromosome3d', 'gem']
LAB = dict(LABELS, minimds='miniMDS', hsa1='HSA', chromsde='ChromSDE', chromosome3d='Chromosome3D', gem='GEM')
# fixed categorical order: slots 1-8 of the validated palette, then neutral greys for 9-11 (direct labels on axis)
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948',
       '#6b6b66', '#8a8a85', '#a9a9a3']


def load():
    files = glob.glob(os.path.join(ROOT, 'results', 'real', 'evals', '*.csv'))
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d['CHROMOSOME'] = d['CHROMOSOME'].astype(str)
    d['RES'] = d['RESOLUTION'].map({1000000: '1 Mb', 400000: '400 kb'})
    return d


def boxes(ax, d, metric, methods):
    for k, m in enumerate(methods):
        v = d[d.METHOD == m][metric].dropna().values
        if len(v) == 0:
            continue
        c = COL[METHODS.index(m)]
        bp = ax.boxplot(v, positions=[k], widths=0.55, patch_artist=True, showfliers=False,
                        medianprops=dict(color=INK, linewidth=1.2), whiskerprops=dict(color=MUTED),
                        capprops=dict(color=MUTED))
        bp['boxes'][0].set(facecolor=c + '55', edgecolor=c, linewidth=1.2)
        ax.scatter(k + np.random.default_rng(k).uniform(-0.18, 0.18, len(v)), v, s=7, color=c,
                   edgecolor='white', linewidth=0.3, zorder=3)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels([LAB[m] for m in methods], rotation=45, ha='right')


def paired(d, metric, better='higher', key=('CELL', 'CHROMOSOME')):
    rows = []
    for res, sub in d.groupby('RES'):
        g = sub[sub.METHOD == 'GAHIC3D'].set_index(list(key)).sort_index()
        for m in METHODS[1:]:
            o = sub[sub.METHOD == m].set_index(list(key)).sort_index()
            k = g.index.intersection(o.index)
            if len(k) < 5:
                continue
            a, b = g.loc[k, metric], o.loc[k, metric]
            ok = ~(a.isna() | b.isna())
            a, b = a[ok], b[ok]
            win = (a > b).sum() if better == 'higher' else (a < b).sum()
            p = stats.wilcoxon(a, b).pvalue if len(a) >= 5 else np.nan
            rows.append([res, LAB[m], len(a), round(a.median(), 3), round(b.median(), 3), int(win), '{0:.2g}'.format(p)])
    return pd.DataFrame(rows, columns=['RES', 'vs', 'n', 'median_GAHIC3D', 'median_other', 'GAHIC3D_better', 'wilcoxon_p'])


def main():
    d = load()
    d.to_csv(os.path.join(ROOT, 'results', 'real', 'real_all.csv'), index=False)
    A, B, C = d[d.ANALYSIS == 'A'], d[d.ANALYSIS == 'B'], d[d.ANALYSIS == 'C']
    ceil = pd.read_csv(os.path.join(ROOT, 'results', 'real', 'model_agreement.csv'))
    betw = pd.read_csv(os.path.join(ROOT, 'results', 'real', 'between_cell_agreement.csv'))

    # Fig 5: three panels per resolution
    for res in ['1 Mb', '400 kb']:
        ms = [m for m in METHODS if m in set(d[d.RES == res].METHOD)]
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        boxes(axes[0], A[A.RES == res], 'PCoeff', ms)
        r = 1000000 if res == '1 Mb' else 400000
        axes[0].axhline(betw[betw.RESOLUTION == r].PCoeff.median(), color=MUTED, linestyle='--', linewidth=1)
        axes[0].set_ylabel('Pearson r of pairwise distances')
        axes[0].set_title('A  Agreement with the published structure', fontsize=9, color=INK, loc='left')
        boxes(axes[1], B[B.RES == res], 'AUC', [m for m in ms if m in set(B.METHOD)])
        axes[1].axhline(0.5, color=MUTED, linestyle='--', linewidth=1)
        axes[1].set_ylabel('AUC (held-out contacts)')
        axes[1].set_title('B  Prediction of held-out contacts (10%)', fontsize=9, color=INK, loc='left')
        boxes(axes[2], C[C.RES == res], 'SCoeff_mean', ms)
        axes[2].set_ylabel('Spearman ρ with mean single-cell distances')
        axes[2].set_title('C  Population Hi-C', fontsize=9, color=INK, loc='left')
        fig.tight_layout()
        tag = '1mb' if res == '1 Mb' else '400kb'
        fig.savefig(os.path.join(FIG, 'fig5_realdata_{0}.png'.format(tag)), dpi=200)
        fig.savefig(os.path.join(FIG, 'fig5_realdata_{0}.pdf'.format(tag)))
        plt.close(fig)

    tA = A.groupby(['RES', 'METHOD'])[['RMSD', 'PCoeff', 'SCoeff']].median().round(3)
    tB = B.groupby(['RES', 'METHOD'])[['AUC']].agg(['median', 'mean']).round(3)
    tC = C.groupby(['RES', 'METHOD'])[['SCoeff_mean', 'PCoeff_mean', 'SCoeff_cells']].median().round(3)
    tA.to_csv(os.path.join(TAB, 'table6a_real_singlecell.csv'))
    tB.to_csv(os.path.join(TAB, 'table6b_real_heldout.csv'))
    tC.to_csv(os.path.join(TAB, 'table6c_real_population.csv'))
    pA, pB = paired(A, 'PCoeff'), paired(B, 'AUC')
    pC = paired(C, 'SCoeff_mean', key=('CHROMOSOME',)) if len(C) else None
    pA.to_csv(os.path.join(TAB, 'table6d_real_singlecell_tests.csv'), index=False)
    pB.to_csv(os.path.join(TAB, 'table6e_real_heldout_tests.csv'), index=False)
    ce = ceil.groupby('RESOLUTION')[['RMSD_models', 'PCoeff_models']].median()
    be = betw.groupby('RESOLUTION')[['RMSD', 'PCoeff']].median()
    print(tA, tB, tC, pA, pB, pC, ce, be, sep='\n\n')


if __name__ == '__main__':
    main()
