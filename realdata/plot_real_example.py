#Example reconstruction from experimental single-cell Hi-C (GSE80280): published structure vs reconstructions.
#Row 1: pairwise-distance matrices (published structure and each reconstruction; Pearson r to the published one).
#Row 2: reconstructions superposed on the published structure (Kabsch, reflection allowed), drawn as points coloured
#       by position along the chromosome (no connecting lines, so bins with few contacts appear as isolated points).
#usage: python realdata/plot_real_example.py [cell] [chrom] [res]     (default: cell 7, chr9, 400 kb; Fig. 6)
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from GAHIC3D import getTransformation
from plot_results import COLORS, METHODS, LABELS, MUTED, INK, FIG

DATA = os.path.join(ROOT, 'data', 'real')
POSCMAP = 'plasma'
OUT = os.path.join(ROOT, 'results', 'real', 'structures', 'A')


def load(f):
    a = np.loadtxt(f)
    return a[:, 0].astype(int), a[:, 1:4]


def ranked(D):
    """pairwise distances as percentile ranks (0 = closest pair, 1 = farthest), so that patterns are comparable"""
    iu = np.triu_indices(len(D), 1)
    R = np.zeros_like(D)
    R[iu] = rankdata(D[iu]) / len(iu[0])
    return R + R.T


def main():
    ce = sys.argv[1] if len(sys.argv) > 1 else '7'
    ch = sys.argv[2] if len(sys.argv) > 2 else '9'
    r = sys.argv[3] if len(sys.argv) > 3 else '400000'
    show = ['GAHIC3D', 'shrec3d', 'pastis-pm1']
    ri, R = load(os.path.join(DATA, 'ref', 'ref_Cell_{0}_{1}_{2}_model1.txt'.format(ce, ch, r)))
    ests = {m: load(os.path.join(OUT, '{0}_{1}_{2}'.format(ce, ch, r), m + '_struct.txt')) for m in show}
    common = ri
    for ei, _ in ests.values():
        common = np.intersect1d(common, ei)
    Rc = R[np.searchsorted(ri, common)]
    DR = squareform(pdist(Rc))
    iu = np.triu_indices(len(common), 1)

    fig = plt.figure(figsize=(11, 5.6))
    ax = fig.add_subplot(2, 4, 1)
    ax.imshow(ranked(DR), cmap='viridis_r', vmin=0, vmax=1)
    ax.set_title('Published structure', fontsize=9, color=INK)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_ylabel('pairwise distance (rank)', fontsize=9)
    pos = np.linspace(0, 1, len(common))                      # colour = position along the chromosome
    ax3 = fig.add_subplot(2, 4, 5, projection='3d')
    sc = ax3.scatter(*Rc.T, c=pos, cmap=POSCMAP, s=4, depthshade=False)
    ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
    ax3.set_title('Published structure', fontsize=8, color=INK)
    for k, m in enumerate(show):
        ei, E = ests[m]
        Ec = E[np.searchsorted(ei, common)]
        X, Y = getTransformation(Ec, Rc, reflection=True)     # X: prediction superposed on Y (published, scaled)
        D = squareform(pdist(X)) * (DR[iu].mean() / pdist(X).mean())
        rr = np.corrcoef(D[iu], DR[iu])[0, 1]
        ax = fig.add_subplot(2, 4, k + 2)
        im = ax.imshow(ranked(D), cmap='viridis_r', vmin=0, vmax=1)
        ax.set_title('{0}  (r = {1:.2f})'.format(LABELS[m], rr), fontsize=9, color=INK)
        ax.set_xticks([]); ax.set_yticks([])
        ax3 = fig.add_subplot(2, 4, k + 6, projection='3d')
        ax3.scatter(*Y.T, color=MUTED, s=2, alpha=0.25, depthshade=False, label='published')
        ax3.scatter(*X.T, c=pos, cmap=POSCMAP, s=4, depthshade=False, label=LABELS[m])
        ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
        ax3.set_title('{0} on published (grey)'.format(LABELS[m]), fontsize=8, color=INK)
    cb = fig.colorbar(im, cax=fig.add_axes([0.93, 0.56, 0.01, 0.3]))
    cb.set_ticks([0, 1]); cb.set_ticklabels(['close', 'far'])
    cb2 = fig.colorbar(sc, cax=fig.add_axes([0.93, 0.08, 0.01, 0.3]))
    cb2.set_ticks([0, 1]); cb2.set_ticklabels(['start', 'end'])
    cb2.set_label('position on chromosome', fontsize=8)
    fig.suptitle('Cell {0}, chromosome {1}, {2} kb: {3} bins with contacts'.format(ce, ch, int(r) // 1000, len(common)),
                 fontsize=10, color=INK)
    fig.tight_layout(rect=[0, 0, 0.9, 1])
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'figS_real_example_{0}_{1}_{2}.{3}'.format(ce, ch, r, ext)), dpi=200)


if __name__ == '__main__':
    main()
