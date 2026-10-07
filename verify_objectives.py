"""Sanity check for the benchmark: is the failure of other methods on sparse data a property of their
objective, or a bug in how we ran them?
For each simulated input we score the TRUE structure and each reconstruction under
  (i)  the Poisson likelihood on all pairs (zeros included; GAHIC3D)
  (ii) the same likelihood on non-zero pairs only (the PASTIS-PM1 model)
and compare with random-walk structures (RMSD of a structure that carries no information).
If a method's wrong structure scores better than the truth under the objective it optimizes,
the failure lies in the objective, not in the implementation.
output: results/tables/verify_objectives.csv
usage:  python verify_objectives.py"""
import os
import numpy as np
import pandas as pd
from GAHIC3D import Objective, makeSymmetric, remove_nan_col, findEuclideanDist, create_pop3
from simulateHiC import load_groundtruth
from RMSD import compare
from iced import normalization

ROOT = os.path.dirname(os.path.abspath(__file__))
METHODS = ['GAHIC3D', 'pastis-pm1', 'pastis-pm2', 'pastis-mds', 'shrec3d', 'lordg', 'minimds', 'chromsde']
CONDS = [('19', '400000', '0.01'), ('19', '400000', '1.0'), ('19', '1000000', '0.01'), ('X', '1000000', '0.01'),
         ('19', '400000', '0.1'), ('1', '1000000', '0.003')]


def score(obj, X):
    E = findEuclideanDist(X[None], obj.I, obj.J)
    return obj.error(E)[0]


def main():
    rng = np.random.default_rng(0)
    rows = []
    for ch, r, cov in CONDS:
        for ce in range(1, 9):
            m = np.loadtxt(os.path.join(ROOT, 'data', 'simulated', 'alpha2.0_cov' + cov,
                                        'simHiC_Cell_{0}_{1}_{2}_matrix.txt'.format(ce, ch, r)))
            m, idx = remove_nan_col(makeSymmetric(m))
            C = normalization.ICE_normalization(m)
            gi, g = load_groundtruth(os.path.join(ROOT, 'data', 'groundtruth', 'groundtruth_{0}_{1}_{2}.txt'.format(ce, ch, r)), ch, r)
            structs = {'truth': (gi, g)}
            for meth in METHODS:
                f = os.path.join(ROOT, 'results', 'structures', 'cov' + cov, '{0}_{1}_{2}'.format(ce, ch, r), meth + '_struct.txt')
                if os.path.exists(f):
                    a = np.loadtxt(f)
                    structs[meth] = (a[:, 0].astype(int), a[:, 1:4])
            for k in range(3):
                structs['random walk %d' % k] = (idx, create_pop3(1, len(idx), rng)[0])
            common = idx
            for si, _ in structs.values():
                common = np.intersect1d(common, si)
            sel = np.searchsorted(idx, common)
            Cs = C[np.ix_(sel, sel)]
            objs = {'all pairs': Objective(Cs, 2.0, 'poisson'), 'non-zero pairs': Objective(Cs, 2.0, 'poissonnz')}
            for name, (si, X) in structs.items():
                Xc = X[np.searchsorted(si, common)]
                res = compare(common, Xc, gi, g)
                rows.append([ch, r, cov, ce, name, len(common), score(objs['all pairs'], Xc), score(objs['non-zero pairs'], Xc),
                             res['RMSD'], res['PCoeff'], (Cs[np.triu_indices(len(common), 1)] == 0).mean()])
    t = pd.DataFrame(rows, columns=['CHROMOSOME', 'RESOLUTION', 'DEPTH', 'CELL', 'STRUCTURE', 'nbins',
                                    'deviance_all_pairs', 'deviance_nonzero_pairs', 'RMSD', 'PCoeff', 'fraction_zero'])
    t.to_csv(os.path.join(ROOT, 'results', 'tables', 'verify_objectives.csv'), index=False)
    t['STRUCTURE'] = t.STRUCTURE.str.replace(r'random walk \d', 'random walk', regex=True)
    print(t.groupby(['CHROMOSOME', 'RESOLUTION', 'DEPTH', 'STRUCTURE'])[
        ['deviance_all_pairs', 'deviance_nonzero_pairs', 'RMSD', 'PCoeff', 'fraction_zero']].median().round(3).to_string())
    # how often does a reconstruction beat the truth under each objective?
    w = t.pivot_table(index=['CHROMOSOME', 'RESOLUTION', 'DEPTH', 'CELL'], columns='STRUCTURE', values='deviance_nonzero_pairs', aggfunc='min')
    print('\nnon-zero objective: structures scoring better than the truth (cells out of 8)')
    print(pd.DataFrame({m: (w[m] < w['truth']).groupby(level=[0, 1, 2]).sum() for m in METHODS + ['random walk'] if m in w}).to_string())
    w = t.pivot_table(index=['CHROMOSOME', 'RESOLUTION', 'DEPTH', 'CELL'], columns='STRUCTURE', values='deviance_all_pairs', aggfunc='min')
    print('\nall-pairs objective: structures scoring better than the truth (cells out of 8)')
    print(pd.DataFrame({m: (w[m] < w['truth']).groupby(level=[0, 1, 2]).sum() for m in METHODS + ['random walk'] if m in w}).to_string())


if __name__ == '__main__':
    main()
