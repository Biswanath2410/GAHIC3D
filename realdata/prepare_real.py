#build the real-data inputs from GSE80280 (Stevens et al. 2017)
#   data/real/sc/       single-cell contact matrices (cells 1-8; chr1, 9, 19, X; 1 Mb, 400 kb)
#   data/real/holdout/  same, with 10 % of the contacted pairs (|i-j| >= 2) removed + the list of removed pairs
#   data/real/pop/      population Hi-C matrices
#   data/real/ref/      published structures binned to the same bins, all 10 models: bin x y z
#usage: python realdata/prepare_real.py
import os
import sys
import glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_gse80280 import read_pdb_models, bin_structure, read_sc_contacts, read_population
from simulateHiC import get_chrom_length

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, 'data', 'real', 'raw')
OUT = os.path.join(ROOT, 'data', 'real')
cells = [1, 2, 3, 4, 5, 6, 7, 8]
chroms = ['19', '9', 'X', '1']
resolutions = [1000000, 400000]
HOLDOUT = 0.10


def main():
    for d in ['sc', 'holdout', 'pop', 'ref']:
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    rng = np.random.default_rng(2017)
    stats = []
    for r in resolutions:
        for ch in chroms:
            nbins = int(get_chrom_length(ch) / r) + 1
            pop = read_population(os.path.join(OUT, 'population_hic_{0}.txt'.format(r)), ch, nbins)
            np.savetxt(os.path.join(OUT, 'pop', 'pop_{0}_{1}_matrix.txt'.format(ch, r)), pop, fmt='%d', delimiter='\t')
            for ce in cells:
                pairs = glob.glob(os.path.join(RAW, '*_Cell_{0}_contact_pairs.txt.gz'.format(ce)))[0]
                pdb = glob.glob(os.path.join(RAW, '*_Cell_{0}_genome_structure_model.pdb.gz'.format(ce)))[0]
                M = read_sc_contacts(pairs, ch, r, nbins)
                tag = 'Cell_{0}_{1}_{2}'.format(ce, ch, r)
                np.savetxt(os.path.join(OUT, 'sc', 'sc_' + tag + '_matrix.txt'), M, fmt='%d', delimiter='\t')
                # hold-out set
                I, J = np.triu_indices(nbins, k=2)
                hit = M[I, J] > 0
                cand = np.where(hit)[0]
                k = max(1, int(round(HOLDOUT * len(cand))))
                pick = rng.choice(cand, size=k, replace=False)
                T = M.copy()
                T[I[pick], J[pick]] = 0
                T[J[pick], I[pick]] = 0
                np.savetxt(os.path.join(OUT, 'holdout', 'train_' + tag + '_matrix.txt'), T, fmt='%d', delimiter='\t')
                np.savetxt(os.path.join(OUT, 'holdout', 'heldout_' + tag + '.txt'), np.column_stack((I[pick], J[pick])), fmt='%d')
                # reference structures (10 models)
                for m, (p, x) in enumerate(read_pdb_models(pdb, ch)):
                    bi, bx = bin_structure(p, x, r)
                    np.savetxt(os.path.join(OUT, 'ref', 'ref_{0}_model{1}.txt'.format(tag, m + 1)),
                               np.column_stack((bi, bx)), fmt=['%d', '%.4f', '%.4f', '%.4f'])
                stats.append([ce, ch, r, int(M.sum() / 2), int((M.sum(1) > 0).sum()), len(cand), k])
    with open(os.path.join(OUT, 'sc_matrix_stats.tsv'), 'w') as f:
        f.write('cell\tchrom\tres\tcontacts\tbins_with_contacts\tcontacted_pairs_sep2\theld_out\n')
        for s in stats:
            f.write('\t'.join(map(str, s)) + '\n')
    print('done')


if __name__ == '__main__':
    main()
