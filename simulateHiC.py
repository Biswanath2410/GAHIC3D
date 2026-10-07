#simulate Hi-C contact maps from ground-truth single-cell structures (GSE80280, Stevens et al. 2017)
#usage: python simulateHiC.py <alpha> <coverage> [seed] [nb_dispersion]
#   alpha    : contact-distance exponent, c_ij ~ Poisson(beta * d_ij^alpha)   (main: -2.0)
#   coverage : fraction of the reference sequencing depth                      (main: 1.0)
#   nb_dispersion : if given, negative binomial counts with var = mu + mu^2/r (robustness test)
import sys
import os
import numpy as np
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

cells = [1, 2, 3, 4, 5, 6, 7, 8]
chroms = ['1', '9', '19', 'X']
resolutions = ['1000000', '400000']

# reference depth = mean contacts per bin pair in the original 'simulate4' matrices for cell 1
ref_depth = {'1000000': 130.0, '400000': 20.0}


def get_chrom_length(chrom):    # mm10
    chrom_lengths={
    "1":195471971,
    "9":124595110,
    "19":61431566,
    "X":171031299,
    }
    return chrom_lengths[chrom]


def load_groundtruth(gstfile, chromosome, resolution):
    """Read groundtruth_<cell>_<chr>_<res>.txt -> bin index, xyz.
    Rows starting beyond the chromosome end are dropped (cells 1 and 2, chr19 1 Mb had
    a manually appended 62-63 Mb row in the 2022 backup)."""
    leng = get_chrom_length(chromosome)
    idx = []
    cord = []
    with open(gstfile, 'r') as fh:
        for line in fh:
            l = line.strip().split()
            if len(l) < 6:
                continue
            start = int(float(l[1]))
            if start >= leng:
                continue
            idx.append(start // int(resolution))
            cord.append([float(l[3]), float(l[4]), float(l[5])])
    return np.array(idx, dtype=int), np.array(cord)


def simulate(cord, idx, nbins, alpha, depth, rng, nb_r=None):
    """Poisson contacts from Euclidean distances; full nbins x nbins symmetric matrix."""
    D = squareform(pdist(cord))
    n = len(idx)
    iu = np.triu_indices(n, k=1)
    lam = D[iu] ** alpha
    beta = depth * len(lam) / lam.sum()
    mu = beta * lam
    if nb_r is None:
        cnt = rng.poisson(mu)
    else:
        cnt = rng.poisson(rng.gamma(nb_r, mu / nb_r))      # gamma-Poisson = negative binomial
    mat = np.zeros((nbins, nbins))
    mat[idx[iu[0]], idx[iu[1]]] = cnt
    mat = mat + mat.T
    return mat


def main():

    alpha = float(sys.argv[1]) if len(sys.argv) > 1 else -2.0
    coverage = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2022
    nb_r = float(sys.argv[4]) if len(sys.argv) > 4 else None

    curr_dir = Path(__file__).resolve().parent
    gst_dir = curr_dir / "data" / "groundtruth"
    out_dir = curr_dir / "data" / "simulated" / ("alpha{0}_cov{1}".format(abs(alpha), coverage) + ("" if nb_r is None else "_nb{0:g}".format(nb_r)))
    os.makedirs(out_dir, exist_ok=True)

    rng = np.random.default_rng(seed)
    for ce in cells:
        for ch in chroms:
            for r in resolutions:
                gstfile = gst_dir / ("groundtruth_" + str(ce) + "_" + ch + "_" + r + ".txt")
                idx, cord = load_groundtruth(gstfile, ch, r)
                nbins = int(get_chrom_length(ch) / int(r)) + 1
                mat = simulate(cord, idx, nbins, alpha, ref_depth[r] * coverage, rng, nb_r)
                fname = "simHiC_Cell_" + str(ce) + "_" + ch + "_" + r + "_matrix.txt"
                np.savetxt(out_dir / fname, mat, fmt='%d', delimiter='\t')
    print('saved in', out_dir)


if __name__ == '__main__':
    main()
