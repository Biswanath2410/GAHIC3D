#evaluation of a predicted structure against the ground-truth structure
# RMSD after Kabsch superposition with scaling,
#RMSE between the two Euclidean distance matrices, Pearson correlation of distances.
#Bins are matched by genomic index.
#usage: python RMSD.py <est_struct.txt> <groundtruth.txt> <chromosome> <resolution> <METHOD> <CELL> <out.csv> [extra tag]
import sys
import os
import csv
import math
import numpy as np
from scipy import stats
from sklearn.metrics.pairwise import euclidean_distances
from simulateHiC import load_groundtruth
from GAHIC3D import getTransformation


def rmsd(X, Y):
    n, _ = X.shape
    return (((X - Y) ** 2).sum() / n) ** 0.5


def pearson(mat1, mat2):
    """Pearson correlation of distance matrices, ignoring zeroes (as 2022)"""
    vec1 = mat1.flatten()
    vec2 = mat2.flatten()
    nonzero = (vec1 != 0) & (vec2 != 0)
    r, p = stats.pearsonr(vec1[nonzero], vec2[nonzero])
    return r


def spearman(mat1, mat2):
    iu = np.triu_indices(len(mat1), k=1)
    return stats.spearmanr(mat1[iu], mat2[iu])[0]


def load_est(estfile):
    est = np.loadtxt(estfile)
    return est[:, 0].astype(int), est[:, 1:4]


def compare(est_idx, est, gst_idx, gst):
    common, ie, ig = np.intersect1d(est_idx, gst_idx, return_indices=True)
    X, Y = est[ie], gst[ig]
    # Hi-C cannot tell a structure from its mirror image, so the mirror is allowed
    # (2022 evaluation used reflection=False, which makes RMSD depend on a random chirality)
    Xa, Ya = getTransformation(X, Y, reflection=True)   # scale X to Y, rotate/reflect
    Xn, Yn = getTransformation(X, Y, reflection=False)
    dx = euclidean_distances(Xa)
    dy = euclidean_distances(Ya)
    MSE = np.square(dy - dx).mean()
    return dict(RMSD=rmsd(Xa, Ya), RMSE=math.sqrt(MSE), PCoeff=pearson(dx, dy),
                SCoeff=spearman(dx, dy), NBINS=len(common), RMSD_noreflect=rmsd(Xn, Yn))


def main():
    estfile, gstfile, chromosome, resolution, method, cell, outcsv = sys.argv[1:8]
    tag = sys.argv[8] if len(sys.argv) > 8 else ''
    est_idx, est = load_est(estfile)
    gst_idx, gst = load_groundtruth(gstfile, chromosome, resolution)
    res = compare(est_idx, est, gst_idx, gst)
    header = ['METHOD', 'CELL', 'CHROMOSOME', 'RESOLUTION', 'TAG', 'RMSD', 'RMSE', 'PCoeff', 'SCoeff', 'NBINS', 'RMSD_noreflect']
    data = [method, cell, chromosome, resolution, tag, res['RMSD'], res['RMSE'], res['PCoeff'], res['SCoeff'], res['NBINS'], res['RMSD_noreflect']]
    with open(outcsv, 'a') as f:
        writer = csv.writer(f)
        if f.tell() == 0:
            writer.writerow(header)
        writer.writerow(data)


if __name__ == '__main__':
    main()
