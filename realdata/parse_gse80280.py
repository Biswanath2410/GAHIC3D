#parsers for GSE80280 (Stevens et al. 2017) files
import gzip
import numpy as np

CHAIN = {'1': 'a', '9': 'i', '19': 's', 'X': 't'}      # PDB chain letter: chr1=a, chr2=b, ..., chr19=s, chrX=t


def read_pdb_models(pdbgz, chromosome):
    """Returns list of (positions, xyz) per model for one chromosome (100 kb particles)."""
    models, pos, xyz = [], [], []
    with gzip.open(pdbgz, 'rt') as fh:
        for line in fh:
            if line.startswith('MODEL'):
                pos, xyz = [], []
            elif line.startswith('ENDMDL'):
                models.append((np.array(pos), np.array(xyz)))
            elif line.startswith('HETATM') and line[21] == CHAIN[chromosome] and line[17:20] == 'chr':
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                pos.append(int(line.split()[-1]))
    return models


def bin_structure(pos, xyz, resolution):
    """Average particle coordinates within each genomic bin -> bin index, xyz"""
    b = pos // resolution
    ub = np.unique(b)
    return ub, np.array([xyz[b == u].mean(axis=0) for u in ub])


def read_sc_contacts(pairsgz, chromosome, resolution, nbins):
    """Single-cell contact pairs -> intra-chromosomal binned count matrix"""
    c = 'chr' + chromosome
    M = np.zeros((nbins, nbins))
    with gzip.open(pairsgz, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            a, pa, b, pb = line.split()[:4]
            if a == c and b == c:
                i, j = int(pa) // resolution, int(pb) // resolution
                if i != j and i < nbins and j < nbins:
                    M[i, j] += 1
                    M[j, i] += 1
    return M


def read_population(binfile, chromosome, nbins):
    """population_hic_<res>.txt (chr, bin_i, bin_j, count; aggregated from the 50 kb GEO file)"""
    c = 'chr' + chromosome
    M = np.zeros((nbins, nbins))
    with open(binfile) as fh:
        for line in fh:
            ch, i, j, v = line.split()
            if ch == c:
                i, j = int(i), int(j)
                if i != j and i < nbins and j < nbins:
                    M[i, j] += float(v)
                    M[j, i] += float(v)
    return M
