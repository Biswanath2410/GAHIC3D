#run existing reconstruction methods on one contact matrix 
#usage: python run_baselines.py <matrix.txt> <outdir> <method> [alpha] [seed]
#   method : pastis-mds | pastis-pm1 | pastis-pm2 | shrec3d | lordg | minimds | hsa0 | hsa1 | gem | chromsde | chromosome3d
#output : <outdir>/<method>_struct.txt  (bin x y z), <outdir>/<method>_param.txt (seconds)
import sys
import os
import time
import shutil
import subprocess
import numpy as np
from pathlib import Path
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from iced import normalization
import warnings
warnings.filterwarnings('ignore')

TOOLS = Path(__file__).resolve().parent / "tools"
sys.path.insert(0, str(TOOLS / "pastis-0.4.0"))


def remove_nan_col(hic):
    hic = np.nan_to_num(hic)
    col_sum = np.sum(hic, axis=1)
    idx = np.where(col_sum > 0)[0]
    return hic[np.ix_(idx, idx)], idx


def run_pastis(nmat, method, alpha, seed):
    from pastis.optimization import MDS, PM1, PM2
    counts = sparse.coo_matrix(np.triu(nmat, k=1))
    if method == 'pastis-mds':
        X = MDS(alpha=-alpha, beta=1., random_state=seed, max_iter=5000).fit(counts)
    elif method == 'pastis-pm1':
        X = PM1(alpha=-alpha, beta=1., random_state=seed, max_iter=5000).fit(counts)
    elif method == 'pastis-pm2':
        X = PM2(alpha=-alpha, beta=1., random_state=seed, max_iter=5000).fit(counts)
    return X


def run_pastis_native(rawmat, method, alpha, seed):
    """PASTIS run as its own pipeline does it (pastis-pm1 / pastis-pm2 scripts): raw counts,
    4 % lowest-coverage bins filtered, ICE biases passed to the Poisson model"""
    import iced
    from pastis.optimization import PM1, PM2
    counts = iced.filter.filter_low_counts(rawmat.copy(), sparsity=False, percentage=0.04)
    counts = np.nan_to_num(counts)                  # filtered bins come back as NaN
    keep = np.where(counts.sum(axis=0) > 0)[0]
    counts = counts[np.ix_(keep, keep)]
    _, bias = iced.normalization.ICE_normalization(counts.copy(), max_iter=300, output_bias=True)
    bias = np.asarray(bias, dtype=float).flatten()
    c = sparse.coo_matrix(np.triu(counts, k=1))
    model = PM1 if method == 'pastis-pm1-native' else PM2
    Xk = model(alpha=-alpha, beta=1., random_state=seed, max_iter=5000, bias=bias).fit(c)
    X = np.full((len(rawmat), 3), np.nan)
    X[keep] = Xk
    return X


def run_shrec3d(nmat, alpha):
    """ShRec3D (Lesne et al. 2014): wish distances, shortest-path completion, classical MDS"""
    n = len(nmat)
    D = np.zeros_like(nmat)
    nz = nmat > 1e-12 * nmat.max()             # ICE can leave denormal entries on very sparse (single-cell) maps
    D[nz] = nmat[nz] ** (-1. / alpha)
    D = shortest_path(sparse.csr_matrix(D), directed=False)
    D[np.isinf(D)] = D[~np.isinf(D)].max()
    J = np.identity(n) - np.ones((n, n)) / n
    B = -0.5 * J.dot(D ** 2).dot(J)
    val, vec = np.linalg.eigh(B)
    order = np.argsort(val)[::-1][:3]
    return vec[:, order] * np.sqrt(np.maximum(val[order], 0))


def run_lordg(nmat, alpha, workdir):
    workdir = os.path.abspath(workdir)
    """LorDG (Trieu & Cheng 2017), java jar from github.com/BDM-Lab/LorDG"""
    os.makedirs(workdir, exist_ok=True)
    n = len(nmat)
    I, J = np.triu_indices(n, k=1)
    keep = nmat[I, J] > 0
    infile = os.path.join(workdir, 'counts.txt')
    np.savetxt(infile, np.column_stack((I[keep], J[keep], nmat[I, J][keep])), fmt='%d %d %.4f')
    outdir = os.path.join(workdir, 'out')
    os.makedirs(outdir, exist_ok=True)
    config = os.path.join(workdir, 'parameters.txt')
    with open(config, 'w') as fout:
        fout.write("NUM = 1\nOUTPUT_FOLDER = {0}\nINPUT_FILE = {1}\nCONVERT_FACTOR = {2}\nVERBOSE = false\n"
                   "LEARNING_RATE = 1\nMAX_ITERATION = 10000\n".format(outdir, infile, 1. / alpha))
    jar = TOOLS / "LorDG" / "bin" / "3DDistanceBaseLorentz.jar"
    lib = TOOLS / "LorDG" / "lib" / "commons-math3-3.5.jar"
    subprocess.run(['java', '-cp', '{0}:{1}'.format(jar, lib), 'noisy_mds.StructureGeneratorLorentz', config], cwd=workdir,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)   # list form: paths may contain spaces
    pdb = [f for f in os.listdir(outdir) if f.endswith('.pdb')][0]
    mapping = [f for f in os.listdir(outdir) if 'mapping' in f][0]
    cord = []
    with open(os.path.join(outdir, pdb)) as fh:
        for line in fh:
            if line.startswith('ATOM'):
                cord.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    cord = np.array(cord)
    mp = np.loadtxt(os.path.join(outdir, mapping)).astype(int)   # genomic position -> model index
    X = np.full((n, 3), np.nan)
    X[mp[:, 0]] = cord[mp[:, 1]]
    return X


EXT = TOOLS / "ext"


def run_minimds(nmat, idx, alpha, workdir):
    """miniMDS (Rieber & Mahony 2017), full metric MDS, contact-distance exponent set to alpha"""
    os.makedirs(workdir, exist_ok=True)
    res = 1000000                                   # nominal bin size; only bin order matters
    bed = os.path.join(workdir, 'in.bed')
    with open(bed, 'w') as f:
        for a in range(len(idx)):
            for b in range(a, len(idx)):
                if nmat[a, b] > 0:
                    f.write('chrA\t{0}\t{1}\tchrA\t{2}\t{3}\t{4}\n'.format(
                        idx[a] * res, idx[a] * res + res, idx[b] * res, idx[b] * res + res, nmat[a, b]))
    out = os.path.join(workdir, 'out.tsv')
    subprocess.run(['python3', str(TOOLS / 'miniMDS' / 'minimds.py'), '-a', str(alpha), '-o', out, bed],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        fh.readline(); fh.readline()
        start = int(fh.readline())
    s = np.loadtxt(out, skiprows=3)
    X = np.full((len(idx), 3), np.nan)
    pos = {b: k for k, b in enumerate(idx)}
    for row in s:
        b = int(row[0]) + start // res
        if b in pos:
            X[pos[b]] = row[1:4]
    return X


def run_hsa(rawmat, idx, mk, seed, workdir):
    """HSA (Zou et al. 2016), R package hsa (non-official fork of the authors' scripts); raw counts"""
    os.makedirs(workdir, exist_ok=True)
    infile = os.path.join(workdir, 'in.txt')
    np.savetxt(infile, np.column_stack((idx, idx + 1, rawmat)), fmt='%d')
    outp = os.path.join(workdir, 'out')
    subprocess.run('ulimit -s unlimited; Rscript "{0}" "{1}" "{2}" {3} {4}'.format(EXT / 'hsa_run.R', infile, outp, mk, seed),
                   shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    s = np.loadtxt(outp + '.txt')
    X = np.full((len(idx), 3), np.nan)
    pos = {b: k for k, b in enumerate(idx)}
    for row in s:
        if int(row[0]) in pos:
            X[pos[int(row[0])]] = row[2:5]
    return X


def run_gem(nmat, idx, workdir):
    """GEM (Zhu et al. 2018), authors' MATLAB code run in GNU Octave; one conformation (M = 1)"""
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    np.savetxt(os.path.join(workdir, 'HiC.txt'), nmat, delimiter='\t')
    np.savetxt(os.path.join(workdir, 'loci.txt'), idx * 1000000, fmt='%d')
    cmd = "pkg load statistics; addpath('{0}'); GEM('HiC.txt','loci.txt',1E4,1,5E12,0,-1)".format(EXT / 'GEM')
    subprocess.run(['octave', '--no-gui', '-q', '--eval', cmd], cwd=workdir,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return np.loadtxt(os.path.join(workdir, 'conformation1.txt'))


def run_chromsde(nmat, alpha, workdir):
    """ChromSDE (Zhang et al. 2013), authors' MATLAB code in GNU Octave, quadratic SDP, exponent known"""
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    infile, outfile = os.path.join(workdir, 'in.txt'), os.path.join(workdir, 'out.txt')
    np.savetxt(infile, nmat, delimiter='\t')
    cmd = "addpath('{0}'); chromsde_run('{1}','{2}')".format(EXT, infile, outfile)
    subprocess.run(['octave', '--no-gui', '-q', '--eval', cmd], cwd=workdir,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return np.loadtxt(outfile)


def run_chromosome3d(nmat, alpha, workdir, models=5):
    """Chromosome3D (Adhikari et al. 2016), perl + CNS distance-geometry simulated annealing; best of 5 models"""
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    infile = os.path.join(workdir, 'in.txt')
    np.savetxt(infile, nmat, delimiter='\t')
    out = os.path.join(workdir, 'c3d')
    subprocess.run(['perl', str(EXT / 'Chromosome3D' / 'chromosome3D.pl'), '-i', infile, '-o', out,
                    '-a', str(1. / alpha), '-m', str(models)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cord = []
    with open(os.path.join(out, 'in_model1.pdb')) as fh:
        for line in fh:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                cord.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(cord)


def run_nucdynamics(matrix, idx, seed, workdir):
    """NucDynamics (Stevens et al. 2017) via tools/ext/nucdynamics_run.py; the resolution is read from the file name.
    The authors' Cython code needs numpy<2 / Cython<3: set NUCDYN_PYTHON to such an interpreter."""
    import re
    import subprocess
    res = int(re.findall(r'_(\d{5,})(?=[_.])', os.path.basename(matrix))[-1])
    os.makedirs(workdir, exist_ok=True)
    out = os.path.join(workdir, 'nd_struct.txt')
    py = os.environ.get('NUCDYN_PYTHON', '/opt/ndenv/bin/python' if os.path.exists('/opt/ndenv/bin/python') else sys.executable)
    subprocess.run([py, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools', 'ext', 'nucdynamics_run.py'), matrix, str(res), out, str(seed)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    a = np.loadtxt(out)
    X = np.full((len(idx), 3), np.nan)
    pos = {int(b): k for k, b in enumerate(a[:, 0])}
    for k, b in enumerate(idx):
        if int(b) in pos:
            X[k] = a[pos[int(b)], 1:4]
    return X


def main():
    matrix, outdir, method = sys.argv[1:4]
    alpha = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    os.makedirs(outdir, exist_ok=True)

    rawmat = np.loadtxt(matrix)
    rawmat, idx = remove_nan_col(rawmat)
    nmat = normalization.ICE_normalization(rawmat)

    start = time.time()
    if method.endswith('-native'):
        X = run_pastis_native(rawmat, method, alpha, seed)
    elif method.startswith('pastis'):
        X = run_pastis(nmat, method, alpha, seed)
    elif method == 'shrec3d':
        X = run_shrec3d(nmat, alpha)
    elif method == 'lordg':
        wd = os.path.join(outdir, 'lordg_tmp')
        X = run_lordg(nmat, alpha, wd)
        shutil.rmtree(wd, ignore_errors=True)
    elif method == 'nucdynamics':
        wd = os.path.join(outdir, 'nucdynamics_tmp')
        X = run_nucdynamics(os.path.abspath(matrix), idx, seed, wd)
        shutil.rmtree(wd, ignore_errors=True)
    elif method in ('minimds', 'hsa0', 'hsa1', 'gem', 'chromsde', 'chromosome3d'):
        wd = os.path.join(outdir, method + '_tmp')
        if method == 'minimds':
            X = run_minimds(nmat, idx, alpha, wd)
        elif method in ('hsa0', 'hsa1'):
            X = run_hsa(rawmat, idx, int(method[-1]), seed, wd)
        elif method == 'gem':
            X = run_gem(nmat, idx, wd)
        elif method == 'chromsde':
            X = run_chromsde(nmat, alpha, wd)
        else:
            X = run_chromosome3d(nmat, alpha, wd)
        shutil.rmtree(wd, ignore_errors=True)
    elapsed = time.time() - start

    ok = ~np.isnan(X).any(axis=1)
    np.savetxt(os.path.join(outdir, method + '_struct.txt'), np.column_stack((idx[ok], X[ok])),
               fmt=['%d', '%.6f', '%.6f', '%.6f'])
    with open(os.path.join(outdir, method + '_param.txt'), 'w') as f:
        f.write('method\tseconds\n{0}\t{1:.1f}\n'.format(method, elapsed))


if __name__ == '__main__':
    main()
