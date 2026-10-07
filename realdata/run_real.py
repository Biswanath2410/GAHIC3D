#real-data validation on GSE80280 (run realdata/prepare_real.py first)
#  A  single-cell contacts of cell c  -> compare with the published structure of cell c (model 1)
#  B  single-cell contacts with 10 % of contacted pairs held out -> are held-out pairs close in 3D? (AUC)
#  C  population Hi-C -> compare with the mean distance matrix of the 8 published single-cell structures
#usage: python realdata/run_real.py [ncores]      (finished jobs are skipped)
import os
import sys
import csv
import glob
import itertools
import subprocess
import numpy as np
from multiprocessing import Pool
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from RMSD import compare
from sklearn.metrics.pairwise import euclidean_distances

DATA = os.path.join(ROOT, 'data', 'real')
OUT = os.path.join(ROOT, 'results', 'real')
cells = [1, 2, 3, 4, 5, 6, 7, 8]
chroms = ['19', '9', 'X', '1']
resolutions = ['1000000', '400000']
fast = ['GAHIC3D', 'shrec3d', 'lordg', 'pastis-mds', 'pastis-pm1', 'minimds', 'hsa1', 'chromsde']


def methods_for(ch, r, analysis):
    m = list(fast)
    if r == '400000':
        m.remove('hsa1')                        # HSA: up to 10 min per 400 kb structure -> 1 Mb only
        if ch != '19':
            m.remove('chromsde')                # ChromSDE: >10 min or no convergence on sparse 400 kb input beyond chr19
    if r == '1000000':
        m += ['chromosome3d']
        if analysis != 'B':
            m += ['pastis-pm2']
        if ch == '19' and analysis != 'B':
            m += ['gem']
    return m


def job_list():
    jobs = []
    for r, ch in itertools.product(resolutions, chroms):
        for m in methods_for(ch, r, 'C'):
            jobs.append(('C', m, 0, ch, r))
        for ce in cells:
            for m in methods_for(ch, r, 'A'):
                jobs.append(('A', m, ce, ch, r))
            if r == '400000' and ch == '1':
                continue                        # hold-out test at 400 kb: chr19, 9, X (runtime)
            for m in methods_for(ch, r, 'B'):
                jobs.append(('B', m, ce, ch, r))
    return jobs


def call(cmd, timeout=900):
    """run a reconstruction; kill the whole process group after `timeout` s (counted as a failure)"""
    import signal
    p = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    try:
        p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()


def reconstruct(method, matrix, outdir, seed):
    os.makedirs(outdir, exist_ok=True)
    if method == 'GAHIC3D':
        call(['python3', 'GAHIC3D.py', matrix, os.path.join(outdir, method), '--seed', str(seed)])
    else:
        call(['python3', 'run_baselines.py', matrix, outdir, method, '2.0', str(seed)])
    f = os.path.join(outdir, method + '_struct.txt')
    if not os.path.exists(f):
        return None
    e = np.loadtxt(f)
    return e[:, 0].astype(int), e[:, 1:4]


def load_ref(ce, ch, r, model=1):
    a = np.loadtxt(os.path.join(DATA, 'ref', 'ref_Cell_{0}_{1}_{2}_model{3}.txt'.format(ce, ch, r, model)))
    return a[:, 0].astype(int), a[:, 1:4]


def auc_heldout(idx, X, M_full, heldout):
    """stratified AUC: each held-out pair vs all non-contact pairs at the same genomic separation"""
    pos = {b: k for k, b in enumerate(idx)}
    D = euclidean_distances(X)
    scores = []
    for i, j in heldout:
        if i not in pos or j not in pos:
            continue
        s = j - i
        d_pos = D[pos[i], pos[j]]
        neg = []
        for a in range(M_full.shape[0] - s):
            b = a + s
            if M_full[a, b] == 0 and a in pos and b in pos:
                neg.append(D[pos[a], pos[b]])
        if len(neg) == 0:
            continue
        neg = np.array(neg)
        scores.append(((neg > d_pos).sum() + 0.5 * (neg == d_pos).sum()) / len(neg))
    return (np.mean(scores) if scores else np.nan), len(scores)


def ensemble_mean_distance(ch, r):
    refs = [load_ref(ce, ch, r) for ce in cells]
    allb = np.unique(np.concatenate([b for b, _ in refs]))
    S = np.full((len(cells), len(allb), len(allb)), np.nan)
    for k, (b, x) in enumerate(refs):
        D = euclidean_distances(x)
        ii = np.searchsorted(allb, b)
        S[k][np.ix_(ii, ii)] = D
    return allb, np.nanmean(S, axis=0), refs


def run(job):
    try:
        return run_job(job)
    except Exception as e:
        name = '{0}_{1}_{2}_{3}_{4}'.format(*job)
        with open(os.path.join(OUT, 'evals', name + '.failed'), 'w') as f:
            f.write(repr(e))
        return name, 'failed: ' + repr(e)[:80]


def run_job(job):
    analysis, method, ce, ch, r = job
    name = '{0}_{1}_{2}_{3}_{4}'.format(analysis, method, ce, ch, r)
    evalfile = os.path.join(OUT, 'evals', name + '.csv')
    if os.path.exists(evalfile):
        return name, 'skip'
    tag = 'Cell_{0}_{1}_{2}'.format(ce, ch, r)
    if analysis == 'A':
        matrix = os.path.join(DATA, 'sc', 'sc_' + tag + '_matrix.txt')
    elif analysis == 'B':
        matrix = os.path.join(DATA, 'holdout', 'train_' + tag + '_matrix.txt')
    else:
        matrix = os.path.join(DATA, 'pop', 'pop_{0}_{1}_matrix.txt'.format(ch, r))
    outdir = os.path.join(OUT, 'structures', analysis, '{0}_{1}_{2}'.format(ce, ch, r))
    res = reconstruct(method, matrix, outdir, max(ce, 1))
    if res is None:
        open(os.path.join(OUT, 'evals', name + '.failed'), 'w').close()
        return name, 'failed'
    idx, X = res
    row = dict(ANALYSIS=analysis, METHOD=method, CELL=ce, CHROMOSOME=ch, RESOLUTION=r)
    if analysis == 'A':
        gi, g = load_ref(ce, ch, r)
        c = compare(idx, X, gi, g)
        row.update(RMSD=c['RMSD'], PCoeff=c['PCoeff'], SCoeff=c['SCoeff'], NBINS=c['NBINS'])
    elif analysis == 'B':
        M_full = np.loadtxt(os.path.join(DATA, 'sc', 'sc_' + tag + '_matrix.txt'))
        held = np.loadtxt(os.path.join(DATA, 'holdout', 'heldout_' + tag + '.txt'), dtype=int, ndmin=2)
        auc, n = auc_heldout(idx, X, M_full, held)
        row.update(AUC=auc, NHELD=n)
    else:
        allb, Dmean, refs = ensemble_mean_distance(ch, r)
        common, ie, ig = np.intersect1d(idx, allb, return_indices=True)
        De = euclidean_distances(X[ie])
        Dm = Dmean[np.ix_(ig, ig)]
        iu = np.triu_indices(len(common), k=1)
        ok = ~np.isnan(Dm[iu])
        row.update(SCoeff_mean=stats.spearmanr(De[iu][ok], Dm[iu][ok])[0],
                   PCoeff_mean=stats.pearsonr(De[iu][ok], Dm[iu][ok])[0], NBINS=len(common))
        per = []
        for b, x in refs:
            c = compare(idx, X, b, x)
            per.append(c['SCoeff'])
        row.update(SCoeff_cells=np.mean(per))
    with open(evalfile, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)
    return name, 'ok'


def ceilings():
    """agreement among the 10 published models of each cell (upper bound for A)"""
    rows = []
    for r, ch, ce in itertools.product(resolutions, chroms, cells):
        gi, g = load_ref(ce, ch, r, 1)
        v = [compare(*load_ref(ce, ch, r, m), gi, g) for m in range(2, 11)]
        rows.append([ce, ch, r, np.mean([x['RMSD'] for x in v]), np.mean([x['PCoeff'] for x in v])])
    with open(os.path.join(OUT, 'model_agreement.csv'), 'w') as f:
        f.write('CELL,CHROMOSOME,RESOLUTION,RMSD_models,PCoeff_models\n')
        for x in rows:
            f.write(','.join(map(str, x)) + '\n')
    # between-cell agreement (how different are two cells?)
    rows = []
    for r, ch in itertools.product(resolutions, chroms):
        for a, b in itertools.combinations(cells, 2):
            c = compare(*load_ref(a, ch, r), *load_ref(b, ch, r))
            rows.append([a, b, ch, r, c['RMSD'], c['PCoeff']])
    with open(os.path.join(OUT, 'between_cell_agreement.csv'), 'w') as f:
        f.write('CELL_A,CELL_B,CHROMOSOME,RESOLUTION,RMSD,PCoeff\n')
        for x in rows:
            f.write(','.join(map(str, x)) + '\n')


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    os.makedirs(os.path.join(OUT, 'evals'), exist_ok=True)
    ceilings()
    jobs = job_list()
    size = {'19': 1, '9': 2, 'X': 3, '1': 4}
    jobs.sort(key=lambda j: (j[4] == '400000', size[j[3]], j[1] in ('gem', 'chromosome3d', 'pastis-pm2')))
    print(len(jobs), 'jobs', flush=True)
    with Pool(ncores) as pool:
        for k, (name, status) in enumerate(pool.imap_unordered(run, jobs)):
            print(k, name, status, flush=True)


if __name__ == '__main__':
    main()
