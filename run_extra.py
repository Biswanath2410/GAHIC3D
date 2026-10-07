"""Additional experiments :
   (1) depth series (100, 30, 10, 3, 1, 0.3% depth): 1 Mb, all chromosomes; 400 kb, chr19
   (2) parameter sensitivity of GAHIC3D (chr19, 400 kb, full and 1% depth)
   (3) convergence: GAHIC3D vs multi-start with 120 relaxations (chr19, 400 kb)
   (4) GAHIC3D on raw counts with ICE biases (--raw-bias; Supplementary Note 1)
usage: python run_extra.py [ncores]"""
import os
import sys
import subprocess
from multiprocessing import Pool
import run_benchmark as rb

GA = ['--variant', 'memetic', '--fitness', 'poisson']
DEPTH_METHODS = ['pastis-pm1', 'pastis-pm2', 'pastis-mds', 'lordg', 'shrec3d', 'minimds', 'chromsde']
SENS = {'GAHIC3D-pop5': ['--memetic-pop', '5'], 'GAHIC3D-pop20': ['--memetic-pop', '20'],
        'GAHIC3D-relax20': ['--n-relax', '20'], 'GAHIC3D-relax80': ['--n-relax', '80'],
        'GAHIC3D-pivot0': ['--pivot-rate', '0'], 'GAHIC3D-pivot0.6': ['--pivot-rate', '0.6'],
        'GAHIC3D-iter100': ['--local-iter', '100'], 'GAHIC3D-iter400': ['--local-iter', '400'],
        'GAHIC3D-a1.5': ['--alpha', '1.5'], 'GAHIC3D-a2.5': ['--alpha', '2.5'], 'GAHIC3D-a3': ['--alpha', '3']}
CONV = {'GAHIC3D-relax120': GA + ['--n-relax', '120'],
        'MultiStart-relax120': ['--variant', 'multistart', '--fitness', 'poisson', '--n-relax', '120']}


def job_list():
    jobs = []
    for cov in ['0.3', '0.1', '0.03', '0.003']:
        for ch in ['19', '9', 'X', '1']:
            for ce in rb.cells:
                jobs.append(('GAHIC3D', ce, ch, '1000000', cov, GA))
                jobs += [(m, ce, ch, '1000000', cov, None) for m in DEPTH_METHODS]
    for cov in ['0.3', '0.1', '0.03', '0.003']:
        for ce in rb.cells:
            jobs.append(('GAHIC3D', ce, '19', '400000', cov, GA))
            jobs += [(m, ce, '19', '400000', cov, None) for m in DEPTH_METHODS]
    for cov in ['1.0', '0.01']:
        for ce in rb.cells:
            jobs += [(m, ce, '19', '400000', cov, GA + a) for m, a in SENS.items()]
            jobs += [(m, ce, '19', '400000', cov, a) for m, a in CONV.items()]
    rawbias = GA + ['--raw-bias']
    for cov in ['1.0', '0.01']:
        jobs += [('GAHIC3D-rawbias', ce, ch, '1000000', cov, rawbias) for ch in ['19', '9', 'X', '1'] for ce in rb.cells]
        jobs += [('GAHIC3D-rawbias', ce, '19', '400000', cov, rawbias) for ce in rb.cells]
    return jobs


_run = subprocess.run
def run_with_timeout(*a, **k):
    k.setdefault('timeout', 1800)
    try:
        return _run(*a, **k)
    except subprocess.TimeoutExpired:
        return subprocess.CompletedProcess(a[0], 1, '', 'timeout')
rb.subprocess.run = run_with_timeout


if __name__ == '__main__':
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    os.makedirs(rb.EVAL, exist_ok=True)
    jobs = job_list()
    print(len(jobs), 'jobs', flush=True)
    with Pool(ncores) as pool:
        for k, (name, status) in enumerate(pool.imap_unordered(rb.run, jobs, chunksize=1)):
            print(k, name, status, flush=True)
