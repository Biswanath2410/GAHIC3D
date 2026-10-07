#benchmark driver 
#usage: python run_benchmark.py [ncores]
#every job: reconstruct -> evaluate with RMSD.py -> results/evals/<job>.csv ; finished jobs are skipped
import sys
import os
import glob
import subprocess
import itertools
from multiprocessing import Pool

cells = [1, 2, 3, 4, 5, 6, 7, 8]
chroms = ['19', '9', 'X', '1']
resolutions = ['1000000', '400000']
coverages = ['1.0', '0.01']
baselines = ['shrec3d', 'lordg', 'pastis-mds', 'pastis-pm1', 'pastis-pm2',
             'minimds', 'hsa1', 'chromsde', 'chromosome3d', 'gem']

ROOT = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(ROOT, 'data', 'simulated')
GST = os.path.join(ROOT, 'data', 'groundtruth')
OUT = os.path.join(ROOT, 'results', 'structures')
EVAL = os.path.join(ROOT, 'results', 'evals')


def job_list():
    jobs = []
    for r, cov, ch, ce in itertools.product(resolutions, coverages, chroms, cells):
        # main benchmark: GAHIC3D (memetic, Poisson) + existing methods
        jobs.append(('GAHIC3D', ce, ch, r, cov, ['--variant', 'memetic', '--fitness', 'poisson']))
        for m in baselines:
            # PASTIS-PM2 (alpha estimated) is the slowest method (6 min chr9, 19 min chrX per 400 kb
            # structure); at 400 kb it is run on chr19 (both depths) and chr9 (full depth) only.
            # PASTIS-PM1 (same Poisson model, alpha given) is run everywhere.
            if m == 'pastis-pm2' and r == '400000' and (ch in ('X', '1') or (ch == '9' and cov == '0.01')):
                continue
            # Chromosome3D (CNS annealing, up to 15 min) and HSA (up to 10 min) per 400 kb structure -> 1 Mb only
            if m in ('chromosome3d', 'hsa1') and r == '400000':
                continue
            # GEM in Octave: 6 min for chr19 at 1 Mb, >30 min for chr9 at 1 Mb -> chr19, 1 Mb only
            if m == 'gem' and not (ch == '19' and r == '1000000'):
                continue
            jobs.append((m, ce, ch, r, cov, None))
        # optimiser control at 400 kb: same objective, same number of relaxations, no evolution
        if r == '400000':
            jobs.append(('MultiStart', ce, ch, r, cov, ['--variant', 'multistart', '--fitness', 'poisson']))
        # objective / operator ablation on chr19 400 kb
        if r == '400000' and ch == '19':
            jobs.append(('GA-poissonNZ', ce, ch, r, cov, ['--variant', 'memetic', '--fitness', 'poissonnz']))
            jobs.append(('GA-wstress', ce, ch, r, cov, ['--variant', 'memetic', '--fitness', 'wstress']))
            jobs.append(('GA-stress', ce, ch, r, cov, ['--variant', 'memetic', '--fitness', 'stress']))
            jobs.append(('GA-generational', ce, ch, r, cov, ['--variant', 'ga', '--fitness', 'poisson']))
    # robustness: data simulated with alpha=-3 and negative-binomial noise (r=20), methods keep alpha=2
    for cov, ch, ce in itertools.product(['1.0_nb20', '0.01_nb20'], chroms, cells):
        rs = '1000000'
        jobs.append(('GAHIC3D', ce, ch, rs, cov, ['--variant', 'memetic', '--fitness', 'poisson']))
        for m in baselines[:5]:                 # robustness set: the original five baselines
            jobs.append((m, ce, ch, rs, cov, None))
    return jobs


def run(job):
    method, ce, ch, r, cov, args = job
    tag = 'cov' + cov
    name = '{0}_{1}_{2}_{3}_{4}'.format(method, ce, ch, r, tag)
    evalfile = os.path.join(EVAL, name + '.csv')
    if os.path.exists(evalfile):
        return name, 'skip'
    simdir = ('alpha3.0_cov' if 'nb' in cov else 'alpha2.0_cov') + cov
    matrix = os.path.join(SIM, simdir, 'simHiC_Cell_{0}_{1}_{2}_matrix.txt'.format(ce, ch, r))
    outdir = os.path.join(OUT, tag, '{0}_{1}_{2}'.format(ce, ch, r))
    os.makedirs(outdir, exist_ok=True)
    if args is not None:
        prefix = os.path.join(outdir, method)
        cmd = ['python3', 'GAHIC3D.py', matrix, prefix, '--seed', str(ce)] + args
        struct = prefix + '_struct.txt'
    else:
        cmd = ['python3', 'run_baselines.py', matrix, outdir, method, '2.0', str(ce)]
        struct = os.path.join(outdir, method + '_struct.txt')
    p = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if not os.path.exists(struct):
        with open(os.path.join(EVAL, name + '.failed'), 'w') as f:
            f.write(p.stderr[-3000:])
        return name, 'failed'
    gst = os.path.join(GST, 'groundtruth_{0}_{1}_{2}.txt'.format(ce, ch, r))
    subprocess.run(['python3', 'RMSD.py', struct, gst, ch, r, method, str(ce), evalfile, tag],
                   cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return name, 'ok'


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    os.makedirs(EVAL, exist_ok=True)
    jobs = job_list()
    # small problems first so that partial results are usable
    size = {'19': 1, '9': 2, 'X': 3, '1': 4}
    jobs.sort(key=lambda j: (j[3] == '400000', size[j[2]], j[0] == 'GA-2022'))
    print(len(jobs), 'jobs', flush=True)
    with Pool(ncores) as pool:
        for k, (name, status) in enumerate(pool.imap_unordered(run, jobs)):
            print(k, name, status, flush=True)
    # merge
    files = sorted(glob.glob(os.path.join(EVAL, '*.csv')))
    with open(os.path.join(ROOT, 'results', 'benchmark_all.csv'), 'w') as fout:
        for i, f in enumerate(files):
            lines = open(f).read().strip().split('\n')
            fout.write('\n'.join(lines if i == 0 else lines[1:]) + '\n')


if __name__ == '__main__':
    main()
