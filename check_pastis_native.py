# PASTIS PM1/PM2 run on raw counts with their own bias model (as the pastis-pm1/pm2 scripts do)
#instead of ICE normalised input. Results: results/evals/pastis-pm*-native_*.csv
import itertools
from multiprocessing import Pool
import run_benchmark as rb

jobs = []
for cov, ch, ce in itertools.product(['0.01', '1.0'], rb.chroms, rb.cells):
    jobs.append(('pastis-pm1-native', ce, ch, '1000000', cov, None))
    jobs.append(('pastis-pm2-native', ce, ch, '1000000', cov, None))
for ch, ce in itertools.product(rb.chroms, rb.cells):
    jobs.append(('pastis-pm1-native', ce, ch, '400000', '0.01', None))

if __name__ == '__main__':
    with Pool(2) as p:
        for k, r in enumerate(p.imap_unordered(rb.run, jobs)):
            print(k, r, flush=True)
