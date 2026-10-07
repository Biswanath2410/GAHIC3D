#!/bin/bash
# Reproduce all results, tables and figures of the paper.
# Usage: bash run_all.sh [ncores]
set -e
cd "$(dirname "$0")"
source env.sh
N=${1:-2}

# data
[ -s data/real/raw/GSM2123564_Haploid_mESC_population_hic.txt.gz ] || bash data/real/download_gse80280.sh

python3 tests/test_core.py

# simulated Hi-C (the 10% depth set is provided in data/simulated/alpha2.0_cov0.1)
python3 simulateHiC.py -2.0 1.0 2022
python3 simulateHiC.py -2.0 0.01 2024
python3 simulateHiC.py -3.0 1.0 3031 20
python3 simulateHiC.py -3.0 0.01 3032 20
python3 simulateHiC.py -2.0 0.3 2022
python3 simulateHiC.py -2.0 0.03 2022
python3 simulateHiC.py -2.0 0.003 2022

# benchmark on simulated data
python3 run_benchmark.py $N
python3 check_pastis_native.py
python3 run_extra.py $N

# experimental data (GSE80280)
bash realdata/bin_population.sh data/real/raw/GSM2123564_Haploid_mESC_population_hic.txt.gz
python3 realdata/prepare_real.py
python3 realdata/run_real.py $N

# figures and tables
python3 plot_results.py
python3 plot_extra.py
python3 realdata/plot_real.py
python3 realdata/plot_real_example.py 7 9 400000
python3 verify_objectives.py
python3 manuscript_numbers.py > results/manuscript_numbers.txt
