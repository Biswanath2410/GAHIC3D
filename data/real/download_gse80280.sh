#!/bin/bash
# Download GSE80280 (Stevens et al. 2017) into data/real/raw/ (~1.4 GB)
set -e
RAW="$(cd "$(dirname "$0")" && pwd)/raw"
mkdir -p "$RAW"
get() { [ -s "$RAW/$2" ] || curl -fL --retry 3 -o "$RAW/$2" "https://ftp.ncbi.nlm.nih.gov/geo/samples/$1/$2"; }
get GSM2123nnn/GSM2123564/suppl GSM2123564_Haploid_mESC_population_hic.txt.gz
n=1
for gsm in 2219497 2219498 2219499 2219500 2219501 2219502 2219503 2219504; do
  get GSM2219nnn/GSM$gsm/suppl GSM${gsm}_Cell_${n}_contact_pairs.txt.gz
  get GSM2219nnn/GSM$gsm/suppl GSM${gsm}_Cell_${n}_genome_structure_model.pdb.gz
  n=$((n + 1))
done
ls -la "$RAW"
