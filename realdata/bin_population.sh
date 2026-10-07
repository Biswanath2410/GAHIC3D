#!/bin/bash
# Bin the GSE80280 population Hi-C (GSM2123564, 50 kb contact list) to 1 Mb and 400 kb,
# intra-chromosomal contacts of chr1, 9, 19 and X only.
# output: data/real/population_hic_<res>.txt   (chrom  binA  binB  count;  binA <= binB, bin = start // res)
# usage:  bash realdata/bin_population.sh <path to GSM2123564_Haploid_mESC_population_hic.txt[.gz]>
set -e
IN=${1:-../GSE80280/GSM2123564_Haploid_mESC_population_hic.txt}
cd "$(dirname "$0")/.."
CAT=cat; case "$IN" in *.gz) CAT="gzip -dc";; esac
for RES in 1000000 400000; do
  $CAT "$IN" | awk -v r=$RES 'BEGIN{OFS="\t"}
    $1==$4 && ($1=="chr1"||$1=="chr9"||$1=="chr19"||$1=="chrX") {
      a=int($2/r); b=int($5/r); if (a>b) {t=a; a=b; b=t}
      n[$1 OFS a OFS b]+=$7 }
    END{for (k in n) print k, n[k]}' | sort -k1,1 -k2,2n -k3,3n > data/real/population_hic_${RES}.txt
done
echo done
