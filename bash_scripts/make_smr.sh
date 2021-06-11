
#! /bin/bash

PATH_LD=$1
PATH_SUMSTATS=$2
PATH_EQTL=$3
PATH_OUT=$4

mkdir -p $PATH_OUT
set -e

gwas_base_name=$(basename $PATH_SUMSTATS)
gwas_name=${gwas_base_name%%\.tsv}

LD_name=$(basename $PATH_LD)

for file in $PATH_EQTL/*.summary
do
	base_name=$(basename $file)
	tissue=${base_name%%\.summary}
	./smr_Linux --bfile $PATH_LD --gwas-summary $PATH_SUMSTATS --beqtl-summary "$PATH_EQTL/$tissue" --diff-freq-prop 0.1 --out "$PATH_OUT/${gwas_name}_${LD_name}_${tissue}"

done
