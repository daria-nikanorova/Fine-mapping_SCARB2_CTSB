#! /bin/bash

PATH_GWAS=$1
PATH_EQTL=$2
PATH_OUT=$3

set -e

for file in $PATH_EQTL/*.txt
do
	base_name=$(basename $file)
        eqtl_name=${base_name%%\.txt}
	output_name=$"$PATH_OUT/${eqtl_name}"
        if [ ! -e $output_name ]; then
		echo $output_name
                echo "Start $eqtl_name!"
                Rscript make_coloc.R $PATH_GWAS $file $output_name
	fi
done


