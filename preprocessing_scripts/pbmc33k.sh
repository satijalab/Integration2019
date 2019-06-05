#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/pbmc/pbmc33k/
DLNAME="pbmc33k.tar"

if [ $1 = "orig" ]; then
    URL="http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc33k/pbmc33k_filtered_gene_bc_matrices.tar.gz"
elif [ $1 = "dropbox" ]; then
    URL="https://www.dropbox.com/s/yviil2pk1kvmrq6/pbmc33k_filtered_gene_bc_matrices.tar?dl=0"
else
	echo "bad option"; exit;
fi

curl -L $URL --create-dirs -o $RD/$DLNAME
tar -xvf $RD/$DLNAME -C $RD
