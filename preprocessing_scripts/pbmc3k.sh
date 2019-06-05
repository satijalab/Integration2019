#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/pbmc/pbmc3k/
DLNAME="pbmc3k.tar"

if [ $1 = "orig" ]; then
    URL="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
elif [ $1 = "dropbox" ]; then
    URL="https://www.dropbox.com/s/shkmmbjf0ia191j/pbmc3k_filtered_gene_bc_matrices.tar.gz?dl=0"
else
	echo "bad option"; exit;
fi

curl -L $URL --create-dirs -o $RD/$DLNAME
tar -xvf $RD/$DLNAME -C $RD
