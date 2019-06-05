#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/ica
DLNAME="ica_bone_marrow.h5"
DLNAME2="ica_scrublet_results.tar.gz"

if [ $1 = "orig" ]; then
	URL="https://s3.amazonaws.com/preview-ica-expression-data/ica_bone_marrow_h5.h5"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/xe5tithw1xjxrfs/ica_bone_marrow.h5?dl=0"
else
	echo "bad option"; exit;
fi

URL2="https://www.dropbox.com/s/q0715fq5j9m2x1c/ICA_scrublet_results.tar.gz?dl=0"

curl -L $URL --create-dirs -o $RD/$DLNAME
curl -L $URL2 --create-dirs -o $RD/$DLNAME2

tar -xvf $RD/$DLNAME2 -C $RD

