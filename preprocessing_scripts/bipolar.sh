#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/bipolar
DLNAME="bipolar.txt.gz"

if [ $1 = "orig" ]; then
	URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81904/suppl/GSE81904_BipolarUMICounts_Cell2016.txt.gz"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/wwrgrnluvpuk35r/bipolar.txt.gz?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

