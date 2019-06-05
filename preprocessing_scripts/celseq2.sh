#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/pancreas
DLNAME="celseq2.csv.gz"

if [ $1 = "orig" ]; then
	URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/ga1qn8njcb7qxcx/celseq2.csv.gz?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME
