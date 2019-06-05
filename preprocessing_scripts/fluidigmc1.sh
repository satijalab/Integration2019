#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/pancreas
DLNAME="fluidigmc1.csv.gz"

if [ $1 = "orig" ]; then
	URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86469/suppl/GSE86469%5FGEO%2Eislet%2Esingle%2Ecell%2Eprocessed%2Edata%2ERSEM%2Eraw%2Eexpected%2Ecounts%2Ecsv%2Egz"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/lww6ngfaq4gw26b/fluidigmc1.csv.gz?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME
