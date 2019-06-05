#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/pancreas
DLNAME="celseq.csv.gz"

if [ $1 = "orig" ]; then
	URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81076/suppl/GSE81076%5FD2%5F3%5F7%5F10%5F17%2Etxt%2Egz"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/v6j4h4y6gtx9i0r/celseq.csv.gz?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME
