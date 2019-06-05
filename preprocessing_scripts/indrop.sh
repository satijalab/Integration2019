#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/pancreas
DLNAME="indrop.tar"

if [ $1 = "orig" ]; then
	URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/6lobg6iwngy8s8s/indrop.tar?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

tar -xvf $RD/indrop.tar -C $RD
gunzip $RD/*human*_umifm_counts.csv.gz
head -1 $RD/GSM2230757_human1_umifm_counts.csv > $RD/indrop.csv
tail -n +2 -q $RD/*human*_umifm_counts.csv >> $RD/indrop.csv
gzip $RD/indrop.csv
rm $RD/*_umifm_counts.csv*
rm $RD/indrop.tar
