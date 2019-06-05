#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/mca
DLNAME="mca_tm_facs.zip"
DLNAME2="annotations_facs.csv"

if [ $1 = "orig" ]; then
	URL="https://ndownloader.figshare.com/articles/5715040/versions/1"
	URL2="https://github.com/czbiohub/tabula-muris/raw/e5f4ccc6bb613fdbacb02d92ae8751ead8b21f47/00_data_ingest/18_global_annotation_csv/annotations_facs.csv"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/u136m7hrltrulln/mca_tm_facs.zip?dl=0"
	URL2="https://www.dropbox.com/s/39d0tvkqe3o0n5h/annotations_facs.csv?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

unzip $RD/$DLNAME -d $RD
unzip $RD/FACS.zip -d $RD

curl -L $URL2 -o $RD/$DLNAME2