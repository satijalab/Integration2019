#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/mca
DLNAME1="MCA_500more_dge.tar.gz"
DLNAME2="MCA_CellAssignments.csv"

if [ $1 = "orig" ]; then
	URL1="https://ndownloader.figshare.com/files/10756798"
	URL2="https://ndownloader.figshare.com/files/11083451"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/8hoxr38b2u8eszc/MCA_500more_dge.tar.gz?dl=0"
	URL2="https://www.dropbox.com/s/55ofq23xiiflmlj/MCA_CellAssignments.csv?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL1 --create-dirs -o $RD/$DLNAME1
curl -L $URL2 --create-dirs -o $RD/$DLNAME2

tar -xvzf $RD/$DLNAME1 -C $RD
