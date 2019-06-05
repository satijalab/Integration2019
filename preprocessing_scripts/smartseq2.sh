#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/pancreas
DLNAME="smartseq2.zip"

if [ $1 = "orig" ]; then
	URL="https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/ksl63ljwil08wmx/smartseq2.zip?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

unzip $RD/smartseq2.zip -d $RD
mv $RD/pancreas_refseq_rpkms_counts_3514sc.txt $RD/smartseq2.txt
rm $RD/smartseq2.zip
