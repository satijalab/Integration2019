#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/allen_brain
DLNAME="allen_brain.zip"

if [ $1 = "orig" ]; then
	URL="http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/uwrs4k5tc77tqkv/allen_brain.zip?dl=0"
else
	echo "bad option"; exit;
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

unzip $RD/allen_brain.zip -d $RD
rm $RD/allen_brain.zip
