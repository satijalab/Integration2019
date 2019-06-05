#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/spatial/starmap
DLNAME="starmap.zip"

if [ $1 = "orig" ]; then
	URL="https://www.dropbox.com/sh/f7ebheru1lbz91s/AABYSSjSTppBmVmWl2H4s_K-a?dl=0"
elif [ $1 == "dropbox"]; then
	URL=""
else
	echo "bad option"; exit;
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

unzip $RD/starmap.zip -d $RD
rm $RD/starmap.zip
