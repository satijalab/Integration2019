#!/bin/bash

# $1 download source - orig versus dropbox 

RD=raw_data/spatial
DLNAME="osm_fish.loom"

if [ $1 = "orig" ]; then
	URL="http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom"
elif [ $1 == "dropbox"]; then
	URL="https://www.dropbox.com/s/py2jtgonubflk89/osm_fish.loom?dl=0"
else
	echo "bad option"; exit;	
fi

curl -L $URL --create-dirs -o $RD/$DLNAME

