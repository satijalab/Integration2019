#!/bin/bash

# $1 download source - orig versus dropbox

RD=raw_data/atac
DLNAME1="activity_scores.quantitative.rds"
DLNAME2="cell_metadata.txt"
DLNAME3="atac_matrix.binary.qc_filtered.rds"
DLNAME4="PreFrontalCortex_62216.bam"
DLNAME5="PreFrontalCortex_62216.bam.bai"


if [ $1 = "orig" ]; then
	URL1="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/activity_score_matrices/activity_scores.quantitative.rds"
	URL2="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_metadata.txt"
	URL3="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/matrices/atac_matrix.binary.qc_filtered.rds"
	URL4="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/bams/PreFrontalCortex_62216.bam"
	URL5="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/bams/PreFrontalCortex_62216.bam.bai"
elif [ $1 == "dropbox"]; then
	URL1="https://www.dropbox.com/s/7rasn05clfhqzx9/activity_scores.quantitative.rds?dl=0"
	URL2="https://www.dropbox.com/s/4nc5j92sbgd23st/atac_pfc_cell_metadata.txt?dl=0"
	URL3="https://www.dropbox.com/s/4lu2ltp3nudp0hw/atac_matrix.binary.qc_filtered.rds?dl=0"
else
	echo "bad option"; exit;
fi

curl -L $URL1 --create-dirs -o $RD/$DLNAME1
curl -L $URL2 --create-dirs -o $RD/$DLNAME2
curl -L $URL3 --create-dirs -o $RD/$DLNAME3
curl -L $URL4 --create-dirs -o $RD/$DLNAME4
curl -L $URL5 --create-dirs -o $RD/$DLNAME5
