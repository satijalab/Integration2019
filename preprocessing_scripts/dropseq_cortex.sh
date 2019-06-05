#!/bin/bash

RD=raw_data/dropseq_cortex

curl https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz --create-dirs -o $RD/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz
curl https://storage.googleapis.com/dropviz-downloads/static/regions/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.cell_cluster_outcomes.RDS --create-dirs -o $RD/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.cell_cluster_outcomes.RDS
