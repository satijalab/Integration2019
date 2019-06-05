#!/bin/bash

RD=raw_data/immune/citeseq/

URL_RNA_A="https://www.dropbox.com/s/nyoqrrl3ilqvo65/raw_gene_bc_matrices_h5_mnca.h5?dl=0"
URL_RNA_B="https://www.dropbox.com/s/buksfjyant9o7u8/raw_gene_bc_matrices_h5_mncb.h5?dl=0"
URL_HTO_A="https://www.dropbox.com/s/l3rq5bfikjp7eph/hto_demux_mnc_a.tsv?dl=0"
URL_HTO_B="https://www.dropbox.com/s/elveow6hwccdaxt/hto_demux_mnc_b.tsv?dl=0"
URL_ADT_A="https://www.dropbox.com/s/7hc1h1e3u54bxwr/adt_demux_mnc_a.tsv?dl=0"
URL_ADT_B="https://www.dropbox.com/s/zcagxdqg0lnc64k/adt_demux_mnc_b.tsv?dl=0"

curl -L $URL_RNA_A --create-dirs -o $RD/raw_gene_bc_matrices_h5_mnca.h5
curl -L $URL_RNA_B --create-dirs -o $RD/raw_gene_bc_matrices_h5_mncb.h5

curl -L $URL_HTO_A --create-dirs -o $RD/hto_demux_mnc_a.tsv
curl -L $URL_HTO_B --create-dirs -o $RD/hto_demux_mnc_b.tsv

curl -L $URL_ADT_A --create-dirs -o $RD/adt_demux_mnc_a.tsv
curl -L $URL_ADT_B --create-dirs -o $RD/adt_demux_mnc_b.tsv
