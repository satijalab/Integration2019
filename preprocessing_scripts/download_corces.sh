#!/bin/bash

RD=raw_data/corces_atac

curl https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Mono-7.merge.s20.w150sw.bw --create-dirs -o $RD/bigwig/Mono-7.merge.s20.w150sw.bw
curl https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/CD4-9.merge.s20.w150sw.bw --create-dirs -o $RD/bigwig/CD4-9.merge.s20.w150sw.bw
curl https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/CD8-10.merge.s20.w150sw.bw --create-dirs -o $RD/bigwig/CD8-10.merge.s20.w150sw.bw
curl https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Nkcell-11.merge.s20.w150sw.bw --create-dirs -o $RD/bigwig/Nkcell-11.merge.s20.w150sw.bw
curl https://s3-us-west-1.amazonaws.com/ryancorces-heme/hg19/Bcell-13.merge.s20.w150sw.bw --create-dirs -o $RD/bigwig//Bcell-13.merge.s20.w150sw.bw
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/suppl/GSE74912_ATACseq_All_Counts.txt.gz --create-dirs -o $RD/GSE74912_ATACseq_All_Counts.txt.gz
gzip -d $RD/GSE74912_ATACseq_All_Counts.txt.gz
sed 1d $RD/GSE74912_ATACseq_All_Counts.txt | cut -f1,2,3 > $RD/peaks.bed
