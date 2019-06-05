curl http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5 \
--create-dirs -o raw_data/10x_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5

curl http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv \
--create-dirs -o raw_data/10x_atac/atac_v1_pbmc_10k_singlecell.csv

curl http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_possorted_bam.bam \
--create-dirs -o raw_data/10x_atac/atac_v1_pbmc_10k_possorted_bam.bam

curl http://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_possorted_bam.bam.bai \
--create-dirs -o raw_data/10x_atac/atac_v1_pbmc_10k_possorted_bam.bam.bai

curl ftp://ftp.ensembl.org/pub/grch37/release-82/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz \
--create-dirs -o raw_data/10x_atac/Homo_sapiens.GRCh37.82.gtf

curl -L https://www.dropbox.com/s/iljanfyjl99o4xi/clustery.txt?dl=0 \
--create-dirs -o raw_data/10x_atac/clustery.txt

curl -L https://www.dropbox.com/s/hfnrplpvwu0qrk0/clusterx.txt?dl=0 \
--create-dirs -o raw_data/10x_atac/clusterx.txt