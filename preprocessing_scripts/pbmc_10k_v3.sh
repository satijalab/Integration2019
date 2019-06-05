curl http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5 \
--create-dirs -o raw_data/10x_scrna/pbmc10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5

curl -L https://www.dropbox.com/s/m390vmhnabti4bl/matrix_doublets.tsv?dl=0 \
--create-dirs -o raw_data/10x_scrna/pbmc10k_v3/matrix_doublets.tsv