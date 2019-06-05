bedtools multicov -bams \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/Bcellprogenitor.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/pre-Bcell.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Naive.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD8effector.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/DoublenegativeTcell.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD8Naive.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/NKdim.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD14Monocytes.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD16Monocytes.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/Dendriticcell.bam \
	analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/pDC.bam \
	-bed raw_data/corces_atac/peaks.bed \
	> analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/celltype_coverage.bed
