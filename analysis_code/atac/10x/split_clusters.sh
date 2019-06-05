#! /bin/bash

for filename in analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/*.txt; do
	echo $filename
	outf=${filename%.txt}
	sinto filterbarcodes -b raw_data/10x_atac/atac_v1_pbmc_10k_possorted_bam.bam \
		-c $filename \
		-o "${outf}.bam" \
		-p 8 \
		-m tag
	samtools index -@ 8 "${outf}.bam"
	bamCoverage --bam "${outf}.bam" --binSize 1 --numberOfProcessors 8 --normalizeUsing RPKM -o "${outf}.bw"
done

bamCoverage --bam raw_data/10x_atac/atac_v1_pbmc_10k_possorted_bam.bam --binSize 1 --numberOfProcessors 8 --normalizeUsing RPKM -o analysis_data/atac/10x/celltypes/single_cell/pbmc_all.bw