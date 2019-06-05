#! /bin/bash

for filename in analysis_data/atac/clusters/*.txt; do
	echo $filename
	outf=${filename%.txt}
	sinto filterbarcodes -b raw_data/atac/PreFrontalCortex_62216.bam \
		-c $filename \
		-o "${outf}.bam" \
		-p 10 \
		-m readname
	samtools index "${outf}.bam"
	bamCoverage --bam "${outf}.bam" --binSize 1 --numberOfProcessors 10 --normalizeUsing RPKM -o "${outf}.bw"
done
bamCoverage --bam raw_data/atac/PreFrontalCortex_62216.bam --binSize 1 --numberOfProcessors 10 --normalizeUsing RPKM -o analysis_data/atac/clusters/pfc_all.bw
