#! bin/sh

# map RNA reads with cellranger v2.1.0
cellranger count --id=MNC-A \
--fastqs=MNC-A \
--transcriptome=refdata-cellranger-GRCh38-1.2.0 \
--sample=MNC-A --localcores=10 --localmem=50

cellranger count --id=MNC-B \
--fastqs=MNC-B \
--transcriptome=refdata-cellranger-GRCh38-1.2.0 \
--sample=MNC-B --localcores=10 --localmem=50

# count ADTs and HTOs with CITE-seq-count v1.2
python CITE-seq-count.py -R1 MNC-A-HTO_R1.fastq.gz -R2 MNC-A-HTO_R2.fastq.gz -cbf 1 -cbl 16 -umif 17 -umil 26 -hd 1 \
-tr "^[ATGC]{15}" -o hto_demux_mnc_a.tsv -t hto_mnc.csv

python CITE-seq-count.py -R1 MNC-B-HTO_R1.fastq.gz -R2 MNC-B-HTO_R2.fastq.gz -cbf 1 -cbl 16 -umif 17 -umil 26 -hd 1 \
-tr "^[ATGC]{15}" -o hto_demux_mnc_b.tsv -t hto_mnc.csv

python CITE-seq-count.py -R1 MNC-A-ADT_R1.fastq.gz -R2 MNC-A-ADT_R2.fastq.gz -cbf 1 -cbl 16 -umif 17 -umil 26 -hd 1 \
-tr "^[ATGC]{15}" -o adt_demux_mnc_a.tsv -t adt_mnc.csv

python CITE-seq-count.py -R1 MNC-B-ADT_R1.fastq.gz -R2 MNC-B-ADT_R2.fastq.gz -cbf 1 -cbl 16 -umif 17 -umil 26 -hd 1 \
-tr "^[ATGC]{15}" -o adt_demux_mnc_b.tsv -t adt_mnc.csv

