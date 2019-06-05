library(DESeq2)
library(Seurat)
library(ggplot2)
library(dplyr)

sample_design <- c("CD69-A", "ACAGTGCTCATG", "CD69.negative",
       "CD69-B", "GAACACCATAGC", "CD69.negative",
       "CD69-C", "GTGTAGTGTGGA", "CD69.negative",
       "CD69-D", "GGATAGTTAGCC", "CD69.negative",
       "CD69+A", "GTTGAGCTTGAC", "CD69.positive",
       "CD69+B", "TAAGCGTGGAAC", "CD69.positive",
       "CD69+C", "CAGAGTCGACTT", "CD69.positive",
       "CD69+D", "GTTCCACCTGTA", "CD69.positive")

sample_design <- data.frame(matrix(sample_design, ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)
colnames(sample_design) <- c("Sample", "Sequence", "Group")
rownames(sample_design) <- sample_design$Sample
counts <- read.table("raw_data/bone_marrow/CD69.umi.txt",
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
rownames(counts) <- counts$GENE
counts$GENE <- NULL
counts <- counts[, sample_design$Sequence]
counts <- as.matrix(counts)
colnames(counts) <- rownames(sample_design)

de.genes <- read.table('analysis_data/cd69_markers.tsv', sep = '\t')
common.genes <- rownames(counts)[rownames(counts) %in% rownames(de.genes)]
rownames(sample_design) <- sample_design$Sample

dd <- DESeqDataSetFromMatrix(counts, sample_design, ~ Group)
dd <- dd[ rowSums(counts(dd)) > 1, ]
dd <- DESeq(dd)
results <- results(dd, alpha=0.01)
results$gene <- rownames(results)
results <- results[!(is.na(results$padj)),]

results <- as.data.frame(results) %>% 
  rowwise() %>% 
  mutate(significant = ifelse(abs(log2FoldChange) > 1.5 & padj < 0.01, TRUE, FALSE),
         single.cell = gene %in% rownames(de.genes))

highlight.genes <- c('IFNG', "CCL4L2", "CCL3", "CCL4", "XCL2")

ggplot(results, aes(log2FoldChange, -log2(padj), color = single.cell)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c('black', 'red')) +
  geom_text(aes(label=ifelse(gene %in% highlight.genes, gene, '')), hjust=0, vjust=0, size = 2) +
  theme_classic(base_size = 6) +
  theme(legend.position = 'none') + 
  ylab("-log2(qvalue)") + 
  xlab("log2 fold change expression") +
  ggsave("figures/citeseq/facs_volcano.png", height = 5, width = 5, units = 'cm', dpi = 800)

write.table(results, "analysis_data/facs_cd69_de.tsv", sep = '\t', quote = FALSE)
