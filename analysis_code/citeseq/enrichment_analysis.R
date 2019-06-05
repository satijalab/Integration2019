library(org.Hs.eg.db)
library(GOstats)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)

# thanks to Dave Tang for tutorial: https://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/
# Note RSQLite version 2.1.1 causes a big in GOstats: https://www.biostars.org/p/316760/

# load CD69 marker genes
markers <- read.table("analysis_data/cd69_markers.tsv", sep = "\t")
cd69.cite <- read.table("analysis_data/cd69_high_moran_citeseq.tsv", sep = "\t")
cd69.hca <- read.table("analysis_data/cd69_high_moran_hca.tsv", sep = "\t")

# make go term object
entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object)
entrez_to_go <- as.list(entrez_object[mapped_genes])
go_object <- as.list(org.Hs.egGO2EG)

run_go <- function(genes, ontology){
  geneSymbols <- mapIds(org.Hs.eg.db, keys=as.character(genes),
                        column="ENTREZID", keytype="SYMBOL", multiVals="first")
  params <- new(
    Class = 'GOHyperGParams',
    geneIds = as.character(geneSymbols),
    universeGeneIds = mapped_genes,
    ontology = ontology,
    pvalueCutoff = 0.001,
    conditional = FALSE,
    testDirection = 'over',
    annotation = "org.Hs.eg.db"
  )
  hgOver <- hyperGTest(params)
  results <- summary(hgOver)
  return(results)
}

bp <- run_go(rownames(markers), 'BP')
mf <- run_go(rownames(markers), 'MF')

bp.cite <- run_go(cd69.cite$gene, 'BP')
mf.cite <- run_go(cd69.cite$gene, 'MF')

bp.hca <- run_go(cd69.hca$gene, 'BP')
mf.hca <- run_go(cd69.hca$gene, 'MF')

write.table(bp, 'analysis_data/go_terms_bp_cd69.tsv', quote = FALSE, sep = '\t')
write.table(mf, 'analysis_data/go_terms_mf_cd69.tsv', quote = FALSE, sep = '\t')

ggplot(head(mf, 5), aes(reorder(Term, -Pvalue), -log2(Pvalue))) +
  geom_bar(stat = 'identity') +
  xlab("Molecular function") +
  coord_flip() +
  theme_classic(base_size = 10) +
  ggsave("figures/citeseq/go_terms_cd69_mf.pdf",
         height = 2, width = 4)

ggplot(head(bp, 5), aes(reorder(Term, -Pvalue), -log2(Pvalue))) +
  geom_bar(stat = 'identity') +
  xlab("Biological process") +
  coord_flip() +
  theme_classic(base_size = 10) +
  ggsave("figures/citeseq/go_terms_cd69_bp.pdf",
         height = 2, width = 4)

# save top GO terms for supplementary table
write.table(head(bp, 30), 'analysis_data/top_go_terms_bp_cd69.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(head(mf, 30), 'analysis_data/top_go_terms_mf_cd69.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
