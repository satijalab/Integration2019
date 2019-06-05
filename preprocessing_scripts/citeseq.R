suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

dir.create("figures/citeseq/")

# load data
rna_a <- Read10X_h5("raw_data/immune/citeseq/raw_gene_bc_matrices_h5_mnca.h5")
rna_b <- Read10X_h5("raw_data/immune/citeseq/raw_gene_bc_matrices_h5_mncb.h5")
hto_a <- read.table("raw_data/immune/citeseq/hto_demux_mnc_a.tsv")
hto_b <- read.table("raw_data/immune/citeseq/hto_demux_mnc_b.tsv")
adt_a <- read.table("raw_data/immune/citeseq/adt_demux_mnc_a.tsv")
adt_b <- read.table("raw_data/immune/citeseq/adt_demux_mnc_b.tsv")

rownames(x = rna_a) <- make.unique(names = rownames(x = rna_a))
rownames(x = rna_b) <- make.unique(names = rownames(x = rna_b))

process_cite <- function(rna, hto, adt, orig.ident) {
  # filter cells based on RNA counts
  rownames(x = hto) <- paste0(rownames(x = hto), "-1")
  rownames(x = adt) <- paste0(rownames(x = adt), "-1")
  cellcounts <- Matrix::colSums(x = rna)
  rna <- rna[, cellcounts > 1000]
  hto <- hto[colnames(x = rna), ]
  
  # create object, add HTO assay
  citeseq <- CreateSeuratObject(counts = rna, project = orig.ident, min.cells = 10)
  citeseq[['HTO']] <- CreateAssayObject(counts = t(hto[, 1:10]))
  DefaultAssay(object = citeseq) <- 'HTO'
  citeseq <- NormalizeData(object = citeseq, normalization.method = 'CLR', across = "cells")
  citeseq <- HTODemux(object = citeseq, positive.quantile = 0.99, nsamples = 100)

  # process RNA data
  DefaultAssay(object = citeseq) <- 'RNA'
  Idents(object = citeseq) <- "HTO_classification.global"
  citeseq <- SubsetData(
    object = citeseq, 
    cells = WhichCells(object = citeseq, idents = c('Negative', 'Doublet'), invert = TRUE)
  )
  citeseq <- NormalizeData(object = citeseq)
  citeseq <- FindVariableFeatures(object = citeseq, nfeatures = 2000)
  citeseq <- ScaleData(object = citeseq)

  # process ADT data
  adt <- adt[colnames(x = citeseq), 1:(ncol(x = adt)-3)]
  citeseq[['ADT']] <- CreateAssayObject(counts = t(x = adt))
  DefaultAssay(object = citeseq) <- 'ADT'
  citeseq <- NormalizeData(object = citeseq, normalization.method = 'CLR', across = 'cells')
  return(list(citeseq))
}

rep1 <- process_cite(rna_a, hto_a, adt_a, orig.ident = 'batch1')
rep2 <- process_cite(rna_b, hto_b, adt_b, orig.ident = 'batch2')

# merge objects
citeseq <- merge(x = rep1[[1]], y = rep2[[1]], add.cell.ids = c("a", "b"))
DefaultAssay(object = citeseq) <- "RNA"
citeseq <- NormalizeData(object = citeseq, across = "cells")
citeseq <- FindVariableFeatures(object = citeseq, nfeatures = 2000)
citeseq <- ScaleData(object = citeseq)
citeseq <- RunPCA(object = citeseq, features = VariableFeatures(object = citeseq), npcs = 60)
citeseq <- RunUMAP(object = citeseq, dims = 1:60)
citeseq <- RunTSNE(object = citeseq, dims = 1:60, tsne.method = 'FIt-SNE')
citeseq <- FindNeighbors(object = citeseq, assay = "RNA", reduction = "pca", dims = 1:60)
citeseq <- FindClusters(object = citeseq, resolution = 0.5)
markers <- FindAllMarkers(object = citeseq, features = VariableFeatures(object = citeseq), logfc.threshold = 0.5)
citeseq@misc <- list('markers_res_0.5' = markers)

DefaultAssay(object = citeseq) <- "ADT"
p <- FeaturePlot(
  object = citeseq,
  features = rownames(x = citeseq),
  pt.size = 0.1,
  min.cutoff = 'q10',
  max.cutoff = 'q99',
  combine = FALSE,
  cols = c("grey", "darkgreen")
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + cowplot::theme_cowplot(font_size = 30) + theme(legend.position = 'none')
}
plots <- CombinePlots(p, ncol = 5, legend = 'none')
save_plot(plot = plots, filename = "figures/citeseq/citeseq_featureplot.png", base_height = 30, base_width = 30, dpi = 200)
saveRDS(object = citeseq, file = args[2])
