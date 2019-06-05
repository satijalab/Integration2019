library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

allen <- read.csv("raw_data/allen_brain/VIS_gene_expression_matrix_2016-10-27/fpkm_table.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen$gene_id...analysis_run_id <- NULL
allen <- as.matrix(allen)
genes <- read.table("raw_data/allen_brain/VIS_gene_expression_matrix_2016-10-27/rows-genes.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(allen) <- make.unique(genes$gene_symbol)
meta.data <- read.csv("raw_data/allen_brain/VIS_gene_expression_matrix_2016-10-27/columns-samples.csv",
                      row.names=1, stringsAsFactors=FALSE)
colnames(allen) <- substring(colnames(allen), 2)

al <- CreateSeuratObject(
  counts = allen,
  project = 'VISp2016',
  meta.data = meta.data,
  min.cells = 10,
  min.features = 200
)
al <- NormalizeData(al)
al <- FindVariableFeatures(al, nfeatures = 2000, selection.method = 'dispersion')
al <- ScaleData(al)
al <- RunPCA(al, npcs = 50)
al <- RunUMAP(al, dims = 1:50)
saveRDS(al, file = args[1])
