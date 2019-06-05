suppressMessages(library(Seurat))
suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)

input.file <- paste0(getwd(), "/raw_data/pbmc/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc.data <- Read10X(data.dir = input.file)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = pbmc, slot = "counts")[mito.genes, ]) /
  Matrix::colSums(GetAssayData(object = pbmc, slot = "counts"))
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc <- subset(x = pbmc, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.05)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize")
pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- ProjectDim(object = pbmc, reduction = "pca", verbose = FALSE)
pbmc <- BuildSNN(object = pbmc, reduction = "pca", dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.4, verbose = FALSE)
pbmc <- RunTSNE(object = pbmc, dims = 1:10)
pbmc <- RenameIdents(object = pbmc, "0" = "CD4 T cells", "1" = "CD14+ Monocytes", "2" = "B cells", 
                     "3" = "CD8 T cells", "4" = "FCGR3A+ Monocytes", "5" = "NK cells", 
                     "6" = "Dendritic cells", "7" = "Megakaryocytes")

pbmc$cluster <- Idents(object = pbmc)

saveRDS(object = pbmc, file = args[1])

