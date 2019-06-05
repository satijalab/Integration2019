suppressMessages(library(Seurat))
suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)

input.file <- paste0(getwd(), "/raw_data/pbmc/pbmc33k/filtered_gene_bc_matrices/hg19/")

raw.data <- Read10X(data.dir = input.file)
seurat.object <- CreateSeuratObject(counts = raw.data)
seurat.object <- seurat.object[, WhichCells(object = seurat.object, expression = "nCount_RNA > 1000")]
seurat.object <- NormalizeData(object = seurat.object)
seurat.object <- FindVariableFeatures(object = seurat.object)
seurat.object <- ScaleData(object = seurat.object, features = VariableFeatures(object = seurat.object))
seurat.object <- RunPCA(object = seurat.object, npcs = 25, verbose = FALSE)
seurat.object <- RunUMAP(object = seurat.object, dims = 1:25)
seurat.object <- RunTSNE(object = seurat.object, dims = 1:25, tsne.method = "FIt-SNE")
seurat.object <- BuildSNN(seurat.object, dims = 1:25, k.param = 20)
seurat.object <- FindClusters(seurat.object, resolution = 1.2)

#remove doublet populations, cycling cells, and IFN+ T cells 
seurat.object <- seurat.object[, WhichCells(object = seurat.object, idents = c(18,19,15,12,14), invert = TRUE)]
seurat.object <- ScaleData(object = seurat.object, features = VariableFeatures(object = seurat.object))
seurat.object <- RunPCA(object = seurat.object, npcs = 25)
seurat.object <- RenameIdents(
  object = seurat.object, 
  '0' = 'CD4_Memory', '1' = 'CD4_Naive', '2' = "Mono_CD14", '3' = "B_Pre", '4' = "Mono_CD15", 
  '5' = "CD8_Memory", '6' = "CD8_Naive", '7' = "NK_Dim", '8' = "Mono_CD16", '9' = "CD8_Effector", 
  '10' = 'B_Pro', '11' = 'DC', '13' = 'NK_Bright', '16' = "Mk", '17' = "pDC"
)
seurat.object$cluster <- Idents(object = seurat.object)

# Save rds file
saveRDS(object = seurat.object, file = args[1])
