library(Seurat)
set.seed(1234)

# load atac peak counts and create gene activity matrix
peaks <- Read10X_h5("raw_data/10x_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
activity.matrix <- CreateGeneActivityMatrix(
  peak.matrix = peaks,
  annotation.file = "raw_data/10x_atac/Homo_sapiens.GRCh37.82.gtf"
)

# create object and filter cells
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = 'ATAC', project = '10x_ATAC')
pbmc.atac[['RNA']] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("raw_data/10x_atac/atac_v1_pbmc_10k_singlecell.csv", sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)

# process activity score matrix
DefaultAssay(pbmc.atac) <- 'RNA'
pbmc.atac <- FindVariableFeatures(pbmc.atac, nfeatures = 5000)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.atac <- RunPCA(pbmc.atac, npcs = 50)
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = 'pca', dims = 1:30)
pbmc.atac <- RunUMAP(pbmc.atac, graph = 'RNA_nn')

# process peaks
DefaultAssay(pbmc.atac) <- "ATAC"
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))
pbmc.atac <- RunLSI(pbmc.atac, n = 100)
pbmc.atac <- FindNeighbors(pbmc.atac, dims = 1:30, reduction = 'lsi')
pbmc.atac <- FindClusters(pbmc.atac, graph = 'ATAC_nn', resoluton = 1)
pbmc.atac <- RunTSNE(pbmc.atac, reduction = 'lsi', dims = 1:30, reduction.name = 'tsne.lsi', reduction.key = 'LSItSNE_')
pbmc.atac <- RunUMAP(pbmc.atac, graph = 'ATAC_nn', reduction.name = 'umap.lsi', reduction.key = 'LSIUMAP_')
saveRDS(pbmc.atac, "seurat_objects/pbmc_10k_atac.rds")
