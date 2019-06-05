suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(future))
args <- commandArgs(trailingOnly = TRUE)

#args <- paste0("~/Projects/muir/seurat_objects/mca",c("_tm_droplet.rds","_tm_facs.rds"))
tm_droplet <- readRDS(file = args[1])
tm_facs <- readRDS(file = args[2])

plan(strategy = 'multicore', workers = 6)
options(future.globals.maxSize = 10^12)


tm_droplet[["tissue"]] <- tolower(x = tm_droplet$tissue)
tm_facs[["tissue"]] <- tolower(x = tm_facs$tissue)
tm_facs$tissue <- as.character(tm_facs$tissue)
tm_facs$tissue[names(which(is.na(tm_facs$tissue)))] <- "None"
tm_facs <- subset(tm_facs,tissue != 'None')
tissue_all <- Reduce(
  f = intersect, 
  x = list(
    unique(x = tm_droplet$tissue), 
    unique(x = tm_facs$tissue)
  )
)

ob.list <- list(tm_droplet, tm_facs)
for(i in 1:2) {
  ob.list[[i]] <- NormalizeData(object = ob.list[[i]], verbose = FALSE)
  ob.list[[i]] <- FindVariableFeatures(object = ob.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
}

dimensions <- 100
mca.anchors <- FindIntegrationAnchors(
  object.list = ob.list,
  dims = 1:dimensions,
  eps = 1
)

mca <- IntegrateData(
  anchorset = mca.anchors, 
  dims = 1:dimensions, 
  eps = 1
)

DefaultAssay(object = mca) <- 'integrated'
genes.use <- VariableFeatures(object = mca)
mca <- ScaleData(object = mca, features = genes.use)
mca <- RunPCA(object = mca, features = genes.use, npcs = dimensions)
mca <- FindNeighbors(object = mca, dims = 1:dimensions, k.param = 20)
mca <- FindClusters(object = mca, resolution = 4)
mca <- RunTSNE(object = mca, dims = 1:dimensions, tsne.method = "FIt-SNE")
mca <- RunUMAP(object = mca, dims = 1:dimensions)
saveRDS(object = mca, file = args[3])
