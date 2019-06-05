load_all("~/seurat/")

celseq <- readRDS("seurat_objects/celseq.rds")
celseq2 <- readRDS("seurat_objects/celseq2.rds")
fluidigmc1 <- readRDS("seurat_objects/fluidigmc1.rds")
smartseq2 <- readRDS("seurat_objects/smartseq2.rds")
indrop <- readRDS("seurat_objects/inDrop.rds")
all.pancreas <- c(list(celseq, celseq2, smartseq2, fluidigmc1), indrop)

genes.use <- SelectIntegrationFeatures(obj.list = all.pancreas, nfeatures = 2000)
for(i in 1:length(x = all.pancreas)) {
  all.pancreas[[i]] <- ScaleData(object = all.pancreas[[i]], features = genes.use, verbose = FALSE)
  all.pancreas[[i]] <- RunPCA(object = all.pancreas[[i]], compute.dims = 40, features = genes.use, verbose = FALSE)
}

pancreas.integrated <- MultiIntegrateData(
  obj.list = all.pancreas, 
  assay = "RNA",
  reduction = "cca.cosine", 
  k.nn = 300,
  k.mnn = 5, 
  sd = 0.1, 
  do.cpp = TRUE, 
  verbose = TRUE, 
  dims = 1:30,
  features = genes.use,
  features.integrate = genes.use,
  eps=1
)

# take through clustering and annotate
DefaultAssay(object = pancreas.integrated) <- "integrated"
VariableFeatures(pancreas.integrated) <- genes.use
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated,
                              features = VariableFeatures(pancreas.integrated),
                              verbose = FALSE, compute.dims = 30)
pancreas.integrated <- RunTSNE(object = pancreas.integrated,
                               reduction = "pca",
                               dims = 1:30, tsne.method = "FIt-SNE")
pancreas.integrated <- RunUMAP(object = pancreas.integrated,
                               reduction = "pca",
                               dims = 1:30)
pancreas.integrated <- BuildSNN(object = pancreas.integrated, reduction = "pca", dims = 1:30,k.param = 10)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 1.5)

#quick annotation
initial_ID <- unlist(tapply(pancreas.integrated@meta.data$assigned_cluster,as.numeric(as.character(Idents(pancreas.integrated))),function(x)names(sort(table(x),decreasing = T))[1]))
new_ID <- (initial_ID[as.character(Idents(pancreas.integrated))])
pancreas.integrated$joint_clusterID <- new_ID
pancreas.integrated@meta.data$joint_clusterID[is.na(pancreas.integrated@meta.data$newID)]="acinar"

Idents(pancreas.integrated) <- "joint_clusterID"

