args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

bipolar.integrated <- readRDS(file = args[1])
DefaultAssay(object = bipolar.integrated) <- "RNA"
bipolar.integrated[["integrated"]] <- NULL

genes.use <- bipolar.integrated@misc$integration.features

all.bipolar <- SplitObject(object = bipolar.integrated, split.by = "replicate")
bipolar.counts <- list()
for(i in 1:length(x = all.bipolar)) {
  bipolar.counts[[i]] <- GetAssayData(object = all.bipolar[[i]], slot = "counts")
}

library(devtools)
detach("package:Seurat", unload = TRUE)
library(Seurat)
bipolar.v2 <- list()
for(i in 1:length(x = all.bipolar)) {
  bipolar.v2[[i]] <- CreateSeuratObject(raw.data = bipolar.counts[[i]], meta.data = all.bipolar[[i]]@meta.data)
  bipolar.v2[[i]] <- NormalizeData(object = bipolar.v2[[i]])
  bipolar.v2[[i]] <- ScaleData(object = bipolar.v2[[i]], genes.use = genes.use)
}

dims <- 30

bipolar.integrated <- RunMultiCCA(object.list = bipolar.v2, genes.use = genes.use, num.ccs = dims)


# take through clustering
bipolar.integrated <- AlignSubspace(bipolar.integrated, reduction.type = "cca", grouping.var = "replicate", dims.align = 1:dims)
bipolar.integrated <- FindClusters(bipolar.integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:dims, save.SNN = T, resolution = 0.4)

#bipolar.integrated <- RunTSNE(bipolar.integrated, reduction.use = "cca.aligned", dims.use = 1:dims, tsne.method = "FIt-SNE")
bipolar.integrated <- RunUMAP(bipolar.integrated, reduction.use = "cca.aligned", dims.use = 1:dims)

detach("package:Seurat", unload = TRUE)
devtools::load_all(args[3])
bipolar.integrated <- UpdateSeuratObject(object = bipolar.integrated)

saveRDS(object = bipolar.integrated, file = args[2])
