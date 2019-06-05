suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

pancreas.integrated <- readRDS(file = args[1])
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated[["integrated"]] <- NULL
genes.use <- pancreas.integrated@misc$integration.features

all.pancreas <- SplitObject(object = pancreas.integrated, split.by = "replicate")
pancreas.counts <- list()
for(i in 1:length(x = all.pancreas)) {
  pancreas.counts[[i]] <- GetAssayData(object = all.pancreas[[i]], slot = "counts")
}

library(devtools)
detach("package:Seurat", unload = TRUE)
library(Seurat)
pancreas.v2 <- list()
for(i in 1:length(x = all.pancreas)) {
  pancreas.v2[[i]] <- CreateSeuratObject(raw.data = pancreas.counts[[i]], meta.data = all.pancreas[[i]]@meta.data)
  pancreas.v2[[i]] <- NormalizeData(object = pancreas.v2[[i]])
  pancreas.v2[[i]] <- ScaleData(object = pancreas.v2[[i]], genes.use = genes.use)
}

dims <- 30

pancreas.integrated <- RunMultiCCA(object.list = pancreas.v2, genes.use = genes.use, num.ccs = 30)


# take through clustering
pancreas.integrated <- AlignSubspace(pancreas.integrated, reduction.type = "cca", grouping.var = "replicate", dims.align = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:30, save.SNN = T, resolution = 0.4)

#pancreas.integrated <- RunTSNE(pancreas.integrated, reduction.use = "cca.aligned", dims.use = 1:30)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction.use = "cca.aligned", dims.use = 1:30)

detach("package:Seurat", unload = TRUE)
devtools::load_all(args[3])
pancreas.integrated <- UpdateSeuratObject(object = pancreas.integrated)

saveRDS(object = pancreas.integrated, file = args[2])
