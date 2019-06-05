# umaps for fig2

library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[6])

# args <- c("~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_seuratV3.rds",
#           "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_seuratV2.rds",
#           "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_mnnCorrect.rds",
#           "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_scanorama.rds")

pi.v3 <- readRDS(file = args[1])
pi.v2 <- readRDS(file = args[2])
pi.mnn <- readRDS(file = args[3])
pi.sc <- readRDS(file = args[4])

pi.v3$celltype <- factor(x = pi.v3$celltype)
pi.v2$celltype <- factor(x = pi.v2$celltype, levels = levels(x = pi.v3$celltype))
pi.mnn$celltype <- factor(x = pi.mnn$celltype, levels = levels(x = pi.v3$celltype))
pi.sc$celltype <- factor(x = pi.sc$celltype, levels = levels(x = pi.v3$celltype))

pi.v3$replicate <- factor(x = pi.v3$replicate)
pi.v2$replicate <- factor(x = pi.v2$replicate, levels = levels(x = pi.v3$replicate))
pi.mnn$replicate <- factor(x = pi.mnn$replicate, levels = levels(x = pi.v3$replicate))
pi.sc$replicate <- factor(x = pi.sc$replicate, levels = levels(x = pi.v3$replicate))

dptheme <- NoLegend() +
  theme(legend.position = "top", legend.title = element_blank(), legend.justification = "center")

pa <- DimPlot(object = pi.v3, reduction = "umap", group.by = "replicate") +
  dptheme

colnames(x = pi.v2[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

pb <- DimPlot(object = pi.v2, reduction = "umap", group.by = "replicate") +
  NoLegend()

pc <- DimPlot(object = pi.mnn, reduction = "umap", group.by = "replicate") +
  NoLegend()
  
pd <- DimPlot(object = pi.sc, reduction = "umap", group.by = "replicate") + 
  NoLegend()

pe <- DimPlot(object = pi.v3, reduction = "umap", group.by = "celltype") +
  dptheme

pf <- DimPlot(object = pi.v2, reduction = "umap", group.by = "celltype") +
  NoLegend()

pg <- DimPlot(object = pi.mnn, reduction = "umap", group.by = "celltype") +
  NoLegend()

ph <- DimPlot(object = pi.sc, reduction = "umap", group.by = "celltype") + 
  NoLegend()

save(pa, pb, pc, pd, pe, pf, pg, ph, file = args[5])
