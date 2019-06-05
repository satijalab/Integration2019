library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(Seurat, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_seuratV3.rds",
#           "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_seuratV2.rds",
#           "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_mnnCorrect.rds",
#           "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_scanorama.rds",
#           "~/Projects/muir/figures/supplementary/bipolar_holdout_tsnes.pdf")

bi.v3 <- readRDS(file = args[1])
bi.v2 <- readRDS(file = args[2])
bi.mnn <- readRDS(file = args[3])
bi.sc <- readRDS(file = args[4])

bi.v3$celltype <- factor(x = bi.v3$celltype)
bi.v2$celltype <- factor(x = bi.v2$celltype, levels = levels(x = bi.v3$celltype))
bi.mnn$celltype <- factor(x = bi.mnn$celltype, levels = levels(x = bi.v3$celltype))
bi.sc$celltype <- factor(x = bi.sc$celltype, levels = levels(x = bi.v3$celltype))

bi.v3$replicate <- factor(x = bi.v3$replicate)
bi.v2$replicate <- factor(x = bi.v2$replicate, levels = levels(x = bi.v3$replicate))
bi.mnn$replicate <- factor(x = bi.mnn$replicate, levels = levels(x = bi.v3$replicate))
bi.sc$replicate <- factor(x = bi.sc$replicate, levels = levels(x = bi.v3$replicate))

pa <- DimPlot(object = bi.v3, reduction = "tsne", group.by = "replicate", cells = sample(Cells(bi.v3))) +
  ggtitle(label = "Seurat V3") + labs(color = 'Replicate') + theme(legend.title = element_text(size = 18))

colnames(x = bi.v2[["tsne"]]@cell.embeddings) <- paste0("tSNE_", 1:2)

pb <- DimPlot(object = bi.v2, reduction = "tsne", group.by = "replicate", cells = sample(Cells(bi.v2))) +
  ggtitle(label = "Seurat V2") +
  NoLegend()

pc <- DimPlot(object = bi.mnn, reduction = "tsne", group.by = "replicate", cells = sample(Cells(bi.mnn))) +
  ggtitle(label = "mnnCorrect") + 
  NoLegend()

pd <- DimPlot(object = bi.sc, reduction = "tsne", group.by = "replicate", cells = sample(Cells(bi.sc))) + 
  ggtitle(label = "scanorama") + 
  NoLegend()

pe <- DimPlot(object = bi.v3, reduction = "tsne", group.by = "celltype") +
  ggtitle(label = "Seurat V3") + labs(color = 'Cell Type') + theme(legend.title = element_text(size = 18))

pf <- DimPlot(object = bi.v2, reduction = "tsne", group.by = "celltype") +
  ggtitle(label = "Seurat V2") +
  NoLegend()

pg <- DimPlot(object = bi.mnn, reduction = "tsne", group.by = "celltype") +
  ggtitle(label = "mnnCorrect") + 
  NoLegend()

ph <- DimPlot(object = bi.sc, reduction = "tsne", group.by = "celltype") + 
  ggtitle(label = "scanorama") + 
  NoLegend()

pa <- pa + theme(legend.position = "right")
pe <- pe + theme(legend.position = "right")

l1 <- get_legend(plot = pa)
l2 <- get_legend(plot = pe)

p.final <- plot_grid(
  pa + NoLegend(), pb, pc, pd, l1,
  pe + NoLegend(), pf, pg, ph, l2,
  nrow = 2,
  rel_widths = c(1, 1, 1, 1, 0.3),
  labels = c(LETTERS[1:4], "", LETTERS[5:8]), 
  label_size = 20
)

save_plot(plot = p.final, filename = args[5], base_aspect_ratio = 2/1, base_height = 12.5)
