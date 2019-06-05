suppressMessages(library(Seurat))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)

# args <- c(
#   "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_mnnCorrect.rds",
#   "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_scanorama.rds",
#   "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_seuratV2.rds",
#   "~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_seuratV3.rds",
#   "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_mnnCorrect.rds",
#   "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_scanorama.rds",
#   "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_seuratV2.rds",
#   "~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_seuratV3.rds",
#   "~/Projects/muir/figures/supplementary/integration_comparison_tsnes.png"
#   )

# integration methods comparison figures
panc.mnn.correct <- readRDS(file = args[1])
panc.scanorama <- readRDS(file = args[2])
panc.seurat.v2 <- readRDS(file = args[3])
panc.seurat.v3 <- readRDS(file = args[4])
bipolar.mnn.correct <- readRDS(file = args[5])
bipolar.scanorama <- readRDS(file = args[6])
bipolar.seurat.v2 <- readRDS(file = args[7])
bipolar.seurat.v3 <- readRDS(file = args[8])


p1 <- DimPlot(object = panc.seurat.v3, reduction = "tsne", group.by = "replicate") + ggtitle("Seurat V3")
l1 <- get_legend(plot = p1)
p1 <- p1 + NoLegend()
p2 <- DimPlot(object = panc.mnn.correct, reduction = "tsne", group.by = "replicate") + NoLegend() + ggtitle("MNN Correct")
p3 <- DimPlot(object = panc.scanorama, reduction = "umap", group.by = "replicate") + NoLegend() + ggtitle("Scanorama")
p4 <- DimPlot(object = panc.seurat.v2, reduction = "tsne", group.by = "replicate") + NoLegend() + ggtitle("Seurat V2")

celltype.levels <- names(x = table(panc.seurat.v3$celltype))
panc.seurat.v3$celltype <- factor(x = panc.seurat.v3$celltype, levels = celltype.levels)
panc.mnn.correct$celltype <- factor(x = panc.mnn.correct$celltype, levels = celltype.levels)
panc.scanorama$celltype <- factor(x = panc.scanorama$celltype, levels = celltype.levels)
panc.seurat.v2$celltype <- factor(x = panc.seurat.v2$celltype, levels = celltype.levels)

p5 <- DimPlot(object = panc.seurat.v3, reduction = "tsne", group.by = "celltype") + ggtitle("Seurat V3")
l2 <- get_legend(plot = p5)
p5 <- p5 + NoLegend()
p6 <- DimPlot(object = panc.mnn.correct, reduction = "tsne", group.by = "celltype") + NoLegend() + ggtitle("MNN Correct")
p7 <- DimPlot(object = panc.scanorama, reduction = "umap", group.by = "celltype") + NoLegend() + ggtitle("Scanorama")
p8 <- DimPlot(object = panc.seurat.v2, reduction = "tsne", group.by = "celltype") + NoLegend() + ggtitle("Seurat V2")

r1 <- plot_grid(
  p1, p2, p3, p4, l1, 
  nrow = 1,
  labels = LETTERS[1:4],
  rel_widths = c(1,1,1,1,0.4)
)
r2 <- plot_grid(
  p5, p6, p7, p8, l2, 
  nrow = 1,
  labels = LETTERS[5:8],
  rel_widths = c(1,1,1,1,0.4)
)


p9 <- DimPlot(object = bipolar.seurat.v3, reduction = "tsne", group.by = "replicate") + ggtitle("Seurat V3")
l3 <- get_legend(plot = p9)
p9 <- p9 + NoLegend()
p10 <- DimPlot(object = bipolar.mnn.correct, reduction = "tsne", group.by = "replicate") + NoLegend() + ggtitle("MNN Correct")
p11 <- DimPlot(object = bipolar.scanorama, reduction = "umap", group.by = "replicate") + NoLegend() + ggtitle("Scanorama")
p12 <- DimPlot(object = bipolar.seurat.v2, reduction = "tsne", group.by = "replicate") + NoLegend() + ggtitle("Seurat V2")

celltype.levels <- names(x = table(bipolar.seurat.v3$celltype))
bipolar.seurat.v3$celltype <- factor(x = bipolar.seurat.v3$celltype, levels = celltype.levels)
bipolar.mnn.correct$celltype <- factor(x = bipolar.mnn.correct$celltype, levels = celltype.levels)
bipolar.scanorama$celltype <- factor(x = bipolar.scanorama$celltype, levels = celltype.levels)
bipolar.seurat.v2$celltype <- factor(x = bipolar.seurat.v2$celltype, levels = celltype.levels)

p13 <- DimPlot(object = bipolar.seurat.v3, reduction = "tsne", group.by = "celltype") + ggtitle("Seurat V3")
l4 <- get_legend(plot = p13)
p13 <- p13 + NoLegend()
p14 <- DimPlot(object = bipolar.mnn.correct, reduction = "tsne", group.by = "celltype") + NoLegend() + ggtitle("MNN Correct")
p15 <- DimPlot(object = bipolar.scanorama, reduction = "umap", group.by = "celltype") + NoLegend() + ggtitle("Scanorama")
p16 <- DimPlot(object = bipolar.seurat.v2, reduction = "tsne", group.by = "celltype") + NoLegend() + ggtitle("Seurat V2")

r3 <- plot_grid(
  p9, p10, p11, p12, l3, 
  nrow = 1,
  labels = LETTERS[9:12],
  rel_widths = c(1,1,1,1,0.4)
)
r4 <- plot_grid(
  p13, p14, p15, p16, l4, 
  nrow = 1,
  labels = LETTERS[13:16],
  rel_widths = c(1,1,1,1,0.4)
)

pfinal <- plot_grid(r1, r2, r3, r4, nrow = 4)

save_plot(plot = pfinal, filename = args[9], base_aspect_ratio = 1.5, base_height = 17)
 