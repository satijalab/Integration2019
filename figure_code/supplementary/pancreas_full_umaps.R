suppressMessages(library(Seurat))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)

args <- c("~/Projects/muir/seurat_objects/all_pancreas_no_integration.rds",
          "~/Projects/muir/seurat_objects/integrated_pancreas.rds",
          "~/Projects/muir/figures/supplementary/pancreas_full_umaps.pdf")

pi.noint <- readRDS(file = args[1])
pi <- readRDS(file = args[2])

p1 <- DimPlot(object = pi.noint, reduction = "umap", group.by = "replicate") + 
  theme(legend.position = "top", legend.justification = "center") + 
  ggtitle("Unintegrated Datasets")
p2 <- DimPlot(object = pi, reduction = "umap", group.by = "replicate") + 
  theme(legend.position = "top", legend.justification = "center") + 
  ggtitle("Integrated Datasets")
p3 <- DimPlot(object = pi, reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()

pfinal <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "lt", labels = LETTERS[1:3])

save_plot(plot = pfinal, filename = args[3], base_aspect_ratio = 2/1, base_height = 10)
