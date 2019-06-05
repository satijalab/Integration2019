suppressMessages(library(Seurat))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)
#args <- c("~/Projects/muir/seurat_objects/pfc_atac_coembed.rds", "~/Projects/muir/figures/pfc_atac_coembed_umaps.pdf" )

coembed <- readRDS(file = args[1])

# sample 5000 allen cells so that they do not dominate the plot

cells.display <- c(
  sample(x = WhichCells(object = coembed, expression = dataset == "Allen"), size = 5000),
  WhichCells(object = coembed, expression = dataset == "scATAC")
)
p1 <- DimPlot(object = coembed, cells = cells.display, group.by = "dataset") + 
  theme(legend.position = 'top', legend.justification = "center")
p2 <- DimPlot(object = coembed, cells = cells.display, group.by = 'subclass', label = TRUE) + 
  theme(legend.position = 'top', legend.justification = "center")

pfinal <- plot_grid(p1, p2, align = "h")

save_plot(filename = args[2], plot = pfinal, base_height = 9, base_width = 16)
