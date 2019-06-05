# put figure 2 panels together

library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[8])

# args <- c("~/Projects/muir/figures/integration_dimplots.rda",
#           "~/Projects/muir/figures/anchor_scores.rds",
#           "~/Projects/muir/figures/anchor_barplot.rds",
#           "~/Projects/muir/figures/silhouette.rds",
#           "~/Projects/muir/figures/mixing_ls_metrics.rds",
#           "~/Projects/muir/figures/supplementary/pancreas_full_umaps.pdf",
#           "~/Projects/muir/figures/figure2.pdf")

load(file = args[1])
anchors <- readRDS(file = args[2])
anchors.barplot <- readRDS(file = args[3])
sil <- readRDS(file = args[4])
mixing.ls <- readRDS(file = args[5])

pa <- pa + theme(legend.position = "right") + 
  labs(color = 'Dataset')
pe <- pe + theme(legend.position = "right") + 
  labs(color = 'Cell Type')

l1 <- get_legend(plot = pa)
l2 <- get_legend(plot = pe)

pa <- AugmentPlot(plot = pa, dpi = 200) + ggtitle("Seurat v3") + theme(plot.title = element_text(size = 18))
pb <- AugmentPlot(plot = pb, dpi = 200) + ggtitle("Seurat v2") + theme(plot.title = element_text(size = 18))
pc <- AugmentPlot(plot = pc, dpi = 200) + ggtitle("mnnCorrect") + theme(plot.title = element_text(size = 18))
pd <- AugmentPlot(plot = pd, dpi = 200) + ggtitle("scanorama") + theme(plot.title = element_text(size = 18))

pe <- AugmentPlot(plot = pe, dpi = 200) + ggtitle("Seurat v3") + theme(plot.title = element_text(size = 18))
pf <- AugmentPlot(plot = pf, dpi = 200) + ggtitle("Seurat v2") + theme(plot.title = element_text(size = 18))
pg <- AugmentPlot(plot = pg, dpi = 200) + ggtitle("mnnCorrect") + theme(plot.title = element_text(size = 18))
ph <- AugmentPlot(plot = ph, dpi = 200) + ggtitle("scanorama") + theme(plot.title = element_text(size = 18))


ps <- plot_grid(
  pa + NoLegend(), pb, pc, pd, l1,
  pe + NoLegend(), pf, pg, ph, l2,
  nrow = 2,
  rel_widths = c(1, 1, 1, 1, 0.3)
)

save_plot(plot = ps, filename = args[6], base_aspect_ratio = 2/1, base_height = 12.5)

p <- plot_grid(
  pa + NoLegend(), pb, pc, pd, l1,
  pe + NoLegend(), pf, pg, ph, l2,
  nrow = 2,
  rel_widths = c(1, 1, 1, 1, 0.3)
)

p2 <- plot_grid(
  anchors, anchors.barplot, sil, mixing.ls, 
  ncol = 4,
  rel_heights = c(0.75, 0.75, 1, 1),
  rel_widths = c(0.75, 0.75, 1, 1)
)

p.final <- plot_grid(
  p, p2,
  nrow = 2, 
  rel_heights = c(1, 0.45)
)

save_plot(plot = p.final, filename = args[7], base_width = 21, base_height = 15)
