suppressMessages(library(Seurat))
suppressMessages(library(cowplot))
suppressMessages(library(future))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("seurat_objects/integrated_mca.rds", "figures/figureS3.pdf")

mca <- readRDS(file = args[1])

plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 10^12)

mca$overlap <- 'SMART-Seq + 10X'
mca@meta.data[mca$tissue %in% names(x = which(x = apply(X = table(mca$tech, mca$tissue), MARGIN = 2, FUN = min) == 0)), 'overlap'] <- 'SMART-Seq only'

p1 <- DimPlot(object = mca, reduction = 'tsne', group.by = 'tech', cells = sample(x = Cells(x = mca))) + 
  theme(legend.position="top", legend.justification = "center", legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 5)))

p2 <- DimPlot(object = mca, reduction = 'tsne', group.by = 'tissue') + 
  theme(legend.position = "top", legend.justification = "center", legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 5)))

p3 <- DimPlot(object = mca, reduction = 'tsne', label = FALSE, group.by = 'overlap') + 
  theme(legend.position="top", legend.justification = "center", legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

mix_mca <- MixingMetric(object = mca, dims = 1:100, grouping.var = 'tech', eps = 1, verbose = TRUE, k = 20)
mca$mix <- max(mix_mca) - mix_mca
mca_local <- LocalStruct(object = mca, grouping.var = 'tech', reduced.dims = 1:100, orig.dims = 1:100, reduction = "pca")
mca$local <- unlist(mca_local)

p4a <- VlnPlot(object = mca, features = "mix", group.by = "overlap", pt.size = 0) +
  ylab("Mixing Metric") + ggtitle("Mixing Metric") + theme(legend.position = 'none')

p4b <- VlnPlot(object = mca, features = "local", group.by = "overlap", pt.size = 0) +
  ylab("Local Structure Metric") + ggtitle("Local Structure Metric") + theme(legend.position = 'none')

p4 <- plot_grid(p4a, p4b)

# load("~/xfer/mca_markers.Rda")
markers1 <- c("Msln", 'Upk3b', 'Lrrn4', 'Pkhd1l1', 'Wt1', 'Rspo1', 'Chst4', 'Tmem151a', 'Lgals2', 'Upk1b', 'Cst9', 'Cldn15', 'Ildr2', 'Tm4sf5', 'Calcb')
markers2 <- c('Igj', 'Prg2', 'Derl3', 'Tnfrsf17', 'Cacna1s', 'Eaf2', 'Chst11', 'Lax1', 'Ly6c2', 'Fkbp11', 'BC031353', 'St6gal1', 'Txndc11')

DefaultAssay(object = mca) <- 'RNA'
mcates <- subset(x = mca, downsample = 500)
Idents(object = mcates) <- 'tissue'
cell_label <- intersect(x = WhichCells(mca, idents = 121), y = Cells(mcates))

Idents(mcates, cells = cell_label) <- paste("Mesothelial_", mca$tissue[cell_label], sep = "")

p5 <- DotPlot(object = mcates, features = markers1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cell_label <- intersect(x = WhichCells(object = mca, idents = 114), y = Cells(x = mcates))

Idents(object = mcates) <- 'tissue'
Idents(object = mcates, cells = cell_label) <- paste("pDC_", mca$tissue[cell_label], sep = "")

p6 <- DotPlot(object = mcates, features = markers2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

mcaFig <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
save_plot(plot = mcaFig, filename = args[2], base_height = 20, base_width = 30)

DimPlot(mca, label = TRUE, reduction = 'tsne') +
  ggsave("figures/mca_labelled.pdf")
