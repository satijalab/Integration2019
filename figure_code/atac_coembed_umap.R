library(Seurat)
library(ggplot2)

atac <- readRDS("seurat_objects/pfc_atac_coembed.rds")
atac$combined_ident <- ifelse(is.na(atac$subclass), atac$predicted.id, atac$subclass)

DimPlot(
  object = atac,
  reduction = 'umap',
  group.by = 'dataset'
) + theme(legend.position = 'none') +
  ggsave("figures/atac/pfc/coembed_atac_dataset.png", height = 6, width = 6, dpi = 800)

# reorder celltypes so the inhibitory neurons get very different colors
atac$combined_ident <- factor(atac$combined_ident, levels = c("Astro", "CR", "Lamp5", "L2/3 IT", "L6 IT", "L5 IT",
                                                              "L5 PT", "L6 CT", "Sst", "L6b", "Endo", "Macrophage",
                                                              "Meis2", "NP", "Oligo", "Vip", "Serpinf1", "Pvalb",
                                                              "SMC", "Sncg", "L4", "Peri", "VLMC"))

DimPlot(
  object = atac,
  reduction = 'umap',
  group.by = 'combined_ident',
  label = TRUE,
  repel = TRUE
) + theme(legend.position = 'none') +
  ggsave("figures/atac/pfc/coembed_atac_celltype.png", height = 6, width = 6, dpi = 800)

filtered <- subset(atac, subset = dataset == 'scATAC')

DimPlot(
  object = filtered,
  reduction = 'umap',
  group.by = 'combined_ident',
  label = TRUE,
  repel = TRUE
) + theme(legend.position = 'none') +
  ggsave("figures/atac/pfc/coembed_atac_celltype_only_atac.png", height = 6, width = 6, dpi = 800)
