library(Seurat)
library(ggplot2)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("seurat_objects/allen_brain_upper.rds", "seurat_objects/atac_pfc.rds", "seurat_objects/atac_pfc_coembed.rds", "seurat_objects/atac_pfc_predictions.rds")
upallen <- readRDS(file = args[1])
pfc <- readRDS(file = args[2])
pfc <- RunUMAP(pfc, reduction = 'lsi', dims = 1:30)
DefaultAssay(pfc) <- 'RNA'
upallen <- FindVariableFeatures(object = upallen, nfeatures = 5000)
genes.use <- intersect(x = VariableFeatures(object = upallen)[1:5000], y = rownames(x = pfc))
upallen$dataset <- 'Allen'
pfc$dataset <- 'scATAC'

dir.create(path = "figures/atac/pfc/")
dir.create(path = "analysis_data/atac/clusters/")

pfc.original <- pfc

# find transfer anchors and predict class labels
transfer.anchors <- FindTransferAnchors(
  reference = upallen,
  query = pfc,
  features = genes.use,
  reduction = 'cca'
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = upallen$subclass,
  weight.reduction = pfc[['lsi']]
)

pfc.stash <- pfc

pfc <- AddMetaData(object = pfc, metadata = predictions)
Idents(object = pfc) <- 'predicted.id'
pfc$subclass <- pfc$predicted.id

# impute gene expression for the ATAC cells. 
# Note that this is only necessary for visualization of a co-embedding
refdata <- GetAssayData(
  object = upallen[["RNA"]],
  slot = "data"
)[VariableFeatures(upallen), ]

imputed.expression <- TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = 'cca',
  l2.norm = FALSE
)

pfc[["imputed"]] <- imputed.expression

pfc.embed <- pfc

coembed.assay <- merge(
  x = subset(x = upallen[["RNA"]], features = genes.use), 
  y = subset(x = pfc.embed[["imputed"]], features = genes.use)
)
coembed <- merge(
  x = upallen[genes.use, ],
  y = pfc.embed[genes.use, ]
)

# create new assay due to issue where cell order gets messed up -- need to fix later
new.assay <- new(
  Class = "Assay",
  counts = coembed.assay@counts[, colnames(coembed)],
  data = coembed.assay@data[, colnames(coembed)],
  scale.data = matrix()
)

coembed[['combined']] <- new.assay
DefaultAssay(object = coembed) <- 'combined'

coembed <- ScaleData(object = coembed, features = genes.use)
coembed <- RunPCA(object = coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(object = coembed, dims = 1:30)

# create new identities, using original atac classifications and the cluster 5 identities
Idents(pfc) <- "peaks_snn_res.1"
c7.cells <- Idents(pfc)
c7.cells <- c7.cells[c7.cells == 8]
Idents(pfc) <- 'cell_label'
Idents(pfc, cells = names(c7.cells)) <- "Unknown"


VlnPlot(pfc, c("GNG5", "SPATA1", "PRKG2", "FAM19A1", "THBD", "OLFR1423"), ncol = 2, pt.size = 0.1) +
  ggsave("figures/atac/pfc/vln_plot_unknown_cells.png", height = 10, width = 6)

saveRDS(object = coembed, file = args[3])
saveRDS(object = pfc, file = args[4])
write.table(predictions, 'analysis_data/atac_celltype_prediction.tsv', sep = '\t', quote = FALSE)

DimPlot(coembed, cells.highlight = names(c7.cells), pt.size = 0.1) + NoLegend() +
  ggsave("figures/atac/pfc/nonoverlap_population_dimplot.png", height = 6, width = 6)


confident <- subset(pfc, subset = prediction.score.max > 0.5)
Idents(confident) <- 'predicted.id'
vip <- WhichCells(confident, ident = 'Vip')
pv <- WhichCells(confident, ident = 'Pvalb')
sst <- WhichCells(confident, ident = 'Sst')
lamp5 <- WhichCells(confident, ident = 'Lamp5')
l4 <- WhichCells(confident, ident = 'L4')

writeLines(vip, "analysis_data/atac/clusters/vip.txt")
writeLines(pv, "analysis_data/atac/clusters/pvalb.txt")
writeLines(sst, "analysis_data/atac/clusters/sst.txt")
writeLines(lamp5, "analysis_data/atac/clusters/lamp5.txt")
writeLines(l4, "analysis_data/atac/clusters/l4.txt")
