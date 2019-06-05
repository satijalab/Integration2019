library(Seurat)

starmap.1 <- readRDS("seurat_objects/20180505_BY3_1kgenes.rds")
starmap.2 <- readRDS("seurat_objects/20180410-BY3_1kgenes.rds")
starmap.2 <- RenameCells(starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)
allen <- readRDS("seurat_objects/allen_brain.rds")
dropseq <- readRDS("seurat_objects/dropseq_cortex.rds")

dir.create("seurat_objects/spatial_leaveout/")

# remove HPC from starmap
class_labels.1 <- read.table(
  file = "raw_data/spatial/starmap/visual_1020/20180505_BY3_1kgenes/class_labels.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
class_labels.2 <- read.table(
  file = "raw_data/spatial/starmap/visual_1020/20180410-BY3_1kgenes/class_labels.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

class_labels.1$cellname <- paste0('starmap', rownames(class_labels.1))
class_labels.2$cellname <- paste0('r2_starmap', rownames(class_labels.2))

class_labels.1$ClusterName <- ifelse(is.na(class_labels.1$ClusterName), 'Other', class_labels.1$ClusterName)
class_labels.2$ClusterName <- ifelse(is.na(class_labels.2$ClusterName), 'Other', class_labels.2$ClusterName)

hpc.1 <- class_labels.1[class_labels.1$ClusterName == 'HPC', ]$cellname
hpc.2 <- class_labels.2[class_labels.2$ClusterName == 'HPC', ]$cellname

accept.cells.1 <- setdiff(colnames(starmap.1), hpc.1)
accept.cells.2 <- setdiff(colnames(starmap.2), hpc.2)

starmap.1 <- starmap.1[, accept.cells.1]
starmap.2 <- starmap.2[, accept.cells.2]

starmap.1@misc$spatial <- starmap.1@misc$spatial[starmap.1@misc$spatial$cell %in% accept.cells.1, ]
starmap.2@misc$spatial <- starmap.2@misc$spatial[starmap.2@misc$spatial$cell %in% accept.cells.2, ]

# Integrate the STARmap datasets together
i1 <- FindIntegrationAnchors(
  object.list = c(starmap.1, starmap.2),
  anchor.features = rownames(starmap.1),
  dims = 1:30
)

star <- IntegrateData(
  anchorset = i1,
  dims=1:30,
  k.weight = 100
)

DefaultAssay(star) <- 'integrated'
genes.leaveout <- c('Cux2', 'Zmat4', 'Lamp5', 'Pcsk2', 'Nrsn1', 'Rorb', 'Rab3c', 'Syt6', 'Sox2ot', 'Bsg', 'Sst', 'Pvalb', 'Vip')

run_imputation <- function(ref.obj, query.obj, feature.remove) {
  message(paste0('removing ', feature.remove))
  features <- setdiff(rownames(query.obj), feature.remove)
  anchors <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj,
    features = features,
    dims = 1:30,
    reduction = 'cca'
  )
  refdata <- GetAssayData(
    object = ref.obj,
    assay = 'RNA',
    slot = 'data'
  )
  imputation <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = 'pca'
  )
  query.obj[['seq']] <- imputation
  return(query.obj)
}

for(i in 1:length(genes.leaveout)) {
  imputed.ss2 <- run_imputation(ref.obj = allen, query.obj = star, feature.remove = genes.leaveout[[i]])
  imputed.ds <- run_imputation(ref.obj = dropseq, query.obj = star, feature.remove = genes.leaveout[[i]])
  starmap.1[['ss2']] <- imputed.ss2[, colnames(starmap.1)][['seq']]
  starmap.2[['ss2']] <- imputed.ss2[, colnames(starmap.2)][['seq']]
  starmap.1[['dropseq']] <- imputed.ds[, colnames(starmap.1)][['seq']]
  starmap.2[['dropseq']] <- imputed.ds[, colnames(starmap.2)][['seq']]
  saveRDS(starmap.1, paste0('seurat_objects/spatial_leaveout/starmap_1_', genes.leaveout[[i]], '.rds'))
  saveRDS(starmap.2, paste0('seurat_objects/spatial_leaveout/starmap_2_', genes.leaveout[[i]], '.rds'))
}
