library(Seurat)
library(ggplot2)
library(future)
library(future.apply)
set.seed(1234)

plan(strategy = 'multicore', workers = 5)
options(future.globals.maxSize= 10^12)

dir.create("seurat_objects/citeseq_downsampling/")

args <- commandArgs(trailingOnly = TRUE)
citeseq <- readRDS(file = args[1])

DefaultAssay(object = citeseq) <- "RNA"
citeseq <- FindVariableFeatures(object = citeseq, nfeatures = 2000)

reference <- sample(x = colnames(x = citeseq), size = ncol(x = citeseq) / 2, replace = FALSE)
query <- setdiff(x = colnames(x = citeseq), y = reference)
ref.obj <- citeseq[, reference]
query.obj <- citeseq[, query]

ranking <- abs(x = Loadings(ref.obj[['pca']]) * ref.obj[['pca']]@stdev)
ranking <- sort(x = rowSums(x = ranking), decreasing = TRUE)
query.obj[["ADT"]] <- NULL

dir.create(file.path(paste0(getwd(), "/seurat_objects/citeseq_downsampling/")), showWarnings = FALSE)

future_sapply(seq(10, 1000, 10), function(x) {
  features <- names(x = head(x = ranking, x))
  ref.obj <- ScaleData(object = ref.obj, features = features)
  query.obj <- ScaleData(object = query.obj, features = features)
  dims = min(60, x)
  if (x < 100) {
    approx = FALSE
  } else {
    approx = TRUE
  }
  transfor.anchors <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj,
    features = features,
    dims = 1:dims,
    npcs = dims,
    approx.pca = approx,
    verbose = TRUE
  )
  refdata <- GetAssayData(
    object = ref.obj,
    assay = 'ADT',
    slot = 'data'
  )
  imputed.data <- TransferData(
    anchorset = transfor.anchors,
    refdata = refdata,
    l2.norm = FALSE,
    dims = 1:dims,
    k.weight = 50,
    slot = 'counts'
  )
  save.name <- paste0("seurat_objects/citeseq_downsampling/", as.character(x = x), "_feature_downsampling.rds")
  saveRDS(object = GetAssayData(imputed.data, slot = 'data'), file = save.name)
  rm(imputed.data)
  gc()
})
