library(Seurat)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)

citeseq <- readRDS(file = args[1])
DefaultAssay(object = citeseq) <- "RNA"
citeseq <- FindVariableFeatures(object = citeseq, nfeatures = 2000)
features <- VariableFeatures(object = citeseq)

# split into query/reference
reference <- sample(x = colnames(x = citeseq), size = ncol(x = citeseq) / 2, replace = FALSE)
query <- setdiff(x = colnames(x = citeseq), y = reference)
ref.obj <- citeseq[, reference]
query.obj <- citeseq[, query]
ref.obj <- ScaleData(object = ref.obj, features = features)
query.obj <- ScaleData(object = query.obj, features = features)

# remove measured ADTs from query
query.adt <- query.obj[['ADT']]
query.obj[["ADT"]] <- NULL

transfor.anchors <- FindTransferAnchors(
  reference = ref.obj,
  query = query.obj,
  features = features,
  dims = 1:50,
  npcs = 50,
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
  dims = 1:50,
  k.weight = 50,
  slot = 'counts'
)

query.obj[['imputed']] <- imputed.data

# add back original ADT
query.obj[['ADT']] <- query.adt

# save object
saveRDS(object = query.obj, file = args[2])
