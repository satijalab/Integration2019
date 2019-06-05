
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[4])
# args <- c("seurat_objects/citeseq.rds", "seurat_objects/integrated_bone_marrow.rds")

citeseq <- readRDS(file = args[1])
hca <- readRDS(file = args[2])

DefaultAssay(object = citeseq) <- "RNA"
DefaultAssay(object = hca) <- "integrated"
features <- intersect(x = rownames(x = citeseq), y = rownames(x = hca))

citeseq <- RunPCA(citeseq, npcs = 100)

transfor.anchors <- FindTransferAnchors(
  reference = citeseq,
  query = hca,
  features = features,
  dims = 1:50,
  npcs = 50,
  verbose = TRUE
)

refdata <- GetAssayData(
  object = citeseq,
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

hca[['ADT']] <- imputed.data

# save object
saveRDS(object = hca, file = args[3])
