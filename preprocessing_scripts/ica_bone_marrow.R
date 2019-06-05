library(Matrix)


args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

input.file <- "raw_data/ica/ica_bone_marrow.h5"

bm <- Read10X_h5(input.file)
rownames(bm) <- make.unique(rownames(bm))
ica_bone <- CreateSeuratObject(counts = bm, project = 'ICA',
                               min.cells = 100, min.features = 500,
                               names.field = 1, names.delim = '_')
ica_bone <- NormalizeData(object = ica_bone)
ica_bone <- FindVariableFeatures(object = ica_bone, nfeatures = 2000)
ica_bone <- ScaleData(object = ica_bone)
saveRDS(object = ica_bone, file = args[2])

# save individual matrices for scrublet
bm.list <- SplitObject(ica_bone)
for (i in 1:length(bm.list)) {
  counts <- GetAssayData(bm.list[[i]], assay = 'RNA', slot = 'counts')
  writeMM(obj = counts,
          fil = paste0(getwd(), '/raw_data/ica/', names(bm.list)[[i]], '.mtx'))
}
