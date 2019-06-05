args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[4])

#args <- c("~/Projects/muir/seurat_objects/bipolar.rds", "Rod_photoreceptors", "2")

bipolar <- readRDS(file = args[1])
Idents(object = bipolar) <- 'celltype'
# Rename celltypes to make them easier to work with
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Amacrine cells")) <- "Amacrine_cells"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "BC8/BC9")) <- "BC8_BC9"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Cone photoreceptors")) <- "Cone_photoreceptors"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Müller glia")) <- "Muller_glia"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Rod  photoreceptors")) <- "Rod_photoreceptors"
bipolar[['celltype']] <- Idents(object = bipolar)

bipolar <- SubsetData(object = bipolar, ident.remove = c("Unknown", "Doublets"))

all.bipolar <- SplitObject(object = bipolar, split.by = "replicate")
cell.type.removed <- args[2]
query.name <- args[3]
dataset.names <- as.character(1:6)
query.idx <- which(dataset.names == query.name)

ref.savename <- paste0("analysis_data/bipolar/reference-without-", query.name, "-", cell.type.removed, ".rds")
query.savename <- paste0("analysis_data/bipolar/query-without-", query.name, "-", cell.type.removed, ".rds")

# leave one dataset out of reference, remove given cell type from all but query
ref <- all.bipolar
query <- all.bipolar[[query.idx]]

if (table(query[["celltype"]])[cell.type.removed] < 10  || is.na(table(query[["celltype"]])[cell.type.removed] )) {
  # empty filler
  filler <- numeric()
  saveRDS(object = filler, file = ref.savename)
  saveRDS(object = filler, file = query.savename)
} else {
  ref[[query.idx]] <- NULL
  for(j in 1:length(ref)) {
    if(cell.type.removed %in% unique(ref[[j]][["celltype"]])$celltype) {
      ref[[j]] <- SubsetData(object = ref[[j]], ident.remove = cell.type.removed)
    }
  }
  # gene selection
  for(j in 1:length(x = ref)){
    ref[[j]] <- FindVariableFeatures(
      object = ref[[j]], 
      selection.method = "vst", 
      nfeatures = 2000, 
      verbose = FALSE)
  }
  
  ref.anchors <- FindIntegrationAnchors(
    object.list = ref, 
    anchor.features = 2000, 
    scale = TRUE, 
    l2.norm = TRUE, 
    dims = 1:30, 
    k.anchor = 5, 
    k.filter = 200, 
    k.score = 30, 
    max.features = 200, 
    eps = 0, 
    verbose = TRUE
  )
  
  ref.integrated <- IntegrateData(
    anchorset = ref.anchors, 
    new.assay.name = "integrated", 
    dims = 1:30, 
    k.weight = 100, 
    sd.weight = 1, 
    eps = 0, 
    verbose = TRUE
  )
  
  DefaultAssay(object = ref.integrated) <- "integrated"
  ref.integrated <- ScaleData(object = ref.integrated, verbose = FALSE)
  ref.integrated <- RunPCA(object = ref.integrated,
                           features.use = VariableFeatures(ref.integrated),
                           verbose = FALSE, npcs = 50)
  
  # downsample query to 100 cells (max) per ident with removed pop being 20%
  Idents(object = query) <- "celltype"
  query1 <- SubsetData(object = query, max.cells.per.ident = 100, ident.remove = cell.type.removed)
  newCells <- round(ncol(query1)/4)
  query2 <- SubsetData(query, ident.use = cell.type.removed, max.cells.per.ident = newCells)
  query <- SubsetData(query, cells = c(colnames(query1), colnames(query2)))
  saveRDS(object = ref.integrated, file = ref.savename)
  saveRDS(object = query, file = query.savename)
}
