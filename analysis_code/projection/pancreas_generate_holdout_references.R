args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/celseq.rds", "~/Projects/muir/seurat_objects/celseq2.rds",
#           "~/Projects/muir/seurat_objects/smartseq2.rds", "~/Projects/muir/seurat_objects/fluidigmc1.rds",
#           "~/Projects/muir/seurat_objects/inDrop.rds", "~/Projects/muir/seurat_objects/integrated_pancreas.rds",
#           "alpha", "celseq")
devtools::load_all(args[9])

celseq <- readRDS(file = args[1])
celseq2 <- readRDS(file = args[2])
smartseq2 <- readRDS(file = args[3])
fluidigmc1 <- readRDS(file = args[4])
indrop <- readRDS(file = args[5])
pancreas.integrated <- readRDS(file = args[6])
cell.type.removed <- args[7]
query.name <- args[8]
dataset.names <- c("celseq", "celseq2", "smartseq2", "fluidigmc1", "indrop1", "indrop2", "indrop3", "indrop4")
query.idx <- which(dataset.names == query.name)

ref.master <- c(celseq, celseq2, smartseq2, fluidigmc1, indrop)

for(i in 1:length(ref.master)){
  ref.master[[i]][["celltype"]] <- Idents(object = pancreas.integrated)
  Idents(ref.master[[i]]) <- "celltype"
}

ref.savename <- paste0("analysis_data/pancreas/reference-without-", query.name, "-", cell.type.removed, ".rds")
query.savename <- paste0("analysis_data/pancreas/query-without-", query.name, "-", cell.type.removed, ".rds")

# leave one dataset out of reference, remove given cell type from all but query
ref <- ref.master
query <- ref.master[[query.idx]]
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
