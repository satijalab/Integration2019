suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/software/seurat-private", "~/Projects/muir/seurat_objects/celseq.rds",
#           "~/Projects/muir/seurat_objects/celseq2.rds", "~/Projects/muir/seurat_objects/smartseq2.rds",
#           "~/Projects/muir/seurat_objects/fluidigmc1.rds", "~/Projects/muir/seurat_objects/indrop.rds")

devtools::load_all(args[1])

celseq <- readRDS(file = args[2])
celseq2 <- readRDS(file = args[3])
smartseq2 <- readRDS(file = args[4])
fluidigmc1 <- readRDS(file = args[5])
indrop <- readRDS(file = args[6])
all.pancreas <- c(list(celseq, celseq2, smartseq2, fluidigmc1), indrop)

pancreas.anchors <- FindIntegrationAnchors(
  object.list = all.pancreas, 
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

pancreas.integrated <- IntegrateData(
  anchorset = pancreas.anchors, 
  new.assay.name = "integrated", 
  dims = 1:30, 
  k.weight = 100, 
  sd.weight = 1, 
  eps = 0, 
  verbose = TRUE
)

# take through clustering and annotate
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(
  object = pancreas.integrated, 
  features = VariableFeatures(object = pancreas.integrated)
)
pancreas.integrated <- RunPCA(
  object = pancreas.integrated,
  features = VariableFeatures(object = pancreas.integrated),
  verbose = FALSE, 
  npcs = 30
)
pancreas.integrated <- RunTSNE(
  object = pancreas.integrated,
  reduction = "pca",
  dims = 1:30, 
  tsne.method = "FIt-SNE"
)
pancreas.integrated <- RunUMAP(
  object = pancreas.integrated,
  reduction = "pca",
  dims = 1:30,
  nneighbors = 5
)
pancreas.integrated <- FindNeighbors(
  object = pancreas.integrated, 
  reduction = "pca", 
  dims = 1:30, 
  k.param = 20
)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 1.5)


# annotation
old.ids <- as.vector(x = unique(x = Idents(object = pancreas.integrated)))
new.ids <- unlist(x = sapply(
  X = old.ids,
  FUN = function(x){
    cluster.cells <- WhichCells(object = pancreas.integrated, idents = x)
    names(x = sort(x = table(pancreas.integrated$assigned_cluster[cluster.cells]), decreasing = TRUE)[1])
  }
))
new.ids <- new.ids[as.character(x = Idents(object = pancreas.integrated))]
new.ids[is.na(x = new.ids)] <- "Unknown"
pancreas.integrated[["celltype"]] <- new.ids
Idents(pancreas.integrated) <- "celltype"

pancreas.integrated$replicate[pancreas.integrated$replicate == "human1"] <- "indrop1"
pancreas.integrated$replicate[pancreas.integrated$replicate == "human2"] <- "indrop2"
pancreas.integrated$replicate[pancreas.integrated$replicate == "human3"] <- "indrop3"
pancreas.integrated$replicate[pancreas.integrated$replicate == "human4"] <- "indrop4"

saveRDS(object = pancreas.integrated, file = args[7])

pi.celltypes <- data.frame(cell = rownames(pancreas.integrated[[]]), celltype = pancreas.integrated$celltype, dataset = pancreas.integrated$replicate)
cat("Supplementary Table 2: Cell type annotations for integrated pancreas analysis\n\n", file = args[8])
write.table(x = pi.celltypes, file = args[8], quote = FALSE, append = TRUE, sep = ",", row.names = FALSE)

# For viz/annotation
# DimPlot(pancreas.integrated, reduction = "tsne", label = TRUE)
# DimPlot(pancreas.integrated, reduction.use = "tsne", group.by = "replicate")

# FeaturePlot(object = pancreas.integrated,  features = c("GCG", "INS", "PPY", "SST"))
# FeaturePlot(object = pancreas.integrated,  features = c("PRSS1", "KRT19", "VWF", "PDGFRA"))
#   - alpha (GCG)
#   - beta (INS)
#   - ductal (KRT19)
#   - acinar (PRSS1)
#   - delta (SST)
#   - stellate (PDGFRA/RGS5)
#   - gamma (PPY)
#   - endothelial (VWF)
#   - immune (CD53, SDS) 
