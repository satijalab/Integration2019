args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[4])

pancreas.integrated <- readRDS(file = args[1])
Idents(object = pancreas.integrated) <- "celltype"

# merge the stellates
Idents(object = pancreas.integrated, cells = WhichCells(object = pancreas.integrated, idents = "activated_stellate")) <- "stellate"
Idents(object = pancreas.integrated, cells = WhichCells(object = pancreas.integrated, idents = "quiescent_stellate")) <- "stellate"

cluster.to.drop <- c("alpha", "beta", "endothelial", "ductal", "delta", "acinar", "stellate", "gamma")
datasets <- c("celseq","celseq2","smartseq2","fluidigmc1","indrop1","indrop2","indrop3","indrop4")
pancreas.integrated[["celltype"]] <- Idents(object = pancreas.integrated)
# drop 1 cluster (different) in each sample
for(i in 1:8){
  cells.remove <- intersect(
    x = WhichCells(object = pancreas.integrated, idents = cluster.to.drop[i]), 
    y = WhichCells(object = pancreas.integrated, expression = replicate == datasets[i])
  )
  pancreas.integrated <- SubsetData(object = pancreas.integrated, cells = setdiff(x = colnames(x = pancreas.integrated), y = cells.remove))
}

removal.table <- data.frame(cluster = cluster.to.drop, dataset = datasets)
dir.create(path = paste0(getwd(), "/tables"), showWarnings = FALSE)
write.csv(x = removal.table, file = args[2], quote = FALSE)

pancreas.integrated.save <- pancreas.integrated

# select integration features
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated[["integrated"]] <- NULL

all.pancreas <- SplitObject(object = pancreas.integrated, split.by = "replicate")
for(i in 1:length(x = all.pancreas)) {
  all.pancreas[[i]] <- FindVariableFeatures(object = all.pancreas[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
genes.use <- SelectIntegrationFeatures(object.list = all.pancreas, nfeatures =  2000)

pancreas.integrated.save@misc$integration.features <- genes.use
saveRDS(object = pancreas.integrated.save, file = args[3])
