suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[4])

# args <- c("~/Projects/muir/seurat_objects/integrated_bipolar.rds")

bipolar <- readRDS(file = args[1])
Idents(object = bipolar) <- 'celltype'
# Rename celltypes to make them easier to work with
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Amacrine cells")) <- "Amacrine_cells"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "BC8/BC9")) <- "BC8_BC9"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Cone photoreceptors")) <- "Cone_photoreceptors"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Müller glia")) <- "Muller_glia"
Idents(object = bipolar, cells = WhichCells(object = bipolar, idents = "Rod  photoreceptors")) <- "Rod_photoreceptors"

# Remove Unknown and Doublet clusters
bipolar <- SubsetData(object = bipolar, ident.remove = c("Unknown", "Doublets"))
Idents(object = bipolar) <- factor(Idents(object = bipolar))
bipolar[['celltype']] <- Idents(object = bipolar)

cluster.to.drop <- c("RBC", "Muller_glia", "BC5A", "BC7", "BC6", "BC1A")
datasets <- 1:6

# drop 1 cluster (different) in each sample
for(i in 1:6){
  cells.remove <- intersect(
    x = WhichCells(object = bipolar, idents = cluster.to.drop[i]), 
    y = WhichCells(object = bipolar, expression = replicate == datasets[i])
  )
  bipolar <- SubsetData(object = bipolar, cells = setdiff(x = colnames(x = bipolar), y = cells.remove))
}

removal.table <- data.frame(cluster = cluster.to.drop, dataset = datasets)
dir.create(path = paste0(getwd(), "/tables"), showWarnings = FALSE)
write.csv(x = removal.table, file = args[2], quote = FALSE)

# select integration features
DefaultAssay(object = bipolar) <- "RNA"
all.bipolar <- SplitObject(object = bipolar, split.by = "replicate")
for(i in 1:length(x = all.bipolar)) {
  all.bipolar[[i]] <- FindVariableFeatures(object = all.bipolar[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
genes.use <- SelectIntegrationFeatures(object.list = all.bipolar, nfeatures =  2000)

bipolar@misc$integration.features <- genes.use

saveRDS(object = bipolar, file = args[3])
