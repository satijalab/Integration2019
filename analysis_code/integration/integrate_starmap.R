args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[4])

# args <- c("~/Projects/muir/seurat_objects/20180505_BY3_1kgenes.rds", "~/Projects/muir/seurat_objects/20180410-BY3_1kgenes.rds", 
#           "~/Projects/muir/seurat_objects/integrated_starmap.rds", "~/Projects/seurat-private/")

starmap.1 <- readRDS(file = args[1])
starmap.2 <- readRDS(file = args[2])

starmap.2 <- RenameCells(starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)

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

accept.cells.1 <- setdiff(x = colnames(starmap.1), y = hpc.1)
accept.cells.2 <- setdiff(x = colnames(starmap.2), y = hpc.2)

starmap.1 <- starmap.1[, accept.cells.1]
starmap.2 <- starmap.2[, accept.cells.2]

starmap.1@misc$spatial <- starmap.1@misc$spatial[starmap.1@misc$spatial$cell %in% accept.cells.1, ]
starmap.2@misc$spatial <- starmap.2@misc$spatial[starmap.2@misc$spatial$cell %in% accept.cells.2, ]

# integrate the STARmap datasets together
i1 <- FindIntegrationAnchors(
  object.list = c(starmap.1, starmap.2),
  anchor.features = rownames(x = starmap.1),
  dims = 1:30
)

star <- IntegrateData(anchorset = i1)
DefaultAssay(object = star) <- 'integrated'
star <- ScaleData(object = star)
star <- RunPCA(object = star, npcs = 30, verbose = FALSE)
star <- RunUMAP(object = star, dims = 1:30)

saveRDS(object = star, file = args[3])

