suppressMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)

#args <- c("~/Projects/muir/seurat_objects/pbmc3k.rds", "~/Projects/muir/seurat_objects/pbmc33k.rds")

pbmc <- readRDS(file = args[1])
pbmc33k <- readRDS(file = args[2])

predictions <- ProjectCells(
  query = pbmc,
  reference = pbmc33k,
  dims = 1:25,
  label.name = "cluster",
  normalize.with.lm = FALSE,
  eps = 0
)

saveRDS(object = predictions, file = args[3])
