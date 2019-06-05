suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[2])

pancreas.integrated <- readRDS(file = args[1])
DefaultAssay(object = pancreas.integrated) <- "RNA"
all.pancreas <- SplitObject(object = pancreas.integrated, split.by = "replicate")
dir.create(path = paste0(getwd(), "/analysis_data/pancreas/"), showWarnings = FALSE)
conf.savenames <- c()

for(i in 1:length(x = all.pancreas)) {
  dataset <- unique(x = all.pancreas[[i]]$replicate)
  message("Writing ", dataset)
  savename <- paste0(getwd(), "/analysis_data/pancreas/", dataset ,"_holdout.txt")
  write.table(
   x = as.matrix(x = GetAssayData(object = all.pancreas[[i]], assay = "RNA", slot = "counts")),
   file = savename,
   sep = "\t",
   col.names = NA,
   row.names = TRUE
  )
  conf.savenames <- c(conf.savenames, savename)
}
write.table(
  x = conf.savenames, 
  file = paste0(getwd(), "/analysis_data/pancreas/pancreas_conf.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
conf2.savenames <- unname(sapply(X = conf.savenames, FUN = function(x) Seurat:::ExtractField(string = x, delim = "\\.", field = 1)))
write.table(
  x = conf2.savenames, 
  file = paste0(getwd(), "/analysis_data/pancreas/pancreas_conf2.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
