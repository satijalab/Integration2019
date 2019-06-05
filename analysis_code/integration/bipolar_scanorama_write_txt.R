library(Matrix, quietly = TRUE)
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[2])

bipolar.integrated <- readRDS(file = args[1])
DefaultAssay(object = bipolar.integrated) <- "RNA"
bipolar.integrated[["integrated"]] <- NULL
all.bipolar <- SplitObject(object = bipolar.integrated, split.by = "replicate")
dir.create(path = paste0(getwd(), "/analysis_data/bipolar/"), showWarnings = FALSE)
conf.savenames <- c()

for(i in 1:length(x = all.bipolar)) {
  dataset <- unique(x = all.bipolar[[i]]$replicate)
  message("Writing ", dataset)
  savename <- paste0(getwd(), "/analysis_data/bipolar/", dataset ,"_holdout.txt")
  write.table(
    x = as.matrix(x = GetAssayData(object = all.bipolar[[i]], assay = "RNA", slot = "counts")),
    file = savename,
    sep = "\t",
    col.names = NA,
    row.names = TRUE
  )
  conf.savenames <- c(conf.savenames, savename)
}
write.table(
  x = conf.savenames, 
  file = paste0(getwd(), "/analysis_data/bipolar/bipolar_conf.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
conf2.savenames <- unname(sapply(X = conf.savenames, FUN = function(x) Seurat:::ExtractField(string = x, delim = "\\.", field = 1)))
write.table(
  x = conf2.savenames, 
  file = paste0(getwd(), "/analysis_data/bipolar/bipolar_conf2.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
