library(methods)
library(Seurat)
library(microbenchmark)

mws <- readRDS(file = "seurat_objects/mca_mws.rds")
tm_droplet <- readRDS(file = "seurat_objects/mca_tm_droplet.rds")
mws <- FindVariableFeatures(mws, nfeatures = 2000)
tm_droplet <- FindVariableFeatures(tm_droplet, nfeatures = 2000)
features <- SelectIntegrationFeatures(obj.list = list(mws, tm_droplet), nfeatures = 2000, assay = 'RNA')
mws <- ScaleData(mws, features = features)
tm_droplet <- ScaleData(tm_droplet, features = features)

FindMAD <- function(baseline, query) {
  difference.matrix <- query - baseline
  return(mean(abs(difference.matrix)))
}

integrate <- function(obj1, obj2, ds1, ds2, features, n) {
  timing <- as.data.frame(microbenchmark(RunCCA(
    object1 = obj1,
    object2 = obj2, 
    features = features,
    group1.downsample = ds1,
    group2.downsample = ds2,
    verbose = FALSE, num.cc = 20, use.cpp = FALSE,
    renormalize = FALSE, rescale = FALSE
  ), times = n)
  )
  merged <- RunCCA(
    object1 = obj1,
    object2 = obj2, 
    features = features,
    group1.downsample = ds1,
    group2.downsample = ds2,
    verbose = FALSE, num.cc = 20, use.cpp = FALSE,
    renormalize = FALSE, rescale = FALSE
  )
  merged <- CosineCCA(merged)
  merged <- FindNeighbors(merged, cells1 = colnames(obj1), cells2 = colnames(obj2),
                          reduction = 'cca.cosine',
                          dims = 1:20, k = 200, eps = 5)
  merged <- FindMNN(merged, k = 5)
  top.features <- Seurat:::TopDimFeatures(merged, reduction = 'cca', dims = 1:20)
  merged <- FilterMNN(merged, k = 200, features = top.features)
  merged <- ScoreMNN(merged, sd = 1)
  merged <- FindIntegrationMatrix(merged)
  merged <- FindWeights(merged, k = 500, eps = 5, sd = 1, do.cpp = TRUE)
  merged <- IntegrateData(merged, do.cpp = TRUE)
  integrated <- GetAssayData(merged, assay = 'integrated', slot = 'data')
  return(list('timing' = timing$time, 'matrix' = integrated))
}

baseline <- integrate(obj1 = mws, obj2 = tm_droplet, ds1 = NULL, ds2 = NULL, features = features, n = 3)
results <- data.frame('timing' = baseline$timing$time, 'MAD' = c(0,0,0), downsample = c(1,1,1))
for(i in seq(0.1, 0.95, 0.05)) {
  print(i)
  ds1 <- as.integer(i * ncol(indrop1))
  ds2 <- as.integer(i * ncol(indrop2))
  turbo <- integrate(obj1 = indrop1, obj2 = indrop2, features = features, ds1 = ds1, ds2 = ds2, n = 3)
  mad <- FindMAD(baseline = baseline$matrix, query = turbo$matrix)
  for(x in 1:3) {
    results <- rbind(results, data.frame('timing' = turbo$timing[x], 'MAD' = mad, 'downsample' = i))
  }
}
write.table(results, "analysis_data/cca_benchmarks.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
