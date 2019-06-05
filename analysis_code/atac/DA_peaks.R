library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
atac <- readRDS(args[1])
Idents(atac) <- 'predicted.id'
confident <- subset(atac, subset = prediction.score.max > 0.5)

da_peaks <- function(obj, ident.1, ident.2) {
  cts <- GetAssayData(object = obj, assay = 'peaks', slot = 'data')
  cts_sum <- Matrix::colSums(cts)
  ob_a <- Matrix::rowMeans(cts[, WhichCells(obj, idents = ident.1)])
  ob_b <- Matrix::rowMeans(cts[, WhichCells(obj, idents = ident.2)])
  ob_diff <- sort(ob_a-ob_b)
  return(ob_diff)
}

pv.sst <- da_peaks(confident, 'Sst', 'Pvalb')
vip.lamp5 <- da_peaks(confident, "Vip", "Lamp5")

pv.up <- head(pv.sst, 1000)
sst.up <- tail(pv.sst, 1000)
lamp5.up <- head(vip.lamp5, 1000)
vip.up <- tail(vip.lamp5, 1000)

save_bed <- function(pks, filename) {
  writeLines(sapply(strsplit(pks, '-'), function(x) paste(x[1], x[2], x[3], sep='\t')),
                    con = filename)
}

save_bed(names(pv.up), "analysis_data/atac/da_peaks/pv.bed")
save_bed(names(sst.up), "analysis_data/atac/da_peaks/sst.bed")
save_bed(names(vip.up), "analysis_data/atac/da_peaks/vip.bed")
save_bed(names(lamp5.up), "analysis_data/atac/da_peaks/lamp5.bed")
