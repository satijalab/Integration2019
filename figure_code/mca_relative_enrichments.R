suppressMessages(library(Seurat))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)

#args <- c("seurat_objects/integrated_mca.rds")
mca <- readRDS(file = args[1])

# calculate relative enrichments
pt <- prop.table(x = table(mca$tissue, mca$tech), margin = 1)

relEnrich <- function(x = 0, pcount = 10) {
  cluster.cells <- Cells(x = mca)[which(mca$integrated_snn_res.4 == x)]
  tissue.ct <- table(mca$tissue[cluster.cells])
  tech.ct <- table(mca$tech[cluster.cells])
  expected.10x <- sum(pt[names(x = tissue.ct), , drop = FALSE][, "10X"] * tissue.ct)
  expected.FACS <- sum(pt[names(x = tissue.ct), , drop = FALSE][, "FACS"] * tissue.ct)
  expected.ratio <- (expected.10x + pcount) / (expected.FACS + pcount)
  actual.ratio <- (tech.ct[1] + pcount) / (tech.ct[2] + pcount)
  relative.enrichment <- (tech.ct[1] + pcount) / (expected.10x + pcount)
  return(relative.enrichment)
}

clusters <- names(x = which(x = table(mca$integrated_snn_res.4) > 0))
rel.enrich <- sapply(X = clusters, function(x) relEnrich(x = x, pcount = 20))
top.tissue <- sapply(X = clusters, function(x) {
  names(x = sort(x = table(mca$tissue[Cells(x = mca)[which(mca$integrated_snn_res.4 == x)]]), decreasing = TRUE)[1])
})

data <- data.frame(rel.enrich)
data$Tissue <- top.tissue
data <- data[order(data$rel.enrich), ]

p <- qplot(x = rel.enrich, data = data, geom = "density") + 
  xlim(0,2.2) + 
  xlab("10X vs FACS Cluster Enrichment")
p$data$size <- 1
p$data[p$data$Tissue == "marrow", "size"] <- 2.5
p <- p + geom_point(data = data, aes(x = rel.enrich, colour = Tissue), y = 0, size = p$data$size) +
  theme_classic()

saveRDS(object = p, file = args[2])

