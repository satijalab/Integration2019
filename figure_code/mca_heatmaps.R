suppressMessages(library(Seurat))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

#args <- c("~/Projects/muir/seurat_objects/integrated_mca.rds")
mca <- readRDS(file = args[1])

DefaultAssay(object = mca) <- 'RNA'
mca <- SubsetData(object = mca, max.cells.per.ident = 200)
markers <- FindAllMarkers(object = mca, assay = "RNA", test.use = "LR", latent.vars = "tech", only.pos = TRUE)

markers.to.plot <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC) %>% pull(gene)
mca.subset <- SubsetData(object = mca, max.cells.per.ident = 200)  
DefaultAssay(object = mca.subset) <- 'RNA'
mca.techs <- SplitObject(object = mca.subset, split.by = "tech")
for(i in 1:length(x = mca.techs)) {
  mca.techs[[i]] <- ScaleData(object = mca.techs[[i]], features = markers.to.plot)
}

alldata <- matrix(nrow = length(x = markers.to.plot), ncol = 0)
for(i in 1:length(x = mca.techs)) {
  idents <- sort(x = names(x = which(x = table(Idents(object = mca.techs[[i]])) > 1)))
  newdata <- sapply(X = idents, function(x) {
    rowMeans(x = GetAssayData(object = mca.techs[[i]], slot = "scale.data", assay = "RNA")[, WhichCells(object = mca.techs[[i]], idents = x)])
    })[markers.to.plot, ]
  colnames(x = newdata) <- paste(idents, names(x = mca.techs)[i], sep = "_")
  alldata <- cbind(alldata, newdata)
}
rownames(alldata) <- make.unique(names = rownames(x = alldata))
alldata <- alldata[, sort(x = colnames(x = alldata))]
mca.avg <- CreateSeuratObject(counts = alldata)
mca.avg[["RNA"]]@scale.data <- alldata

hm.theme <- NoLegend() + 
  theme(axis.text.y = element_text(size = 0)) + 
  theme(axis.text.x = element_text(size = 0))

p1 <- DoHeatmap(
  object = mca.avg, 
  features = markers.to.plot, 
  disp.max = 4, 
  cells = grep(pattern = "10X", x = Cells(object = mca.avg), value = TRUE),
  label = FALSE) + hm.theme + ggtitle("Tabula Muris - 10X")
p2 <- DoHeatmap(
  object = mca.avg,
  features = markers.to.plot,
  disp.max = 4,
  cells = grep(pattern = "uWell", x = Cells(object = mca.avg), value = TRUE),
  label = FALSE) + hm.theme + ggtitle("Mouse Cell Atlas - Microwell-seq")
p3 <- DoHeatmap(
  object = mca.avg,
  features = markers.to.plot, 
  disp.max = 4,
  cells = grep(pattern = "FACS", x = Cells(object = mca.avg), value = TRUE),
  label = FALSE) + hm.theme + ggtitle("Tabula Muris - SMARTSeq2")

p.list <- list(p1, p2, p3)

saveRDS(object = p.list, file = args[2])

