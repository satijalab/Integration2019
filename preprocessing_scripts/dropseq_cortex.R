args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])
library(Matrix)

dropseq <- Matrix::readMM(file = "raw_data/dropseq_cortex/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz")
dropseq <- as(Class = 'dgCMatrix', dropseq)

# copied from DropSeq package
strStartsWith = function(theString, thePrefix) {
  return(strtrim(theString, nchar(thePrefix)) == thePrefix)
}
loadSparseDgeNames=function(file) {
  conn=file(file, "r")
  genes=c()
  cell_barcodes=c()
  while(TRUE) {
    line=readLines(con=conn, n=1, ok=FALSE)
    if (!strStartsWith(line, '%')) {
      break
    }
    if (strStartsWith(line, '%%GENES\t')) {
      these_genes=strsplit(line, "\t", fixed=TRUE)[[1]]
      genes = append(genes, these_genes[2:length(these_genes)])
    } else if (strStartsWith(line, '%%CELL_BARCODES\t')) {
      these_cells=strsplit(line, "\t", fixed=TRUE)[[1]]
      cell_barcodes = append(cell_barcodes, these_cells[2:length(these_cells)])
    }
  }
  close(conn)
  return(list(genes=genes, cell_barcodes=cell_barcodes))
}

gene_cell_names <- loadSparseDgeNames("raw_data/dropseq_cortex/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz")
annot <- readRDS("raw_data/dropseq_cortex/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.cell_cluster_outcomes.RDS")
rownames(dropseq) <- make.unique(gene_cell_names$genes)
colnames(dropseq) <- gene_cell_names$cell_barcodes

ds <- CreateSeuratObject(counts = dropseq, meta.data = annot)
ds <- ds[, is.na(ds$reason)]
ds <- FindVariableFeatures(ds)
ds <- NormalizeData(ds)
ds <- ScaleData(ds)
ds <- RunPCA(ds, npcs = 100)
ds <- FindNeighbors(ds, dims = 1:50)
ds <- FindClusters(ds)
ds <- RunUMAP(ds, graph = "RNA_nn")
saveRDS(file = args[2], object = ds)
