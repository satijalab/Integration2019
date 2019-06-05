args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

allen <- read.table(file = "raw_data/allen_brain/mouse_VISp_2018-06-14_exon-matrix.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen$X <- NULL
allen <- as.matrix(x = allen)
genes <- read.table(file = "raw_data/allen_brain/mouse_VISp_2018-06-14_genes-rows.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(x = allen) <- make.unique(names = genes$gene_symbol)
meta.data <- read.csv(file = "raw_data/allen_brain/mouse_VISp_2018-06-14_samples-columns.csv",
                      row.names = 1, stringsAsFactors = FALSE)

al <- CreateSeuratObject(counts = allen, project = 'VISp', meta.data = meta.data, min.cells = 10)
low.q.cells <- rownames(x = meta.data[meta.data$class %in% c('Low Quality', 'No Class'), ])
ok.cells <- rownames(x = meta.data)[!(rownames(x = meta.data) %in% low.q.cells)]
al <- al[, ok.cells]
al <- NormalizeData(object = al)
al <- FindVariableFeatures(object = al, nfeatures = 2000)
al <- ScaleData(object = al)
al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
al <- RunUMAP(object = al, dims = 1:50, nneighbors = 5)
saveRDS(object = al, file = args[2])

# duplicate uppercase version for integration with ATAC 
rownames(x = allen) <- toupper(rownames(x = allen))
al <- CreateSeuratObject(counts = allen, project = 'VISp', meta.data = meta.data, min.cells = 10)
low.q.cells <- rownames(x = meta.data[meta.data$class %in% c('Low Quality', 'No Class'), ])
ok.cells <- rownames(x = meta.data)[!(rownames(x = meta.data) %in% low.q.cells)]
al <- al[, ok.cells]
al <- NormalizeData(object = al)
al <- FindVariableFeatures(object = al, nfeatures = 2000)
al <- ScaleData(object = al)
al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
al <- RunUMAP(object = al, dims = 1:50, nneighbors = 5)
saveRDS(object = al, file = "seurat_objects/allen_brain_upper.rds")

