suppressMessages(library(data.table))
library(biomaRt)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

input.file <- paste0(getwd(), "/raw_data/pancreas/fluidigmc1.csv.gz")

raw.data <- suppressWarnings(expr = fread(
  cmd = paste0("gzip -dc ", input.file),
  data.table = FALSE, 
  fill = FALSE, 
  showProgress = FALSE)
)
rownames(raw.data) <- raw.data[, 1]
raw.data <- raw.data[, -1]

# convert ensemble ids to gene names using biomaRt package
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 94)
genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               mart = ensembl)
rownames(x = genes) <- genes$ensembl_gene_id

# remove rows that don't have a gene name and enforce unique gene names
gene.names <- genes[rownames(x = raw.data), "external_gene_name"]
raw.data <- raw.data[which(x = !is.na(x = gene.names)), ]
rownames(x = raw.data) <- make.unique(names = gene.names[which(x = !is.na(x = gene.names))])
raw.data <- Matrix(data = as.matrix(x = raw.data), sparse = TRUE)
# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = raw.data, project = "FLUIDIGMC1")
# Not necessary to filter here because all cells have > 2500 genes/cell
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(object = seurat.object, verbose = FALSE,
                                      selection.method = "vst", nfeatures = 2000)
seurat.object[["tech"]] <- "fluidigmc1"
seurat.object[["replicate"]] <- "fluidigmc1"
Idents(object = seurat.object) <- "tech"

# Save rds file
saveRDS(object = seurat.object, file = args[2])
