suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])
input.file <- paste0(getwd(), "/raw_data/bipolar/bipolar.txt.gz")

bipolar <- suppressWarnings(fread(cmd = paste0("gzip -dc ", input.file)))
rows <- bipolar$V1
bipolar <- bipolar[,2:ncol(bipolar)]
bipolar <- as.matrix(x = bipolar)
rownames(bipolar) <- rows
metadata <- read.table(file = "analysis_data/bipolar_metadata.tsv", sep = '\t')
bipolar <- bipolar[, rownames(x = metadata)]

celltypes <- c('RBC',
              'Müller glia',
              'BC5A',
              'BC7',
              'BC6',
              'BC5C',
              'BC1A',
              'BC3B',
              'BC1B',
              'BC2',
              'BC5D',
              'BC3A',
              'BC5B',
              'BC4',
              'BC8/BC9',
              'Amacrine cells',
              'Doublets',
              'Doublets',
              'Doublets',
              'Rod  photoreceptors',
              'Doublets',
              'Cone photoreceptors',
              'Unknown',
              'Unknown',
              'Unknown',
              'Unknown')

metadata$celltype <- celltypes[metadata$CLUSTER]

bipolar <- CreateSeuratObject(counts = bipolar, meta.data = metadata, project = 'bipolar')
bipolar <- NormalizeData(object = bipolar, verbose = FALSE)
bipolar <- FindVariableFeatures(object = bipolar, selection.method = "vst", nfeatures = 2000)
bipolar <- ScaleData(object = bipolar, features = VariableFeatures(object = bipolar))
bipolar <- RunPCA(object = bipolar, dims.compute = 40, do.print = FALSE)
bipolar <- RunTSNE(object = bipolar, dims = 1:20, tsne.method = 'FIt-SNE')
saveRDS(object = bipolar, file = args[2])

