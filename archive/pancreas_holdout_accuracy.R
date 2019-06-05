suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scmap))
suppressMessages(library(scPred))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("/home/butlera/Projects/muir/analysis_data/pancreas/reference-without-smartseq2-alpha.rds",
#           "/home/butlera/Projects/muir/analysis_data/pancreas/query-without-smartseq2-alpha.rds",
#           "alpha", "celseq")

ref.integrated <- readRDS(file = args[1])
query <- readRDS(file = args[2])
cell.type.removed <- args[3]
query.name <- args[4]

accuracy_results <- data.frame(query.dataset = character(), celltype.removed = character(), method = character(), accuracy = numeric())

if(is.null(ncol(ref.integrated)) || is.null(ncol(query))) {
  saveRDS(accuracy_results, file = args[5])
  saveRDS(list(), file = args[6])
} else {
  # perform projection 
  seurat.projection <- ProjectCells(
    reference = ref.integrated, query = query, label.name = 'celltype', eps = 0,
    do.cpp = TRUE, verbose = TRUE, sd = 1, k.neighbors = 100, k.mnn = 10, k.filter = 50, 
    k.weights = 50, use.cosine = TRUE, dims = 1:30, normalize.with.lm = FALSE, score = FALSE, score2 = TRUE
  )
  
  seurat.projection$true.label <- as.vector(x = query$celltype)
  seurat.projection <- seurat.projection[order(seurat.projection$prediction.score.max), ]
  seurat.projection$predicted.id.save <- seurat.projection$predicted.id
  seurat.projection[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- cell.type.removed
  seurat.accuracy <- table(seurat.projection$predicted.id == seurat.projection$true.label)[2] / 
    sum(table(seurat.projection$predicted.id == seurat.projection$true.label) )
  seurat.projection[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- "Unassigned"
  
  # predict using scMAP
  sce <- SingleCellExperiment(assays = list(logcounts = ref.integrated[["RNA"]]@data, counts = ref.integrated[["RNA"]]@counts))
  colData(sce) <- S4Vectors::DataFrame(ref.integrated[[]])
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  assays(sce)$counts <- as.matrix(assays(sce)$counts)
  assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
  sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)
  
  sce.query <- SingleCellExperiment(assays = list(logcounts = query[["RNA"]]@data, counts = query[["RNA"]]@counts))
  colData(sce.query) <- S4Vectors::DataFrame(query[[]])
  rowData(sce.query)$feature_symbol <- rownames(sce.query)
  assays(sce.query)$counts <- as.matrix(assays(sce.query)$counts)
  assays(sce.query)$logcounts <- as.matrix(assays(sce.query)$logcounts)
  
  # scmap-cluster
  sce <- indexCluster(sce, cluster_col = 'celltype')
  scmapcluster.results <- scmapCluster(
    projection = sce.query, 
    index_list = list(yan = metadata(sce)$scmap_cluster_index), 
    threshold = -Inf  # force assignments 
  )
  scmapcluster.results <- data.frame(
    predicted.id = as.vector(scmapcluster.results$combined_labs), 
    siml = unname(scmapcluster.results$scmap_cluster_siml), 
    true.label = colData(sce.query)$celltype, stringsAsFactors = FALSE, row.names = rownames(colData(sce.query))
  )
  scmapcluster.results$siml[which(is.na(scmapcluster.results$siml))] <- 0
  scmapcluster.results <- scmapcluster.results[order(scmapcluster.results$siml), ]
  scmapcluster.results$predicted.id.save <- scmapcluster.results$predicted.id
  scmapcluster.results[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- cell.type.removed
  scmapcluster.accuracy <- table(scmapcluster.results$predicted.id == scmapcluster.results$true.label)[2] / 
    sum(table(scmapcluster.results$predicted.id == scmapcluster.results$true.label))
  scmapcluster.results[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- "Unassigned"
  
  # scmap-cell
  sce <- indexCell(sce)
  scmapcell.results <- scmapCell(projection = sce.query, index_list = list(yan = metadata(sce)$scmap_cell_index))
  scmapcell.clusters <- scmapCell2Cluster(scmapcell.results, list(as.character(colData(sce)$celltype)), threshold = -Inf, w = 1)
  scmapcell.results <- data.frame(predicted.id = as.vector(scmapcell.clusters$combined_labs), siml = unname(scmapcell.clusters$scmap_cluster_siml), 
                                  true.label = colData(sce.query)$celltype, stringsAsFactors = FALSE, row.names = rownames(colData(sce.query)))
  scmapcell.results$siml[which(is.na(scmapcell.results$siml))] <- 0
  scmapcell.results <- scmapcell.results[order(scmapcell.results$siml), ]
  scmapcell.results$predicted.id.save <- scmapcell.results$predicted.id
  scmapcell.results[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- cell.type.removed
  scmapcell.accuracy <- table(scmapcell.results$predicted.id == scmapcell.results$true.label)[2] / 
    sum(table(scmapcell.results$predicted.id == scmapcell.results$true.label))
  scmapcell.results[1:(length(which(query$celltype == cell.type.removed))), "predicted.id"] <- "Unassigned"
  
  # scPred
  # genes.use <- rownames(ref.integrated[["integrated"]])
  # genes.use <- intersect(genes.use, rownames(x = query))
  # Idents(object = ref.integrated) <- 'celltype'
  # train.data <- GetAssayData(object = ref.integrated, assay = "RNA", slot = "data")[genes.use, ]
  # scp <- eigenDecompose(t(x = as.matrix(x = train.data)), n = 30)
  # train.metadata <- ref.integrated[[]]
  # train.metadata$celltype <- factor(train.metadata$celltype)
  # scPred::metadata(object = scp) <- train.metadata
  # scp <- getFeatureSpace(object = scp, pVar = "celltype")
  # newData <- t(x = as.matrix(x = GetAssayData(object = query, assay = "RNA", slot = "data")[genes.use, ]))
  # scpred.predictions <- scPredict(object = scp, newData = newData, threshold = 0)
  # scpred.predictions$score <- apply(X = scpred.predictions[1:(ncol(x = scpred.predictions) - 1)], MARGIN = 1, max)
  # scpred.predictions$match <- scpred.predictions$predClass == scpred.predictions$realID
  # scpred.predictions$realID <- query$celltype
  # scpred.accuracy <- table(scpred.predictions$match)[2] / sum(table(scpred.predictions$match))
  # 
  accuracy_results <- data.frame(query.dataset = character(), celltype.removed = character(), method = character(), accuracy = numeric())
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = query.name, 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "seurat", 
                                                         accuracy = seurat.accuracy))
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = query.name, 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "scmap_cell", 
                                                         accuracy = scmapcell.accuracy))
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = query.name, 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "scmap_cluster", 
                                                         accuracy = scmapcluster.accuracy))
  # accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = query.name, 
  #                                                        celltype.removed = cell.type.removed, 
  #                                                        method = "scPred", 
  #                                                        accuracy = scpred.accuracy))
  
  saveRDS(accuracy_results, file = args[5])
  all.results <- list(seurat = seurat.projection, scmapcell = scmapcell.results, scmapcluster = scmapcluster.results)
  saveRDS(all.results, file = args[6])
}

