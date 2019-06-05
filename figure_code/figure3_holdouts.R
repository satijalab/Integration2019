suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
library(viridis, quietly = TRUE)
library(ggalluvial)

args <- commandArgs(trailingOnly = TRUE)

devtools::load_all(args[8])
# args <- c("/home/butlera/Projects/muir_finaltest2/muir/analysis_data/pancreas/indrop4-alpha-projection-results.rds",
#          "/home/butlera/Projects/muir_finaltest2/muir/analysis_data/pancreas/all_accuracy.rds",
#          "/home/butlera/Projects/muir_finaltest2/muir/analysis_data/bipolar/all_accuracy.rds")

# indrop4 + alpha
proj.data <- readRDS(file = args[1])
filter.classes <- names(which(table(proj.data$seurat$true.label) < 2))

seurat.data <- as.data.frame(table(proj.data$seurat$true.label, proj.data$seurat$predicted.id))
colnames(seurat.data) <- c("actual", "prediction", "counts")
seurat.data %>% group_by(actual) %>% mutate(n.pred = sum(counts)) %>% mutate("Normalized Counts" = counts/n.pred) -> seurat.data
seurat.data %>% filter(!actual %in% filter.classes) -> seurat.data
ph1 <- ggplot(seurat.data, aes(x = prediction, y = reorder(actual))) + geom_tile(aes(fill = `Normalized Counts`)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0))  +
  ylab("Cell Type") + xlab("Prediction") +
  scale_fill_viridis(option = "D") 
  #ggtitle("Seurat") 

#proj.data$scmapcluster$predicted.id <- factor(x = proj.data$scmapcluster$predicted.id, levels = levels(x = seurat.data$prediction))
scmapcluster.data <- as.data.frame(table(proj.data$scmapcluster$true.label, proj.data$scmapcluster$predicted.id))
colnames(scmapcluster.data) <- c("actual", "prediction", "counts")
to.fill <- names(which(table(scmapcluster.data$prediction) == 0))
scmapcluster.data %>% group_by(actual) %>% mutate(n.pred = sum(counts)) %>% mutate("Normalized Counts" = counts/n.pred) -> scmapcluster.data
scmapcluster.data %>% filter(!actual %in% filter.classes) -> scmapcluster.data
ph2 <- ggplot(scmapcluster.data, aes(x = prediction, y = reorder(actual))) + geom_tile(aes(fill = `Normalized Counts`)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))  +
  ylab("Cell Type") + xlab("Prediction") +
  scale_fill_viridis(option = "D") 
  #ggtitle("scMap Cluster")

pfinal <- plot_grid(ph1 + NoLegend(), ph2 + NoLegend(), nrow = 1)

ggsave(plot = pfinal, filename = args[4], height = 7, width = 12)


seurat.predictions <- proj.data$seurat[, c("predicted.id", "prediction.score.max", "true.label")]
colnames(seurat.predictions) <- c("Seurat_V3_Prediction", "Seurat_V3_Score", "True_Label")
scmapcluster.predictions <- proj.data$scmapcluster[, c("predicted.id", "siml")]
colnames(scmapcluster.predictions) <- c("scMapCluster_Prediction", "scMapCluster_siml")

prediction.df <- merge(seurat.predictions, scmapcluster.predictions, by = "row.names")
colnames(prediction.df)[1] <- "cell"

cat("Supplementary Table 3: Seurat and scMapCluster cell type predictions for indrop 4 when beta cells left out of reference \n\n", file = args[7])
write.table(x = prediction.df, file = args[7], quote = FALSE, append = TRUE, sep = ",", row.names = FALSE)

# accuracy boxplots
pancreas.accuracy <- readRDS(file = args[2])
bipolar.accuracy <- readRDS(file = args[3])
pancreas.accuracy$dataset <- "pancreas"
bipolar.accuracy$dataset <- "bipolar"

all.accuracy <- suppressWarnings(rbind(pancreas.accuracy, bipolar.accuracy))
all.accuracy$method <- factor(all.accuracy$method, levels = rev(c("seurat", "scmap_cluster", "scmap_cell")))

pb <- ggplot(data = all.accuracy, mapping = aes(x = dataset, y = accuracy, fill = method)) + 
  geom_boxplot() + 
  ylab("Classification Accuracy") +
  xlab("") + coord_flip()

ggsave(plot = pb, filename = args[5], height = 4, width = 10)

# score distributions

proj.data$seurat$predicted.id.match <- proj.data$seurat$predicted.id
proj.data$seurat$predicted.id.match[proj.data$seurat$predicted.id.match == "Unknown"] <- "alpha"

proj.data$seurat$correct.prediction <- proj.data$seurat$predicted.id.match == proj.data$seurat$true.label
proj.data$seurat$correct.prediction[proj.data$seurat$correct.prediction == TRUE] <- "Correct"
proj.data$seurat$correct.prediction[proj.data$seurat$correct.prediction == FALSE] <- "Incorrect"
proj.data$seurat$correct.prediction <- factor(x = proj.data$seurat$correct.prediction, levels = c("Incorrect", "Correct"))

p.scores <- ggplot(data = proj.data$seurat, aes(fill = correct.prediction, x = prediction.score.max)) + 
  geom_histogram(alpha = 0.5, position="identity") +
  theme(legend.title = element_blank()) +
  xlab(label = "Prediction Score") + coord_flip()

ggsave(plot = p.scores, filename = args[6], height = 4, width = 10)

# p2 <- ggplot(data = proj.data$seurat, aes(color = correct.prediction, x = prediction.score.max)) + 
#   geom_histogram(fill = "white", alpha = 0.5, position="identity") +
#   theme(legend.title = element_blank()) +
#   xlab("Prediction Score")+
#   ggtitle("Option 2")