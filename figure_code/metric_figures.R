# metrics plot
library(tidyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(ggrepel, quietly = TRUE)


args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[13])

# 
# args <- c("/home/butlera/Projects/muir/analysis_data/integrated_pancreas_celltype_holdouts_seuratV3_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_bipolar_celltype_holdouts_seuratV3_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_pancreas_celltype_holdouts_seuratV2_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_bipolar_celltype_holdouts_seuratV2_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_pancreas_celltype_holdouts_mnnCorrect_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_bipolar_celltype_holdouts_mnnCorrect_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_pancreas_celltype_holdouts_scanorama_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_bipolar_celltype_holdouts_scanorama_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_pancreas_celltype_holdouts_none_metrics.rds",
#           "/home/butlera/Projects/muir/analysis_data/integrated_bipolar_celltype_holdouts_none_metrics.rds")

silhouettes <- data.frame(sil = numeric(), method = character(), dataset = character())
mixing.metric <- data.frame(mixing.metric = numeric(), method = character(), dataset = character())
local.struct <- data.frame(local.struct = numeric(), method = character(), dataset = character())

for (ii in args[1:10]) {
  method <- Seurat:::ExtractField(
    string = Seurat:::ExtractField(string = basename(ii), field = 5, delim = "_"),
    field = 1,
    delim = "\\."
  )
  dataset.name <- Seurat:::ExtractField(
    string = Seurat:::ExtractField(string = basename(ii), field = 2, delim = "_"),
    field = 1,
    delim = "\\."
  )
  metrics <- readRDS(file = ii)
  silhouettes <- rbind(
    silhouettes, 
    data.frame(sil = metrics$silhouette, method = method, dataset = dataset.name)
  )
  mixing.metric <- rbind(
    mixing.metric, 
    data.frame(mixing.metric = metrics$mixing.metric, method = method, dataset = dataset.name)
  )
  local.struct <- rbind(
    local.struct, 
    data.frame(local.struct = metrics$local.struct, method = method, dataset = dataset.name)
  )
}

mixing.metric$method <- factor(x = mixing.metric$method, levels = c("seuratV3", "seuratV2", "scanorama", "mnnCorrect", "none"))
local.struct$method <- factor(x = local.struct$method, levels = c("seuratV3", "seuratV2", "scanorama", "mnnCorrect", "none"))
silhouettes %>% filter(method != "none") -> silhouettes
silhouettes$method <- factor(x = silhouettes$method, levels = c("seuratV3", "seuratV2", "mnnCorrect", "scanorama"))

p1 <- ggplot(data = silhouettes, aes(x = dataset, y = sil, col = method)) + 
  geom_boxplot(outlier.size = 0.1)  +
  xlab("Method") + ylab("Silhoutte Metric") +
  scale_color_discrete(name = "Method")

mixing.metric %>% filter(method != "none") -> mixing.metric
local.struct %>% filter(method != "none") -> local.struct

# p2a <- ggplot(data = mixing.metric, aes(x = dataset, y = mixing.metric, col = method)) + 
#   geom_boxplot(outlier.size = 0.1) +
#   xlab("Dataset") +
#   ylab("Mixing Metric")
# p2b <- ggplot(data = local.struct, aes(x = dataset, y = local.struct, col = method)) + 
#   geom_boxplot(outlier.size = 0.1) +
#   xlab("Dataset") +
#   ylab("Local Structure Metric")
# 
# p2 <- list(p2a, p2b)

mixing.metric %>% group_by(method, dataset) %>% mutate(mixing.metric = mixing.metric) %>%
  summarize(sd.mixing = sd(mixing.metric), n = n(), mean.mixing = mean(mixing.metric)) %>%
  mutate(se.mixing = sd.mixing/sqrt(n)) -> mixing.stats
local.struct %>% group_by(method, dataset) %>%
  summarize(sd.ls = sd(local.struct), n = n(), mean.ls = mean(local.struct)) %>%
  mutate(se.ls = sd.ls/sqrt(n)) -> ls.stats

all.metrics <- merge(mixing.stats, ls.stats)
all.metrics$method <- factor(x = all.metrics$method, levels = c("seuratV3", "seuratV2", "mnnCorrect", "scanorama"))
# p2 <- ggplot(data = all.metrics, aes(x = mean.mixing, y = mean.ls)) + geom_point(aes(color = method, shape = dataset), size = 3) +
#   geom_errorbar(aes(ymin=mean.ls-sd.ls, ymax=mean.ls+sd.ls), width = 0.025) +
#   geom_errorbarh(aes(xmin = mean.mixing-sd.mixing, xmax=mean.mixing+sd.mixing))+
#   xlab("Mixing Metric") + ylab("Local Structure Metric")

p2 <- ggplot(data = all.metrics, aes(x = mean.mixing, y = mean.ls)) +
  geom_point(aes(color = method), size = 5) + facet_wrap(~dataset)+
  xlab("Mixing Metric") + ylab("Local Structure Metric") +
  geom_text_repel(aes(label = method), size = 5, box.padding = 0.9, segment.size = 0, seed = 80, min.segment.length = Inf) +
  ylim(c(0, max(ls.stats$mean.ls))) + 
  NoLegend()

saveRDS(object = p1, file = args[11])
saveRDS(object = p2, file = args[12])


