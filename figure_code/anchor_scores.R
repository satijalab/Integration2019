library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[11])

# args <- c("~/Projects/muir/seurat_objects/celseq.rds",
#           "~/Projects/muir/seurat_objects/celseq2.rds",
#           "~/Projects/muir/seurat_objects/smartseq2.rds",
#           "~/Projects/muir/seurat_objects/fluidigmc1.rds",
#           "~/Projects/muir/seurat_objects/indrop.rds",
#           "~/Projects/muir/seurat_objects/integrated_pancreas.rds")

celseq <- readRDS(file = args[1])
celseq2 <- readRDS(file = args[2])
smartseq2 <- readRDS(file = args[3])
fluidigmc1 <- readRDS(file = args[4])
indrop <- readRDS(file = args[5])
pi <- readRDS(file = args[6])

all.pancreas <- c(list(celseq, celseq2, fluidigmc1, smartseq2), indrop)

for(ii in 1:length(x = all.pancreas)) {
  all.pancreas[[ii]][["celltype"]] <- pi$celltype
}

anchors <- Seurat:::AnnotateAnchors(
  object = pi, 
  toolname = "Integration",
  annotation = "celltype",
  object.list = all.pancreas
)
anchors$match <- anchors$cell1.celltype == anchors$cell2.celltype

fun_length <- function(x){
  return(data.frame(y=0.025, label= paste0("n=", length(x))))
}

anchors$match[anchors$match == TRUE] <- "Correct"
anchors$match[anchors$match == FALSE] <- "Incorrect"
anchors$match <- factor(x = anchors$match, levels = c("Incorrect", "Correct"))

p <- ggplot(data = anchors, mapping = aes(x = match, y = score)) + geom_boxplot() +
  ylab("Anchor Score") + xlab("") + # stat_summary(fun.data = fun_length, geom = "text", vjust = -2) +
  coord_flip() 

p2 <- ggplot(data = anchors, mapping = aes(x = match)) + geom_bar() + xlab("") + ylab("Number of Anchors") 

# m <- subset(x = anchors, cell1.celltype %in% c("alpha", "beta", "acinar", "ductal"))
# 
# p <- Seurat:::SingleExIPlot(data = m[, "score", drop = F], idents = m$cell1.celltype, split = m$match, type = "violin") + 
#   coord_flip() + 
#   ggtitle("Distribution of Anchor Scores") +
#   ylab("Anchor Score") +
#   xlab("") 
# 
# ntable <- data.frame(table(m$cell1.celltype, m$match))
# for(i in 1:nrow(ntable)) {
#   offset <- 0.15
#   if (ntable[i, "Var2"] == "FALSE") offset <- -.1
#   #print(paste(i, ((i-1)%%5), offset, ntable[i,3]))
#   p <- p + annotate(
#     geom = "text",
#     x = ((i - 1) %% 4) + 1 + offset,
#     y = 0.35,
#     label = paste0("n = ", ntable[i, 3]),
#     size = 3
#   )
# }

ggsave(plot = p, filename = args[7], width = 3, height = 8)
saveRDS(object = p, file = args[8])
ggsave(plot = p2, filename = args[9], width = 3, height = 8)
saveRDS(object = p2, file = args[10])
