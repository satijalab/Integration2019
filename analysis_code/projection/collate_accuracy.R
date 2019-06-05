suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

accuracy <- data.frame(query.dataset = character(), celltype.removed = character(), method = character(), accuracy = numeric())

for(i in 2:length(args)) {
  dat <- readRDS(file = args[i])
  accuracy <- rbind(accuracy, as.data.frame(dat))
}
accuracy %>% group_by(method, celltype.removed) %>% mutate(sd = sd(accuracy)) %>% mutate(mean.accuracy = mean(accuracy)) -> accuracy
saveRDS(accuracy, file = args[1])
