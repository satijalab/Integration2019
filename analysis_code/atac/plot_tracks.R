library(rtracklayer)
library(Gviz)
library(biomaRt)

excite.bw <- import("analysis_data/atac/clusters/l4.bw")
vip.bw <- import("analysis_data/atac/clusters/vip.bw")
sst.bw <- import("analysis_data/atac/clusters/sst.bw")
pvalb.bw <- import("analysis_data/atac/clusters/pvalb.bw")
lamp5.bw <- import("analysis_data/atac/clusters/lamp5.bw")
bulk.bw <- import("analysis_data/atac/clusters/pfc_all.bw")

excite.color <- "#FE61CF"
vip.color <- "#8F91FF"
sst.color <- "#00BB4B"
pvalb.color <- "#BE80FF"
lamp5.color <- "#DD8D00"

AT <- GenomeAxisTrack()

plot_region <- function(region, max=200) {
  bm <- BiomartGeneRegionTrack(genome = "mm9", 
                               chromosome = region[[1]],
                               start = region[[2]],
                               end = region[[3]], 
                               name = "ENSEMBL",
                               shape = 'box',
                               collapseTranscripts = 'longest',
                               cex = 1/4, lex = 1/4, fontsize.group=6)
  
  dTrack0 <- DataTrack(range = excite.bw, genome = "mm9",
                       type = "h", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "cluster0", ylim = c(0, max), col = excite.color)
  dTrack1 <- DataTrack(range = vip.bw, genome = "mm9",
                       type = "h", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "cluster1", ylim = c(0, max), col = vip.color)
  dTrack2 <- DataTrack(range = lamp5.bw, genome = "mm9",
                       type = "h", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "cluster2",  ylim = c(0, max), col = lamp5.color)
  dTrack3 <- DataTrack(range = sst.bw, genome = "mm9",
                       type = "h", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "cluster4", ylim = c(0, max), col = sst.color)
  dTrack4 <- DataTrack(range = pvalb.bw, genome = "mm9",
                       type = "h", chromosome = region[[1]],
                       start = region[[2]], end = region[[3]],
                       name = "cluster7", ylim = c(0, max), col = pvalb.color)
  dTrack5 <- DataTrack(range = bulk.bw, genome = 'mm9',
                       type = 'h', chromosome = region[[1]],
                       start = region[[2]], end = region[[3]],
                       name = 'bulk', ylim = c(0, max), col = 'black')
  return(plotTracks(list(AT, dTrack0, dTrack1, dTrack2, dTrack3, dTrack4, dTrack5, bm),
                    from = region[[2]],
                    to = region[[3]],
                    scale = 0.2,
                    labelPos = "below",
                    cex = 1/2,
                    transcriptAnnotation = "symbol",
                    window = "auto",
                    cex.title = 1/2, fontsize = 5))
}

neurod6 <- list("chr6", 55624918, 55635412)
pv <- list('chr15', 78018321, 78037234)
sst <- list("chr16", 23884468, 23900790)
gad2 <- list('chr2', 22467475, 22554443)
vip <- list('chr10', 4692696, 4710921)
lamp5 <- list('chr2', 135874572, 135896567)
id2 <- list('chr12', 25774650, 25786427)
lhx6 <- list('chr2', 35930945, 35965047)
synpr <- list('chr14', 14116095, 14447844)
nek7 <- list('chr1', 140370795, 140526503)
tac1 <- list('chr6', 7501840, 7516179)

png("figures/atac/pfc/neurod6.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(neurod6, max = 1500)
dev.off()

png("figures/atac/pfc/gad2.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(gad2, max = 1000)
dev.off()

png("figures/atac/pfc/pvalb.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(pv, max = 1000)
dev.off()

png("figures/atac/pfc/sst.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(sst, max = 1000)
dev.off()

png("figures/atac/pfc/vip.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(vip, max = 1000)
dev.off()

png('figures/atac/pfc/lamp5.png', height = 2, width = 1.5, res = 500, units = 'in')
plot_region(lamp5, max = 1000)
dev.off()

png("figures/atac/pfc/lhx6.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(lhx6, max = 600)
dev.off()

png("figures/atac/pfc/synpr.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(synpr, max = 500)
dev.off()

png("figures/atac/pfc/nek7.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(nek7, max = 500)
dev.off()

png("figures/atac/pfc/tac1.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(tac1, max = 500)
dev.off()

png("figures/atac/pfc/id2.png", height = 2, width = 1.5, res = 500, units = 'in')
plot_region(id2, max = 1000)
dev.off()

