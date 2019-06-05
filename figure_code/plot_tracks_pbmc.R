library(rtracklayer)
library(Gviz)
library(biomaRt)


b.pro <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/Bcellprogenitor.bw")
b.pre <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/pre-Bcell.bw")
nk <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/NKdim.bw")
cd8.naive <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD8Naive.bw")
cd8.memory <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD8Memory.bw")
cd8.effector <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD8effector.bw")
cd4.naive <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD4Naive.bw")
cd4.memory <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD4Memory.bw")
cd16.mono <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD16Monocytes.bw")
cd14.mono <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/CD14Monocytes.bw")
dc <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/Dendriticcell.bw")
pdc <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac_body_2kb/pDC.bw")
merged <- import("analysis_data/atac/10x/celltypes/single_cell/pbmc_all.bw")

bulk.cd8 <- import("raw_data/corces_atac/bigwig/CD8-10.merge.s20.w150sw.bw")
bulk.cd4 <- import("raw_data/corces_atac/bigwig/CD4-9.merge.s20.w150sw.bw")
bulk.b <- import("raw_data/corces_atac/bigwig/Bcell-13.merge.s20.w150sw.bw")
bulk.nk <- import("raw_data/corces_atac/bigwig/Nkcell-11.merge.s20.w150sw.bw")
bulk.mono <- import("raw_data/corces_atac/bigwig/Mono-7.merge.s20.w150sw.bw")

b.pro <- keepSeqlevels(b.pro, paste0('chr', 1:22), pruning.mode = 'coarse')
b.pre <- keepSeqlevels(b.pre, paste0('chr', 1:22), pruning.mode = 'coarse')
nk <- keepSeqlevels(nk, paste0('chr', 1:22), pruning.mode = 'coarse')
cd8.naive <- keepSeqlevels(cd8.naive, paste0('chr', 1:22), pruning.mode = 'coarse')
cd8.memory <- keepSeqlevels(cd8.memory, paste0('chr', 1:22), pruning.mode = 'coarse')
cd8.effector <- keepSeqlevels(cd8.effector, paste0('chr', 1:22), pruning.mode = 'coarse')
cd4.naive <- keepSeqlevels(cd4.naive, paste0('chr', 1:22), pruning.mode = 'coarse')
cd4.memory <- keepSeqlevels(cd4.memory, paste0('chr', 1:22), pruning.mode = 'coarse')
cd16.mono <- keepSeqlevels(cd16.mono, paste0('chr', 1:22), pruning.mode = 'coarse')
cd14.mono <- keepSeqlevels(cd14.mono, paste0('chr', 1:22), pruning.mode = 'coarse')
dc <- keepSeqlevels(dc, paste0('chr', 1:22), pruning.mode = 'coarse')
pdc <- keepSeqlevels(pdc, paste0('chr', 1:22), pruning.mode = 'coarse')
merged <- keepSeqlevels(merged, paste0('chr', 1:22), pruning.mode = 'coarse')

bulk.cd8 <- keepSeqlevels(bulk.cd8, paste0('chr', 1:22), pruning.mode = 'coarse')
bulk.cd4 <- keepSeqlevels(bulk.cd4, paste0('chr', 1:22), pruning.mode = 'coarse')
bulk.b <- keepSeqlevels(bulk.b, paste0('chr', 1:22), pruning.mode = 'coarse')
bulk.nk <- keepSeqlevels(bulk.nk, paste0('chr', 1:22), pruning.mode = 'coarse')
bulk.mono <- keepSeqlevels(bulk.mono, paste0('chr', 1:22), pruning.mode = 'coarse')

AT <- GenomeAxisTrack()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(12)

plot_region <- function(region, max=1000, max.bulk=max) {
  bm <- BiomartGeneRegionTrack(genome = "hg19", 
                               chromosome = region[[1]],
                               start = region[[2]],
                               end = region[[3]], 
                               name = "ENSEMBL",
                               collapseTranscripts = 'longest',
                               shape = 'box',
                               lex = 1/4,
                               fontsize.group=6)
  
  prob <- DataTrack(range = b.pro, genome = "hg19",
                       fill.mountain = c(cols[[1]], cols[[1]]), col.mountain = cols[[1]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "B cell progenitor", ylim = c(0, max))
  preb <- DataTrack(range = b.pre, genome = "hg19",
                       fill.mountain = c(cols[[2]], cols[[2]]), col.mountain = cols[[2]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "pre-B cell", ylim = c(0, max))
  sc.nk <- DataTrack(range = nk, genome = "hg19",
                       fill.mountain = c(cols[[3]], cols[[3]]), col.mountain = cols[[3]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "NK", ylim = c(0, max))
  sc.cd8.naive <- DataTrack(range = cd8.naive, genome = "hg19",
                       fill.mountain = c(cols[[4]], cols[[4]]), col.mountain = cols[[4]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD8 naive", ylim = c(0, max))
  sc.cd8.mem <- DataTrack(range = cd8.memory, genome = "hg19",
                       fill.mountain = c(cols[[5]], cols[[5]]), col.mountain = cols[[5]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD8 memory", ylim = c(0, max))
  sc.cd8.eff <- DataTrack(range = cd8.effector, genome = "hg19",
                       fill.mountain = c(cols[[6]], cols[[6]]), col.mountain = cols[[6]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD8 effector", ylim = c(0, max))
  sc.cd4.naive <- DataTrack(range = cd4.naive, genome = "hg19",
                       fill.mountain = c(cols[[7]], cols[[7]]), col.mountain = cols[[7]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD4 naive", ylim = c(0, max))
  sc.cd4.mem <- DataTrack(range = cd4.memory, genome = "hg19",
                       fill.mountain = c(cols[[8]], cols[[8]]), col.mountain = cols[[8]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD4 memory", ylim = c(0, max))
  sc.cd16.mono <- DataTrack(range = cd16.mono, genome = "hg19",
                       fill.mountain = c(cols[[9]], cols[[9]]), col.mountain = cols[[9]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD16+ monocyte", ylim = c(0, max))
  sc.cd14.mono <- DataTrack(range = cd14.mono, genome = "hg19",
                       fill.mountain = c(cols[[10]], cols[[10]]), col.mountain = cols[[10]],
                       type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "CD14+ monocyte", ylim = c(0, max))
  sc.dc <- DataTrack(range = dc, genome = "hg19",
                        fill.mountain = c(cols[[11]], cols[[11]]), col.mountain = cols[[11]],
                        type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "DC", ylim = c(0, max))
  sc.pdc <- DataTrack(range = pdc, genome = "hg19",
                        fill.mountain = c(cols[[12]], cols[[12]]), col.mountain = cols[[12]],
                        type = "polygon", chromosome = region[[1]], 
                       start = region[[2]], end = region[[3]],
                       name = "pDC", ylim = c(0, max))
  sc.merged <- DataTrack(range = merged, genome = "hg19",
                        fill.mountain = c('black', 'black'), col.mountain = 'black',
                        type = "polygon", chromosome = region[[1]], 
                        start = region[[2]], end = region[[3]],
                        name = "All cells", ylim = c(0, max))
  
  blk.b <- DataTrack(range = bulk.b, genome = "hg19",
                     fill.mountain = c(cols[[2]], cols[[2]]), col.mountain = 'black',
                     type = "polygon", chromosome = region[[1]], 
                     start = region[[2]], end = region[[3]],
                     lwd = 1/2,
                     name = "B bulk", ylim = c(0, max.bulk))
  blk.nk <- DataTrack(range = bulk.nk, genome = "hg19",
                      fill.mountain = c(cols[[3]], cols[[3]]), col.mountain = 'black',
                      type = "polygon", chromosome = region[[1]], 
                      lwd = 1/2,
                      start = region[[2]], end = region[[3]],
                      name = "NK bulk", ylim = c(0, max.bulk))
  blk.cd8 <- DataTrack(range = bulk.cd8, genome = "hg19",
                         fill.mountain = c(cols[[6]], cols[[6]]), col.mountain = 'black',
                         type = "polygon", chromosome = region[[1]], 
                       lwd = 1/2,
                         start = region[[2]], end = region[[3]],
                         name = "CD8 bulk", ylim = c(0, max.bulk))
  blk.cd4 <- DataTrack(range = bulk.cd4, genome = "hg19",
                       fill.mountain = c(cols[[8]], cols[[8]]), col.mountain = 'black',
                       type = "polygon", chromosome = region[[1]], 
                       lwd = 1/2,
                       start = region[[2]], end = region[[3]],
                       name = "CD5 bulk", ylim = c(0, max.bulk))
  blk.mono <- DataTrack(range = bulk.mono, genome = "hg19",
                       fill.mountain = c(cols[[10]], cols[[10]]), col.mountain = 'black',
                       type = "polygon", chromosome = region[[1]], 
                       lwd = 1/2,
                       start = region[[2]], end = region[[3]],
                       name = "Mono bulk", ylim = c(0, max.bulk))
  
  return(plotTracks(list(AT,
                         prob,
                         preb,
                         blk.b,
                         
                         sc.nk,
                         blk.nk,
                         
                         sc.cd8.naive,
                         sc.cd8.mem,
                         sc.cd8.eff,
                         blk.cd8,
                         
                         sc.cd4.naive,
                         sc.cd4.mem,
                         blk.cd4,
                         
                         sc.cd16.mono,
                         sc.cd14.mono,
                         blk.mono,
                         
                         sc.dc,
                         sc.pdc,
                         sc.merged,
                         bm),
                    from = region[[2]],
                    to = region[[3]],
                    scale = 0.2,
                    labelPos = "below",
                    cex = 1/2,
                    transcriptAnnotation = "symbol",
                    window = "auto",
                    cex.title = 1/2, fontsize = 5))
}

# regions
blc11b <- list("chr14", 99620230, 99747227)
blk <- list("chr8", 11343934, 11427243)
cd8 <- list("chr2", 87009034, 87037873)
cd4 <- list("chr12", 6893486, 6929992)
ms4a <- list("chr11", 60219830, 60238946)
cd3e <- list("chr11", 118174453, 118188955)
cd14 <- list("chr5", 140009839, 140015360)
lyz <- list("chr12", 69741059, 69748200)
gnly <- list("chr2", 85918895, 85940608)

png("figures/atac/pbmc/tracks/blc11b.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(blc11b, max=5000, max.bulk = 1000)
dev.off()

png("figures/atac/pbmc/tracks/cd8.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(cd8, max=5000, max.bulk = 2000)
dev.off()

png("figures/atac/pbmc/tracks/cd4.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(cd4, max=5000, max.bulk = 1000)
dev.off()

png("figures/atac/pbmc/tracks/ms4a.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(ms4a, max=10000, max.bulk = 2000)
dev.off()

png("figures/atac/pbmc/tracks/gnly.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(gnly, max=2000, max.bulk = 1000)
dev.off()

png("figures/atac/pbmc/tracks/cd3e.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(cd3e, max=5000, max.bulk = 2000)
dev.off()

png("figures/atac/pbmc/tracks/cd14.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(cd14, max=5000, max.bulk = 1000)
dev.off()

png("figures/atac/pbmc/tracks/lyz.png", height = 6, width = 1.5, units = 'in', res = 500)
plot_region(lyz, max=5000, max.bulk = 1000)
dev.off()
