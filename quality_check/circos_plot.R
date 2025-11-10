# install package and load the library
# install.packages("circlize")
library(circlize)

# sed the working directory
setwd("P:/Yutangchen/Fabian/Sci_Data_paper/source_data/figure_1")

# import chromosome information and make the required input file for circlize
myfai <- read.table("Tuka_dip_chr.fasta.fai", header = F, stringsAsFactors = F)

mychrinfo <- cbind.data.frame(myfai$V1, rep(0, nrow(myfai)), myfai$V2)
colnames(mychrinfo) <- c("chr", "start", "end")
# mychrinfo$chr <- factor(mychrinfo$chr, levels = mychrinfo$chr)
head(mychrinfo)

# reverse the order of chrs in h2
# mychrinfo$chr <- factor(mychrinfo$chr, levels = c(paste("chr", 1:8, "_h1", sep = ""), paste("chr", 8:1, "_h2", sep = "")))
h2 <- mychrinfo[9:16, ]
h2 <- h2[order(h2$chr, decreasing = T), ]
mychrinfo <- rbind.data.frame(mychrinfo[1:8, ], h2)
str(mychrinfo)
mychrinfo$end <- as.numeric(mychrinfo$end)

# import gene and repeat 
gene <- read.table("gene_dip.cov", header = F, stringsAsFactors = F)
gene <- gene[, c(1, 2, 3, 7)]
colnames(gene) <- c("chr", "start", "end", "prop")
gene$pos <- floor((gene$start + gene$end)/2)
gene$prop <- gene$prop*100

REPEAT <- read.table("TE_dip.cov", header = F, stringsAsFactors = F)
REPEAT <- REPEAT[, c(1, 2, 3, 7)]
colnames(REPEAT) <- c("chr", "start", "end", "prop")
REPEAT$pos <- floor((REPEAT$start + REPEAT$end)/2)
REPEAT$prop <- REPEAT$prop*100

# GC content
GC<- read.table("Tuka_dip_chr.nuc", header = F, stringsAsFactors = F)
GC <- GC[, c(1,2,3,5)]
colnames(GC) <- c("chr", "start", "end", "prop")
GC$pos <- floor((GC$start + GC$end)/2)

# short read alignment depth
depth <- read.table("Tuka_dip_chr_depth.bed", header = F, stringsAsFactors = F)
colnames(depth) <- c("chr", "start", "end", "cov")
depth$pos <- floor((depth$start + depth$end)/2)

# gap
gap <- read.table('Tuka_dip_chr_gap.bed', header = F, stringsAsFactors = F)
gap <- gap[, c(1, 2, 3)]
colnames(gap) <- c("chr", "start", "end")
gap$pos <- floor((gap$start + gap$end)/2)

# alignments
link <- read.table("h1_vs_h2_chr.bed", header = F, stringsAsFactors = F)
colnames(link) <- c("chr_h1", "start_h1", "end_h1", "chr_h2", "start_h2", "end_h2")
# only show alignments >= 10kb
link$distance_h1 <- link$end_h1 - link$start_h1
link$distance_h2 <- link$end_h2 - link$start_h2
link <- link[link$distance_h1 >= 20000 & link$distance_h2 >= 20000, ]

# telomere
telo <- read.table("tuka_telo.tsv", header = T, stringsAsFactors = F)
telo <- telo[-c(33:35), ]
telo_bed <- NULL
for(i in 1:nrow(telo)){
  
  if(telo[i, 3] >= 100 & telo[i, 4] < 100){
    telo_bedx <- data.frame(chr = telo[i, 1],
                           start = 0,
                           end = 3000000)
  }else if(telo[i, 3] < 100 & telo[i, 4] >= 100){
    telo_bedx <- data.frame(chr = telo[i, 1],
                            start = mychrinfo[mychrinfo$chr == telo[i, 1], 3] - 3000000,
                            end = mychrinfo[mychrinfo$chr == telo[i, 1], 3])
  }
  telo_bed <- rbind.data.frame(telo_bed, telo_bedx)
}
telo_bed$end <- as.numeric(telo_bed$end)

# non-rep regions difined by k-mers
nonrep <- read.table("dip-non_repetitive.cov", header = F, stringsAsFactors = F)
nonrep <- nonrep[, c(1, 2, 3, 7)]
colnames(nonrep) <- c("chr", "start", "end", "prop")
nonrep$pos <- floor((nonrep$start + nonrep$end)/2)
nonrep$prop <- nonrep$prop*100

# define centromere as the lowest 5 non-rep windows in the middle
cent <- NULL
for(x in unique(nonrep$chr)){
  
  xx <- nonrep[nonrep$chr == x & nonrep$start >= (mychrinfo[mychrinfo$chr == x, 3]/2 - 10000000) & nonrep$end <= (mychrinfo[mychrinfo$chr == x, 3]/2 + 10000000), ]
  xx <- xx[order(xx$prop, decreasing = F), ]
  yy <- data.frame(chr = xx[1, 1], 
                   start = xx[1, 2] - 1000000,
                   end = xx[1, 3] + 1000000)
  cent <- rbind.data.frame(cent, yy)
  
}



#-------------------------------------------------------------------------------------------------------#

# make the circos plot

# add gap between chromosomes
circos.par(gap.after = c(rep(3, 15), 15), 
           track.margin = c(0.05, 0.05),
           start.degree = 90, 
           cell.padding = c(0, 0, 0, 0))

# initialize the circos plot
circos.initialize(sectors = mychrinfo$chr, xlim = mychrinfo[, 2:3])

# first track, highlight each assembly
circos.track(sector = mychrinfo$chr, ylim = c(0, 1), track.height = mm_h(3),
             bg.col = "white",
             bg.border = "black") 

for (x in mychrinfo$chr) {
  circos.text((mychrinfo[mychrinfo$chr == x, 3] - mychrinfo[mychrinfo$chr == x, 2])/2, mm_y(10), 
              labels = x, 
              sector.index = x,
              track.index = 1,
              facing = "bending.inside",
              cex = 1)
}

print("draw axis")

for (x in mychrinfo$chr) {
  circos.axis(h = "top", 
              major.at = seq(0, mychrinfo[mychrinfo$chr == x, 3], by = 50*10^6),
              labels = seq(0, mychrinfo[mychrinfo$chr == x, 3]/10^6, by = 50),
              sector.index = x,
              track.index = 1,
              labels.cex = 0.5,
              lwd = 1,
              minor.ticks = 0,
              major.tick.length = mm_y(0.2),
              direction = "outside",
              labels.facing = "outside",
              labels.pos.adjust = T)
  
}

print("draw gaps")

for (x in mychrinfo$chr) {
  # plot gap
  gapx <- gap[gap$chr == x, ]
  if(nrow(gapx) >= 1){

  for(j in 1:nrow(gapx)){    
      circos.rect(xleft = gapx[j, 2],
                  xright = gapx[j, 3],
                  ybottom = 0,
                  ytop = 1,
                  col = "black",
                  border = "black",
                  track.index = 1,
                  sector.index = x)
    }
  }
}

print("draw telomeres")

for (x in mychrinfo$chr) {

  telo_bedx <- telo_bed[telo_bed$chr == x, ]
  if(nrow(telo_bedx) >= 1){
    
    for(j in 1:nrow(telo_bedx)){    
      circos.rect(xleft = telo_bedx[j, 2],
                  xright = telo_bedx[j, 3],
                  ybottom = 0,
                  ytop = 1,
                  col = "deeppink3",
                  border = "deeppink3",
                  track.index = 1,
                  sector.index = x)
    }
  }
}

print("draw low non-rep regions")

for (x in mychrinfo$chr) {
  
  centx <- cent[cent$chr == x, ]
  if(nrow(centx) >= 1){
    
    for(j in 1:nrow(centx)){    
      circos.rect(xleft = centx[j, 2],
                  xright = centx[j, 3],
                  ybottom = 0,
                  ytop = 1,
                  col = "slateblue3",
                  border = "slateblue3",
                  track.index = 1,
                  sector.index = x)
    }
  }
}


# add label to each track
circos.text(x = 0, y = mm_y(5), labels = 'A', sector.index = mychrinfo$chr[1], track.index = 1,
            facing = "bending.inside", pos = 2, offset = 1.5, cex = 1)

# set gap size between tracks
set_track_gap(mm_h(0))

# second track, distribution of gene, heat map or hist gram

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), track.height = mm_h(10), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  genex <- gene[gene$chr == mychrinfo$chr[i], ]
  
  # plot bar
  circos.rect(xleft = genex[, 2],
              xright = genex[, 3],
              ybottom = rep(0, nrow(genex)),
              ytop = genex[, 4],
              col = "lightblue3",
              border = "lightblue3",
              track.index = 2,
              sector.index = mychrinfo$chr[i])
  
  nonrepx <- nonrep[nonrep$chr == mychrinfo$chr[i], ]
  
  circos.lines(x = nonrepx[, 5],
               y = nonrepx[, 4],
               col = "slateblue3",
               track.index = 2,
               sector.index = mychrinfo$chr[i])
  
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", mychrinfo$chr[1], 2),
  sector.index = mychrinfo$chr[1])

# add label
circos.text(x = 0, y = mm_y(12), labels = 'B', sector.index = mychrinfo$chr[1], track.index = 2,
            facing = "bending.inside", pos = 2, offset = 1.5, cex = 1)

# set gap size between tracks
set_track_gap(mm_h(3))

# third track, distribution of repeat

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), track.height = mm_h(8), bg.border = NA)

# use a line to show the abundance of repeats

for(i in 1:nrow(mychrinfo)){
  
  repeatx <- REPEAT[REPEAT$chr == mychrinfo$chr[i], ]
  
  circos.lines(x = repeatx[, 5],
               y = repeatx[, 4],
               col = "purple1",
               track.index = 3,
               sector.index = mychrinfo$chr[i])
  
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", mychrinfo$chr[1], 3),
  sector.index = mychrinfo$chr[1])

# add label
circos.text(x = 0, y = mm_y(10), labels = 'C', sector.index = mychrinfo$chr[1], track.index = 3,
            facing = "bending.inside", pos = 2, offset = 1.5, cex = 1)

# set gap size between tracks
set_track_gap(mm_h(3))

# fourth track, HiFi read alignment coverage

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), track.height = mm_h(8), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  depthx <- depth[depth$chr == mychrinfo$chr[i], ]
  
  # plot line 
  circos.lines(x = depthx[, 5],
               y = depthx[, 4],
               sector.index = mychrinfo$chr[i],
               col = "orange3",
               track.index = 4)
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", mychrinfo$chr[1], 4),
  sector.index = mychrinfo$chr[1])

# add lable
circos.text(x = 0, y = mm_y(7), labels = 'D', sector.index = mychrinfo$chr[1], track.index = 4,
            facing = "bending.inside", pos = 2, offset = 1.5, cex = 1)

for(i in 1:nrow(link)){
  circos.link(link[i, 1], link[i, 2], link[i, 4], link[i, 5], col = "darkseagreen")
}


circos.clear()



# # add labels to each track
# text(0, 0.95, labels = "A", pos = 2, offset = 1.5)
# text(0, 0.83, labels = "B", pos = 2, offset = 1.5)
# text(0, 0.67, labels = "C", pos = 2, offset = 1.5)
# text(0, 0.53, labels = "D", pos = 2, offset = 1.5)
# text(0, 0.38, labels = "E", pos = 2, offset = 1.5)
# 
# 
# text(-0.2, 0.2, labels = "A: chromosome", cex = 0.8, font = 2, pos = 4)
# text(-0.2, 0.1, labels = "B: gene density", cex = 0.8, font = 2, pos = 4)
# text(-0.2, 0, labels = "C: repeat density", cex = 0.8, font = 2, pos = 4)
# text(-0.2, -0.1, labels = "D: GC content", cex = 0.8, font = 2, pos = 4)
# text(-0.2, -0.2, labels = "E: short-read alignment depth", cex = 0.8, font = 2, pos = 4)
