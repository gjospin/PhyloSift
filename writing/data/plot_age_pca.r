#!/usr/bin/env Rscript

aller<-read.table("all.trans",sep=",")
meta <- read.table("human_micro_meta.csv",sep="\t",head=T)
colors = rainbow(60*12, end=0.7, alpha=0.5)

library(plotrix)

color_legend <- function(x, y, xlen, ylen, main, tiks){
    text(x, y+2*ylen/3, main, adj=c(0,0), cex=0.85)
    color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.7)
}


pdf("age_pca.pdf",width=5,height=5)
plot(aller$V2,aller$V3,pch=16,col=colors[trunc(1 + 12*meta$Age.of.host..years.[match(aller$V1,as.integer(meta$MG.RAST.ID))])],xlab="PC1 (59.9%)",ylab="PC2 (15.2%)")

color_legend( -6.5, 3.5, 3.5, 1.5, "age in years:", c(0,15,30,45,60))
dev.off()

