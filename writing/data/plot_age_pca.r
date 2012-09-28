#!/usr/bin/env Rscript

colors = rainbow(60*12, end=0.7, alpha=0.5)

library(plotrix)

color_legend <- function(x, y, xlen, ylen, main, tiks){
    text(x, y+2*ylen/3, main, adj=c(0,0), cex=0.85)
    color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.7)
}


aller<-read.table("all.trans",sep=",")
meta <- read.table("human_micro_meta.csv",sep="\t",head=T)

ages <- trunc(1 + 12*meta$Age.of.host..years.[match(aller$V1,as.integer(meta$MG.RAST.ID))])
logages <- log(ages)
scaler <- (max(ages[!is.na(ages)])/max(logages[!is.na(ages)]))
logages <- logages * scaler
legvec <- c(0,180,360,540,720)
legvec <- legvec / scaler
legvec <- exp(legvec)

pdf("age_pca_protein_phylosift.pdf",width=5,height=5)
plot(aller$V2,aller$V3,pch=16,col=colors[logages+1],xlab="PC1 (59.9%)",ylab="PC2 (15.2%)")
text(aller$V2,aller$V3,aller$V1,cex=0.5)
color_legend( -6.5, 3.5, 3.5, 1.5, "age in months:", trunc(legvec)-1)
dev.off()

aller<-read.table("16pca.trans",sep=",")
ages <- trunc(1 + 12*meta$Age.of.host..years.[match(aller$V1,as.integer(meta$MG.RAST.ID))])

pdf("age_pca_16s_phylosift.pdf",width=5,height=5)
plot(-aller$V2,aller$V3,pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC2 ")
text(-aller$V2,aller$V3,aller$V1,cex=0.5)
color_legend( -8, 5.5, 3.5, 1.5, "age in months:", trunc(legvec)-1)
dev.off()


pdf("age_pca_pc1pc3.pdf",width=5,height=5)
plot(aller$V2,aller$V4,pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC2 ")
text(aller$V2,aller$V4,aller$V1, cex=0.5)
dev.off()

meta16s <- read.table("metagenome_16s_weighted_unifrac_even_pc.txt",header=T,sep="\t")

ll<-nrow(meta16s)-2

ages <- trunc(1 + 12*meta$Age.of.host..years.[match(meta16s$pc.vector.number,as.integer(meta$MG.RAST.ID))])
logages <- log(ages)
scaler <- (max(ages[!is.na(ages)])/max(logages[!is.na(ages)]))
logages <- logages * scaler
pdf("age_pca_meta16_pc1pc2.pdf",width=5,height=5)
plot(meta16s$X1[1:ll],meta16s$X2[1:ll],pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC2 ")
text(meta16s$X1[1:ll],meta16s$X2[1:ll],meta16s$pc.vector.number, cex=0.5)
color_legend( -0.45, 0.35, 0.3, 0.1, "age in months:", trunc(legvec)-1)
dev.off()

pdf("age_pca_meta16_pc1pc3.pdf",width=5,height=5)
plot(meta16s$X1[1:ll],meta16s$X3[1:ll],pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC3 ")
text(meta16s$X1[1:ll],meta16s$X3[1:ll],meta16s$pc.vector.number, cex=0.5)
color_legend( 0.2, -0.4, 0.2, 0.1, "age in months:", trunc(legvec)-1)
dev.off()

pdf("age_pca_meta16_pc1pc4.pdf",width=5,height=5)
plot(meta16s$X1[1:ll],meta16s$X4[1:ll],pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC4 ")
text(meta16s$X1[1:ll],meta16s$X4[1:ll],meta16s$pc.vector.number, cex=0.5)
color_legend( 0.2, -0.25, 0.2, 0.1, "age in months:", trunc(legvec)-1)
dev.off()

amp16s <- read.table("amplicon_16s_weighted_unifrac_pc.txt",header=T,sep="\t")

ages <- trunc(1 + 12*meta$Age.of.host..years.[match(amp16s$pc.vector.number,meta$Metagenome)])
rastids <- meta$MG.RAST.ID[match(amp16s$pc.vector.number,meta$Metagenome)]
logages <- log(ages)
scaler <- (max(ages[!is.na(ages)])/max(logages[!is.na(ages)]))
logages <- logages * scaler
pdf("age_pca_amp16_pc1pc2.pdf",width=5,height=5)
plot(-amp16s$X1[1:ll],amp16s$X2[1:ll],pch=16,col=colors[logages+1],xlab="PC1 ",ylab="PC2 ")
text(-amp16s$X1[1:ll],amp16s$X2[1:ll],rastids[1:ll], cex=0.5)
color_legend( -0.85, 0.375, 0.5, 0.15, "age in months:", trunc(legvec)-1)
dev.off()

