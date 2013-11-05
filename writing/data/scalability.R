#!/usr/bin/env Rscript

scala<-read.table("scalability.csv",header=T)

pdf("scalability.pdf",width=8,height=5)
par(mar=c(5, 4, 4, 4) + 0.1)
plot(scala$reads,scala$wallclock/3600,type="o",xlab="Number of Illumina reads from gut metagenome", ylab="Single CPU running time (hours)",ylim=c(1,10))
axis(4)
mtext("Peak memory usage (GB)",side=4,line=2)
lines(scala$reads,scala$Maxvmem,col=4,lwd=2)
points(scala$reads,scala$Maxvmem,col=4,pch=2)
legend("bottomright",legend=c("Running time (hours)","Peak memory usage (GB)"),col=c(1,4),lwd=c(1,2),pch=c(1,2))

dev.off()
