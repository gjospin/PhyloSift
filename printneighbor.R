#!/usr/bin/Rscript

printNeighborName <- function(treefile){
	# read the tree
	ttt <- read.tree(treefile)
	# find the target seq
	tnodes <- grep("1_at_", ttt$tip.label)
	# compute all pairwise distances to other nodes
	# find the closest node to our target
	pairdists<-cophenetic(ttt)
	lapply(tnodes, minNeighbor, pairdists)
}

minNeighbor<- function(tnode, pairdists){
	which(pairdists[tnode,pairdists[tnode,]>0]==min(pairdists[tnode,pairdists[tnode,]>0]))
}


options <- commandArgs(trailingOnly = T)
library("ape")

lapply(options, printNeighborName)


