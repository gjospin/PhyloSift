#!/usr/bin/Rscript
#
# (c) 2011 Aaron Darling
# Very simplistic script to read a tree, find any leaf nodes with a name
# like 1_at_*, and report their nearest neighbors 
#

printNeighborName <- function(treefile){
	# read the tree
	ttt <- read.tree(treefile)
	# find the target seq
	tnodes <- grep("1_at_", ttt$tip.label)
	# compute all pairwise distances to other nodes
	# find the closest node to our target
	pairdists<-cophenetic(ttt)
	sapply(tnodes, minNeighbor, pairdists)
}

minNeighbor<- function(tnode, pairdists){
	which(pairdists[tnode,pairdists[tnode,]>0]==min(pairdists[tnode,pairdists[tnode,]>0]))
}


# get command-line options
options <- commandArgs(trailingOnly = T)
# need the ape library, download and install it if
# it's not already available
if(!library("ape", logical.return=TRUE)){
	install.packages("ape", repos="http://cran.r-project.org")
	library("ape")
}

# run printNeighborName on all the tree files
result<-unlist(lapply(lapply(options, printNeighborName), names))
write(result,file="",ncolumns=1)


