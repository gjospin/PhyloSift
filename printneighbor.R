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
	tnodes <- grep('\\d+_at_', ttt$tip.label, perl=TRUE)

	#grep all special nodes

	# compute all pairwise distances to other nodes
	# find the closest node to our target
	pairdists<-cophenetic(ttt)

	#set all pairs with both members in tnodes to something big
	pairdists[tnodes,tnodes]<-max(pairdists)+1
	unlist(sapply(tnodes, minNeighbor, pairdists,ttt))
}

minNeighbor<- function(tnode, pairdists,ttt){
	min_tnode <- which(pairdists[tnode,pairdists[tnode,]>0]==min(pairdists[tnode,pairdists[tnode,]>0]))
	split_txt <-unlist(strsplit(ttt$tip.label[tnode],"_"))
	repeat_node<-as.numeric(split_txt[1])
	rep(min_tnode,repeat_node)     
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


