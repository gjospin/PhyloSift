#!/usr/bin/Rscript
#
# (c) 2011 Aaron Darling
# Very simplistic script to read a tree, find any leaf nodes with a name
# like 1_at_*, and report their nearest neighbors 
#
# function to remove tree nodes with more than some num descendants
depthfilter <- function(ttt, tnodes){
  filterdepth <- 60 # set this to the max descendants
	j <- 0
	filtnodes <- vector()
	for(i in 1:length(tnodes)){
		tedge <- which(ttt$edge[,2]==tnodes[i])
		tttr <- root(ttt, node=ttt$edge[tedge,1], resolve.root=TRUE)
		rootedge <- which( tttr$edge[,2]==max(tttr$edge[,2]) )
		tbalance <- balance(tttr)
		if(min(tbalance[tttr$edge[rootedge,1] - tttr$Nnode,]) < filterdepth){
			#filtnodes[j] <- tnodes[i]
      filtnodes <- c(filtnodes,tnodes[i])
		}
	}
	filtnodes
}



printNeighborName <- function(treefile){
	# read the tree
	ttt <- read.tree(treefile)
	# find the target seq
	tnodes <- grep('\\d+_at_', ttt$tip.label, perl=TRUE)

  	#filter nodes
  	#tnodes<-depthfilter(ttt,tnodes)  
	
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


