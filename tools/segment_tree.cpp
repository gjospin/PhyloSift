/******************************************************************************
** segment_tree.cpp
**   finds groups of leaf nodes in a binary tree where all 
**   connecting branches are less than some maximum length
**
**   input: binary tree with branch lengths in newick format
**   input: a maximum branch length
**
**   output: a tab-delimited table where each line contains a group of taxa
**
**   notes: the algorithm used here is based on Union-find, so compute 
**          time is linear in the number of tree nodes (or very nearly so)
**
** (c) 2012 Aaron Darling
** Licensed under the terms of the GPLv3
*******************************************************************************/

#include <iostream>
#include <PhyloTree.h>
#include <sstream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv){

	if(argc != 3){
		cerr << "Usage: segment_tree <tree> <maximum branch length>\n";
		return -1;
	}
	// read ref tree
	string reftreefile(argv[1]);
	ifstream reftreein(reftreefile.c_str());
	if(!reftreein.is_open()){
		cerr << "Unable to read file " << reftreefile << endl;
		return -2;
	}

	PhyloTree< TreeNode > reftree;
	reftree.readTree( reftreein );
	PhyloTree< TreeNode > original_tree = reftree;
		
	double max_length = atof(argv[2]);

	// cut all branches that have length > max length
	for(int nodeI=0; nodeI<reftree.size(); nodeI++){
		if(reftree[nodeI].distance<max_length)
			continue;
		// branch too long, cut parent->child link
		for(int pI=0; pI<reftree[nodeI].parents.size(); pI++){
			int parent = reftree[nodeI].parents[pI];
			for(int cI=0; cI<reftree[parent].children.size(); cI++){
				if(reftree[parent].children[cI]!=nodeI)
					continue;
				reftree[parent].children.erase(reftree[parent].children.begin() + cI);
				cI--;
			}
		}
		// cut child->parent link
		reftree[nodeI].parents.clear();
	}

	// find groups of nodes that are still connected with a union-find algorithm
	vector<int> node_groups( reftree.size() );
	map<int,int> group_sizes;
	map<int, vector<int> > groups;
	for(int nodeI=0; nodeI<reftree.size(); nodeI++){
		if(reftree[nodeI].parents.size()==0){
			node_groups[nodeI] = nodeI;
			group_sizes[nodeI] = 0;	// init this for later
			groups[nodeI] = vector<int>(0);
			continue;
		}
		// find node ancestor
		int ancestor = reftree[nodeI].parents[0];
		while(reftree[ancestor].parents.size()>0){
			ancestor = reftree[ancestor].parents[0];
		}
		node_groups[nodeI] = ancestor;

		// shortcut the tree to the ancestor to save time later (union-find style)
		int walker = reftree[nodeI].parents[0];
		while(reftree[walker].parents.size()>0){
			int tmp = reftree[walker].parents[0];
			reftree[walker].parents[0] = ancestor;
			walker = tmp;
		}
	}

	// count the number of leaf taxa in each group
	for(int nodeI=0; nodeI<reftree.size(); nodeI++){
		if(original_tree[nodeI].children.size()==0){
			group_sizes[ node_groups[nodeI] ]++;
			groups[ node_groups[nodeI] ].push_back(nodeI);
		}
	}

	// print out names of leaf taxa in each group with > 2 leaves
	map<int,int>::iterator iter = group_sizes.begin();
	for(; iter != group_sizes.end(); iter++){
		if(iter->second < 3)
			continue;
		vector<int>& group = groups[iter->first];
		for(int gI=0; gI<group.size(); gI++){
			if(gI>0)
				cout << "\t";
			cout << reftree[ group[gI] ].name;
		}
		cout << "\n";
	}
	
	return 0;
}

