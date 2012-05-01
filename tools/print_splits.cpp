/******************************************************************************
** print_splits.cpp
**   given a binary tree in newick format this produces a split encoding of
**   the tree topology, e.g. an encoding of the relationships among leaves on the tree
**   The resulting split alignment can be used as a set of constraints for tree
**   inference using tools like FastTree
**
**   input: file containing a binary tree with branch lengths in newick format
**   input: destination file name for the split alignment
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
#include <unordered_map>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace std;

void print_splits( PhyloTree< TreeNode >& reftree, ofstream& output );

unordered_map<string, string> other_map;

int main(int argc, char** argv){

	if(argc < 3){
		cerr << "Usage: printsplits <tree> <split output file>\n";
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
	cout << "The reference tree has " << reftree.size() << " nodes\n";
	

	// read map
	ofstream splitout(argv[2]);
	if(!splitout.is_open()){
		cerr << "Unable to write file " << argv[2] << endl;
		return -3;
	}

	print_splits( reftree, splitout );
	
	return 0;
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property< boost::edge_weight_t, double, boost::property< boost::edge_color_t, boost::default_color_type > > > Graph;
typedef boost::graph_traits< Graph >::vertex_descriptor Vertex;
typedef boost::graph_traits< Graph >::edge_descriptor GraphEdge;


typedef std::pair < size_t, size_t >Edge;
struct PhyloGraph {
	Graph g;
	Edge* edge_array;
	double* weights;
	int V;
	std::size_t E;
	boost::property_map < Graph, boost::edge_weight_t >::type w;
};


struct edge_filter {
  edge_filter(){}
  edge_filter(Edge e, Graph* g) : e(e), g(g){ }

  bool operator()(const GraphEdge& edge) const {
	size_t s = boost::source(edge, *g);
	size_t t = boost::target(edge, *g);
	if (( s == e.first && t == e.second ) ||
		(s == e.second && t == e.first) ){
		return false; 
	}
	return true;
  }
  Edge e;
  Graph* g;
};

typedef boost::filtered_graph<Graph, edge_filter > FilteredGraph;
typedef boost::graph_traits< Graph >::vertex_descriptor FilteredVertex;
typedef boost::graph_traits< Graph >::edge_descriptor FilteredGraphEdge;

class count_weight : public boost::default_dfs_visitor 
{
public:
	count_weight( boost::property_map < Graph, boost::edge_weight_t >::type& pm, double& s ) : pmap(pm), sum(s), recording(true) {};
	count_weight( const count_weight& cw ) : pmap( cw.pmap ), sum( cw.sum ), v(cw.v), recording(cw.recording) {};
	void start_vertex(Vertex v, const FilteredGraph& g){
		this->v = v;
	}
	void finish_vertex(Vertex v, const FilteredGraph& g){
		if(v == this->v)
			recording = false;	
	}
	void tree_edge(GraphEdge e, const FilteredGraph& g){
		if(recording)
			sum += pmap[e];
	}
	boost::property_map < Graph, boost::edge_weight_t >::type& pmap;
	double& sum;
	int v;
	bool recording;
};


class record_split : public boost::default_dfs_visitor 
{
public:
	record_split( boost::dynamic_bitset<>& bs, const vector<Vertex>& vertex_map ) : split(bs), vmap( vertex_map ), recording(true) {};
	record_split( const record_split& cw ) : split( cw.split ), v(cw.v), vmap( cw.vmap ), recording(cw.recording) {};
	void start_vertex(FilteredVertex v, const FilteredGraph& g){
		this->v = v;
	}
	void finish_vertex(FilteredVertex v, const FilteredGraph& g){
		if(recording && out_degree(v, g) <= 1)
			if(vmap[v] == -1){
				// not a leaf, probably internal telephone pole
			}else
				split.set(vmap[v]);
		// stop recording the split when we've gotten past the first connected component
		if(v == this->v)	
			recording = false;
	}
	boost::dynamic_bitset<>& split;
	const vector<Vertex>& vmap;
	int v;
	bool recording;
};

void make_graph( PhyloTree< TreeNode >& tree, PhyloGraph& pg ){
	// use boost's graph algorithms
	pg.V = tree.size();
	pg.E = tree.size()-1;
	pg.edge_array = new Edge[ tree.size() - 1 ];
	pg.weights = new double[ tree.size() - 1 ];
	size_t eI = 0;
	for( size_t vI = 0; vI < pg.V; ++vI )
	{
		if( tree[vI].parents.size() != 0 )
		{
			pg.edge_array[eI] = Edge( vI, tree[vI].parents[0] );
			pg.weights[eI] = tree[vI].distance;
			eI++;
		}
	}

	pg.g = Graph(pg.edge_array, pg.edge_array + pg.E, pg.V);

	// add edge weights
	pg.w = get(boost::edge_weight, pg.g);
	double *wp = pg.weights;
	boost::graph_traits < Graph >::edge_iterator e, e_end;
	for (boost::tie(e, e_end) = boost::edges(pg.g); e != e_end; ++e)
		pg.w[*e] = *wp++;

} 

int make_vertex_map(PhyloGraph& pg, vector< Vertex >& vmap ){
	int leafcount=0;
	for( int k=0; k < pg.V; k++)
		if(out_degree(k, pg.g) == 1){
			vmap[k]=leafcount++;
		}
	return leafcount;
}

void make_reverse_map( vector< Vertex >& vmap, vector<int>& reverse_map, int leafcount ){
	reverse_map.resize(leafcount);
	for(size_t i=0; i<vmap.size(); i++){
		if(vmap[i]!=-1){
			reverse_map[vmap[i]]=i;
		}
	}
}

void enumerate_splits( PhyloGraph& pg, vector< boost::dynamic_bitset<> >& splitlist, vector<Vertex>& vertex_map, vector<int>& reverse_map ){
	// make a leaf vertex map
	vertex_map.resize(pg.V,-1);
	int leafcount = make_vertex_map( pg, vertex_map );
	make_reverse_map( vertex_map, reverse_map, leafcount );
	for( size_t i=0; i < pg.E; i++ ){
		edge_filter ef(pg.edge_array[i], &pg.g);
		FilteredGraph g2(pg.g, ef);
		Vertex dfsroot = pg.edge_array[i].first;
		vector<boost::default_color_type> color(num_vertices( g2 ) );

		// get the split at this node
		boost::dynamic_bitset<> split(leafcount);
		record_split rs(split, vertex_map);
		depth_first_search(g2, rs, &color[0], dfsroot);
		splitlist.push_back(split);
	}
}

void print_splits( PhyloTree< TreeNode >& reftree, ofstream& output ){

//
// construct a boost graph of the tree
//
	PhyloGraph pg;	
	make_graph( reftree, pg );
//
// enumerate the splits present in the tree
//
	vector< boost::dynamic_bitset<> > pg_splitlist;
	vector<Vertex> pg_vertex_map;
	vector<int> reverse_map;
	enumerate_splits( pg, pg_splitlist, pg_vertex_map, reverse_map );
	cerr << "Enumerated " << pg_splitlist.size() << " splits\n";

//
// write the splits as an alignment
//

	vector<string> seqsplits( pg_splitlist[0].size() );
	for(size_t i=0; i < pg_splitlist.size(); i++){
		for(size_t j=0; j<pg_splitlist[i].size(); j++){
			seqsplits[j] +=  pg_splitlist[i].test(j) ? "1" : "0";
		}
	}
	for(size_t i=0; i<seqsplits.size(); i++){
		output << ">" << reftree[ reverse_map[i] ].name << endl;
		output << seqsplits[ i ] << endl;
	}
}


