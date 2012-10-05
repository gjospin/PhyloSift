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

void reconcile( PhyloTree< TreeNode >& reftree, string treefile, unordered_multimap<string, string>& gene_map, string output_fname  );

unordered_map<string, string> other_map;

/** An upper limit on the number of taxonomy mappings that an edge in the gene tree is allowed to have. 
 *  If more than this, the gene tree edge is ignored.  This is ugly, but effectively happens anyway since the mapping gets too diffuse.
 *  We need a real gene tree/species tree reconciliation algorithm to fix this.
 */
const int placement_limit = 30;

int main(int argc, char** argv){

	if(argc < 3){
		cerr << "Usage: readconciler <reference tree> <gene tree> <gene to species map> <mapping output file>\n";
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
	ifstream mapin(argv[3]);
	if(!mapin.is_open()){
		cerr << "Unable to read file " << argv[3] << endl;
		return -3;
	}
	unordered_multimap<string, string> gene_map;
	string line;
	while( getline(mapin, line) ){
		stringstream line_str(line);
		string gene;
		string species;
		getline(line_str, species, '\t');
		getline(line_str, species, '\t');
		getline(line_str, gene);
		gene_map.insert(make_pair(species, gene));
		other_map.insert(make_pair(gene,species));
	}

	// read & reconcile each of the read trees
	string genetreefile = argv[2];
	reconcile( reftree, genetreefile, gene_map, argv[4] );
	
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

void enumerate_splits( PhyloGraph& pg, vector< boost::dynamic_bitset<> >& splitlist, vector<Vertex>& vertex_map ){
	// make a leaf vertex map
	vertex_map.resize(pg.V,-1);
	int leafcount = make_vertex_map( pg, vertex_map );
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

void normalize_split( boost::dynamic_bitset<>& split, const vector< vector< int > >& splitmap, int target_size ){
	boost::dynamic_bitset<> newsplit( target_size );
	for( int i=0; i<splitmap.size(); i++ )
		for(int j=0; j<splitmap[i].size(); j++)
			newsplit.set( splitmap[i][j], split.test(i) );
	split = newsplit;
}

void reconcile( PhyloTree< TreeNode >& reftree, string treefile, unordered_multimap<string, string>& gene_map, string output_fname ){
	// read ref tree
	ifstream treein(treefile.c_str());
	if(!treein.is_open()){
		cerr << "Unable to read file " << treefile << endl;
		return;
	}

//
// read a tree with edge numberings from pplacer
// assume jplace format with treestring on second line
//
	string line;
	string treestring;
	getline( treein, treestring );
	stringstream treestr(treestring);
//	cout << "Trying to read " << treestring << endl;

	PhyloTree< TreeNode > tree;
	tree.readTree( treestr );
	cout << "The read tree has " << tree.size() << " nodes\n";
//
// remove edge numbers
// assume jplace format
//
	std::unordered_map<int,int> edgenum_map;
	for(int i=0; i<tree.size(); i++){
		size_t atpos = tree[i].name.find("{");
		size_t ratpos = tree[i].name.rfind("}");
		int edgenum = -1;		
		if( atpos == string::npos ){
			edgenum = atoi(tree[i].name.c_str());
		}else{
			edgenum = atoi(tree[i].name.substr(atpos+1, ratpos - atpos - 1).c_str());
//			cerr << "node " << i << " edgenum is " << tree[i].name.substr(atpos+1, ratpos - atpos - 1) << " name is " << tree[i].name.substr(0, atpos) << endl;
			tree[i].name = tree[i].name.substr(0, atpos);
		}
//		cerr << "mapping " << i << " to " << edgenum << "\n";
		edgenum_map.insert(make_pair(i,edgenum));
	}
//	cerr << "Done removing edge numbers\n";

//
// construct boost graphs of the trees
//
	PhyloGraph pg;	
	make_graph( tree, pg );

	PhyloGraph refpg;	
	make_graph( reftree, refpg );

//
// Phase 3: construct map to reference tree
//
// a) cut gene tree on each edge
// b) compute splits at cut point
// c) cut species tree on each edge
// d) determine which species tree split matches the gene tree split best
// e) write out the split match


// plan for later...
// c) compute PD on either side of cut point
// d) logical AND splits with reftree splits
// e) compute minimum spanning tree among remaining nodes
// f) compute PD of minimum spanning trees
// 	
	vector< boost::dynamic_bitset<> > pg_splitlist;
	vector<Vertex> pg_vertex_map;
	enumerate_splits( pg, pg_splitlist, pg_vertex_map );
	cout << "Done with gene tree splits\n";
	vector< boost::dynamic_bitset<> > ref_splitlist;
	vector<Vertex> ref_vertex_map;
	enumerate_splits( refpg, ref_splitlist, ref_vertex_map );

	// need a mapping from vertex numbers in refpg to vertex numbers in pg
	cout << "Making gene tree map\n";
	unordered_map< string, int > gtmap;
	for(int i=0; i<tree.size(); i++){
		if(tree[i].children.size()==0){
			gtmap.insert(make_pair(tree[i].name, i));
		}
	}
	cout << gtmap.size() << " genes mapped\n";
	
	cout << "Making species to gene tree map\n";
	vector< vector< int > > species_to_gene_map;	// maps split IDs in species tree to split IDs in gene tree
	for(int i=0; i<refpg.V; i++){
		if(ref_vertex_map[i]==-1)
			continue;
		// which genes does this species contain?
		pair< unordered_multimap<string,string>::iterator, unordered_multimap<string,string>::iterator> iter;
		iter = gene_map.equal_range(reftree[i].name);
		vector<int> curmap;
		if(iter.first ==iter.second){
			cerr << "Error no mapping found for " << reftree[i].name << endl;
		}
		for(; iter.first !=iter.second; iter.first++){			
			if( pg_vertex_map[ gtmap[iter.first->second] ] == -1 )
				continue;
			curmap.push_back( pg_vertex_map[ gtmap[iter.first->second] ] );
//			cout << "mapped ref " << reftree[i].name << "\t" << ref_vertex_map[i] << " to " << curmap.back() << endl;
//			cout << "reverse map to " << tree[gtmap[iter.first->second]].name << " and " << other_map[tree[gtmap[iter.first->second]].name] << endl;
		}
		// add a list of gene vertices for this species
		species_to_gene_map.push_back(curmap);
	}
	cout << species_to_gene_map.size() << " species mapped\n";
	cout << "rs.size() " << ref_splitlist[0].size() << endl;

	cout << "Finding best edges\n";
	ofstream mapout(output_fname.c_str());
	for( size_t i=0; i < pg.E; i++ ){
		// for each reftree edge, calculate mapping quality between this edge and reftree edges
		double scoresum = 0;
		double bestscore = 0;
		vector<double> maxscores;
//		cout << "ts1.count()\t" << pg_splitlist[i].count() << endl;
		if(pg_splitlist[i].count() == 1){
			size_t f = pg_splitlist[i].find_first();
			int qq=0;
			for(int abc=-1; abc<(int)f; qq++)
				if(pg_vertex_map[qq]!=-1)
					abc++;
//			cout << "gene tree " << other_map[ tree[qq].name ] << " treenode " << qq << " split id " << f << " edge " << i << endl;
		}

		boost::dynamic_bitset<> treesplit1 = pg_splitlist[i];
		boost::dynamic_bitset<> treesplit2 = pg_splitlist[i];
		treesplit2.flip();

		for( size_t j=0; j < refpg.E; j++ ){
			// logical AND
			boost::dynamic_bitset<> refsplit1 = ref_splitlist[j];
			boost::dynamic_bitset<> refsplit2 = ref_splitlist[j];
			refsplit2.flip();
//			cout << "rs1.count() " << refsplit1.count() << "\trs2.count() " << refsplit2.count() << endl;
			normalize_split( refsplit1, species_to_gene_map, pg_splitlist[i].size() );
			normalize_split( refsplit2, species_to_gene_map, pg_splitlist[i].size() );
//			cout << "normalized rs1.count() " << refsplit1.count() << "\trs2.count() " << refsplit2.count() << endl;
			
			boost::dynamic_bitset<> and11 = treesplit1 & refsplit1;
			boost::dynamic_bitset<> and21 = treesplit2 & refsplit1;
			boost::dynamic_bitset<> and12 = treesplit1 & refsplit2;
			boost::dynamic_bitset<> and22 = treesplit2 & refsplit2;
			double a11score = (double)and11.count() / (double)treesplit1.count();
			double a22score = (double)and22.count() / (double)treesplit2.count();
			double a1122score = (a11score + a22score) / 2.0;
			double a12score = (double)and12.count() / (double)treesplit1.count();
			double a21score = (double)and21.count() / (double)treesplit2.count();
			double a1212score = (a12score + a21score) / 2.0;
			a1212score = pow( a1212score, 100.0 );
			a1122score = pow( a1122score, 100.0 );
			maxscores.push_back( max(a1122score, a1212score));
			scoresum += maxscores.back();
			bestscore = max(maxscores.back(), bestscore);
		}
		// count the number of nodes with the max score. if it is more than a threshold, ignore this node since it is too hard to reconcile
		int place_count = 0;
		for(size_t j=0; j<maxscores.size(); j++){
			if(maxscores[j] < bestscore)
				continue;
			place_count++;
		}
		if(place_count < placement_limit ){
			for(size_t j=0; j<maxscores.size(); j++){
				if(maxscores[j] < bestscore)
					continue;
				string refnodename = reftree[ refpg.edge_array[j].first ].name;
	//			cout << "gene tree edge " << i << " linking " << other_map[tree[pg.edge_array[i].first].name] << " best reftree edge " << refnodename << endl; 
	//			cout << "found edge " << pg.edge_array[i].first << "\n";
				mapout << edgenum_map[pg.edge_array[i].first] << "\t" << refnodename << endl;
			}
		}else{
			cerr << "Mapping too ambiguous for node " << i << endl;
		}
//		if(pg_splitlist[i].count() == 1)
//			return;
	}
}


