#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __PhyloTree_h__
#define __PhyloTree_h__

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stack>

//typedef unsigned int node_id_t;
typedef size_t node_id_t;
class TreeNode 
{
public:
	TreeNode() : distance(0) {};
	std::string name;	/**< node name */
	double distance;	/**< distance to parent */
	std::vector< node_id_t > parents;	/**< if parents.size() == 0 this is a root node */
	std::vector< node_id_t > children;	/**< if children.size() == 0 this is a leaf node */
};

template< class T >
class PhyloTree 
{
public:
	PhyloTree();
	PhyloTree( const PhyloTree<T>& pt );
	PhyloTree<T>& operator=( const PhyloTree<T>& pt );
	double weight;	/**< Overall tree weight */
	node_id_t root;	/**< root of the tree */
	std::vector< T > nodes;	/**< nodes of the tree */
	void clear();
	/**
	 * Reads a tree in Newick format.  WARNING:  only reads rooted trees correctly
	 */
	void readTree( std::istream& tree_file );
	/**
	 * Writes a tree in Newick format
	 */
	void writeTree( std::ostream& os ) const;
	/**
	 * Determines the height of the tree along the path from the root to the left-most leaf node
	 */
	double getHeight() const;
	/**
	 * Determines the height of the tree along the path from nodeI to its left-most descendant leaf node
	 */
	double getHeight( node_id_t nodeI ) const;

	T& operator[]( const unsigned i ){ return nodes[i]; }
	const T& operator[]( const unsigned i ) const{ return nodes[i]; }
	size_t size() const{ return nodes.size(); }
	void push_back( T& t ){ nodes.push_back(t); }
	T& back() { return nodes.back(); }
	const T& back() const{ return nodes.back(); }
	void resize( const unsigned s ){ nodes.resize(s); }


	void swap( PhyloTree<T>& other )
	{
		std::swap( weight, other.weight );
		std::swap( root, other.root );
		nodes.swap( other.nodes );
	}
protected:
};


template< class T >
PhyloTree<T>::PhyloTree()
{
	weight = 0;
	root = 0;
}

template< class T >
PhyloTree<T>::PhyloTree( const PhyloTree<T>& pt ) :
nodes( pt.nodes ),
weight( pt.weight ),
root( pt.root )
{}

template< class T >
PhyloTree<T>& PhyloTree<T>::operator=( const PhyloTree<T>& pt )
{
	nodes = pt.nodes;
	weight = pt.weight;
	root = pt.root;
	return *this;
}

template< class T >
void PhyloTree<T>::clear()
{
	nodes.clear();
	weight = 0;
	root = 0;
}


/**
 *  readTree version 2.0: read in a phylogenetic tree in the Newick file format.
 *
 */
template< class T >
void PhyloTree<T>::readTree( std::istream& tree_file )
{
	std::string line;
	clear();
	if( !std::getline( tree_file, line ) )
		return;
	// look for either a ; or a matched number of parenthesis, if
	// not found then read another line
	while(true){
		int paren_count = 0;
		for( size_t charI = 0; charI < line.size(); charI++ )
		{
			if( line[charI] == '(' )
				paren_count++;
			if( line[charI] == ')' )
				paren_count--;
		}
		if( paren_count == 0 )
			break;
		if( paren_count != 0 ){
			std::string another_line;
			if( !std::getline( tree_file, another_line ) )
				return;
			line += another_line;
		}
	}

	std::stringstream line_str( line );

	// look for a weight
	std::string::size_type open_bracket_pos = line.find( "[" );
	std::string::size_type bracket_pos = line.find( "]" );
	if( open_bracket_pos != std::string::npos && bracket_pos != std::string::npos && 
		open_bracket_pos < bracket_pos && bracket_pos < line.find( "(" ) ){
		// read in a weight
		getline( line_str, line, '[' );
		getline( line_str, line, ']' );
		std::stringstream weight_str( line );
		weight_str >> weight;
	}
	
	// ready to begin parsing the tree data.
	std::string tree_line;
	std::getline( line_str, tree_line, ';' );
	size_t read_state = 0;	/**< read_state of 0 indicates nothing has been parsed yet */
	size_t section_start = 0;
	std::stack< node_id_t > node_stack;
	std::stringstream blen_str;
	T new_node;
	new_node.distance = 0;	// default the distance to 0
	bool already_read_name = false;
	bool blen_found = false;
	for( size_t charI = 0; charI < tree_line.size(); charI++ ){
		switch( tree_line[ charI ] ){
			// if this is an open parens then simply create a new
			// parent node and push it on the parent stack
			case '(':
				if( node_stack.size() > 0 ){
					new_node.parents.clear();
					new_node.parents.push_back( node_stack.top() );
					(*this)[ node_stack.top() ].children.push_back( (node_id_t)(*this).size() );
				}
				node_stack.push( (node_id_t)(*this).size() );
				nodes.push_back( new_node );
				read_state = 1;
				section_start = charI + 1;
				break;
			case ')':
				if( blen_found )
				{
					// read off a branch length
					blen_str.clear();
					blen_str.str( tree_line.substr( section_start, charI - section_start ) );
					blen_str >> (*this)[ node_stack.top() ].distance;
				}else{
					// read off a name, if possible
					if( read_state == 1 ){
						new_node.parents.clear();
						new_node.parents.push_back( node_stack.top() );
						(*this)[ node_stack.top() ].children.push_back( (node_id_t)(*this).size() );
						node_stack.push( (node_id_t)(*this).size() );
						nodes.push_back( new_node );
						read_state = 2;	// pop this node after reading its branch length
					}
					(*this)[ node_stack.top() ].name = tree_line.substr( section_start, charI - section_start );
				}
				if( read_state == 2 )
					node_stack.pop();
				section_start = charI + 1;
				blen_found = false;

				// pop off the top of the node stack
				read_state = 2;
				break;
			case ',':
				if( blen_found ){
					// read off a branch length
					blen_str.clear();
					blen_str.str( tree_line.substr( section_start, charI - section_start ) );
					blen_str >> (*this)[ node_stack.top() ].distance;
				}else{
					// read off a name, if possible
					if( read_state == 1 ){
						new_node.parents.clear();
						new_node.parents.push_back( node_stack.top() );
						(*this)[ node_stack.top() ].children.push_back( (node_id_t)(*this).size() );
						node_stack.push( (node_id_t)(*this).size() );
						nodes.push_back( new_node );
						read_state = 2;	// pop this node after reading its name
					}
					(*this)[ node_stack.top() ].name = tree_line.substr( section_start, charI - section_start );
				}
				if( read_state == 2 )
					node_stack.pop();
				section_start = charI + 1;
				read_state = 1;	// indicates that we'll be creating a new node when we hit :
				blen_found = false;
				break;
			case ':':
				// read off a name, if possible
				if( read_state == 1 ){
					new_node.parents.clear();
					new_node.parents.push_back( node_stack.top() );
					(*this)[ node_stack.top() ].children.push_back( (node_id_t)(*this).size() );
					node_stack.push( (node_id_t)(*this).size() );
					nodes.push_back( new_node );
					read_state = 2;	// pop this node after reading its branch length
				}
				(*this)[ node_stack.top() ].name = tree_line.substr( section_start, charI - section_start );
				section_start = charI + 1;
				blen_found = true;
				break;
			default:
				break;
		}
	}

}


template< class T >
void PhyloTree<T>::writeTree( std::ostream& os ) const{
	std::stack< node_id_t > node_stack;
	std::stack< size_t > child_stack;
	node_stack.push( root );
	child_stack.push( 0 );
	bool write_branch_lengths = false;
	for( size_t nodeI = 0; nodeI < this->size(); nodeI++ )
	{
		if( (*this)[nodeI].distance != 0 )
		{
			write_branch_lengths = true;
			break;
		}
	}

	if( (*this).weight != 0 )
		os << "[" << weight << "]";
	os << "(";

	while( node_stack.size() > 0 ) {
		if( (*this)[ node_stack.top() ].children.size() != 0 ){
			// this is a parent node
			// if we have scanned all its children then pop it
			if( child_stack.top() == (*this)[ node_stack.top() ].children.size() ){
				os << ")";
				if( node_stack.size() > 1 && write_branch_lengths )
					os << ":" << (*this)[ node_stack.top() ].distance;
				node_stack.pop();
				child_stack.pop();
				continue;
			}
			// try to recurse to its children
			// if the child is a parent as well spit out a paren
			node_id_t child = (*this)[ node_stack.top() ].children[ child_stack.top() ];
			node_stack.push( child );
			child_stack.top()++;
			// print a comma to separate multiple children
			if( child_stack.top() > 1 )
				os << ",";
			if( (*this)[ child ].children.size() > 0 ){
				child_stack.push( 0 );
				os << "(";
			}
			continue;
		}
		
		// this is a leaf node
		os << (*this)[ node_stack.top() ].name;
		if( write_branch_lengths )
			os << ":" << (*this)[ node_stack.top() ].distance;
		
		// pop the child
		node_stack.pop();
	}
	os << ";" << std::endl;
}


template< class T >
double PhyloTree<T>::getHeight() const
{
	return getHeight( root );
}

template< class T >
double PhyloTree<T>::getHeight( node_id_t nodeI ) const
{
	if( (*this)[ nodeI ].children.size() == 0 )
		return (*this)[ nodeI ].distance;
	return (*this)[ nodeI ].distance + getHeight( (*this)[ nodeI ].children[ 0 ] );
}


/** determine which nodes are descendants of a given node */
template< class TreeType >
void getDescendants( TreeType& alignment_tree, node_id_t node, std::vector< node_id_t >& descendants )
{
	// do a depth first search
	std::stack< node_id_t > node_stack;
	node_stack.push( node );
	descendants.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() > 0 )
		{
			node_stack.push(alignment_tree[cur_node].children[0]);
			node_stack.push(alignment_tree[cur_node].children[1]);
		}
		descendants.push_back(cur_node);
	}
}


/** determine which nodes are leaf nodes below a given node */
template< class TreeType >
void getLeaves( TreeType& tree, node_id_t node, std::vector< node_id_t >& leaves )
{
	// do a depth first search
	std::stack< node_id_t > node_stack;
	node_stack.push( node );
	leaves.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( tree[cur_node].children.size() > 0 )
		{
			node_stack.push(tree[cur_node].children[0]);
			node_stack.push(tree[cur_node].children[1]);
		}else
			leaves.push_back(cur_node);
	}
}

namespace std {

template< class T > inline
void swap( PhyloTree<T>& a, PhyloTree<T>& b )
{
	a.swap(b);
}

template<> inline void swap( PhyloTree<TreeNode>& a, PhyloTree<TreeNode>& b){ a.swap(b); }
}

#endif // __PhyloTree_h__


