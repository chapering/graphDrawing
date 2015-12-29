
#include "tree.h"
#include <ostream>
#include <algorithm>
#include <cassert>

int BiRC::treelib::Tree::add_node( const std::string label )
{
    int new_node = parent_.size();

	std::vector<int> tmp;
	child_.push_back( tmp );
    parent_.push_back(-1);
    length_to_parent_.push_back(0.0);
    label_.push_back(label);
    
    return new_node;
}

int BiRC::treelib::Tree::add_node( const std::string label, std::vector<int> tChild, std::vector<double> cLength )
{
	int new_node = parent_.size();

	std::vector<int> tmp;
	child_.push_back( tmp );
	parent_.push_back(-1);
	length_to_parent_.push_back(0.0);
	label_.push_back(label);

	for ( int i = tChild.size()-1; i > -1 ; i -- )
	{
		child_[new_node].push_back( tChild[i] );
		parent_[ tChild[i] ] = new_node;
		length_to_parent_[ tChild[i] ] = cLength[i];
	}

	return new_node;
}

void BiRC::treelib::Tree::dfs_print(std::ostream &os, int node) const
{
    assert(node >= 0);
    
    if (is_leaf(node))
	{
	    // leaf
	    os << label_[node];
	}
    else if (is_inner(node))
	{
	    // inner node
	    os << '(';
	    dfs_print(os, child_[node][0]);

		for ( int i = 1; i < child_[node].size(); i ++ )
		{
			os << ',';
			dfs_print(os, child_[node][i]);
		}

	    os << ')';
	}
    else
	assert(!"Inconsistent tree!");

    os << ':' << length_to_parent_[node];

}

void BiRC::treelib::Tree::print(std::ostream &os) const
{
    dfs_print(os, root_);
    os << ';';
}

const std::string BiRC::treelib::Tree::getMaxlenLabel() const
{
	int len = -1, nodeIndex = -1, tt;
	for ( int i = 0; i < label_.size(); i ++ )
	{
		tt = label_[i].length();
		if ( tt > len )
		{
			len = tt;
			nodeIndex = i;
		}
	}
	if ( nodeIndex == -1)
	{
		return "";
	}
	else
	{
		return label_[nodeIndex];
	}
}

int BiRC::treelib::Tree::leafNum()
{
	int num = 0;
	for ( int i = 0; i < label_.size(); i ++ )
	{
		if ( is_leaf( i ) )
		{
			num ++;
		}
	}
	return num;
}
