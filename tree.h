#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED

//#include <QPoint>
#include <string>
#include <vector>
#include <iosfwd>
#include <cassert>

namespace BiRC { 
	namespace treelib {
		/**
		* \brief Data structure for binary trees.
		*
		* There are essentially no consistency guarantees in this
		* structure except that if the structure is a correct tree before
		* a method call, then it will always be a consistent tree after a
		* method call.  Stronger guarantees about consistency would put
		* too strong restrictions on the tree construction interface.
		*/
		class Tree {

			/// Mapping of a node's index to its parents index.
			std::vector<int> parent_;

			std::vector< std::vector<int> > child_;
			
			/// Mapping of a node's index to its label.
			std::vector<std::string> label_;

			/// Branch lengths (measured from a node to its parent)
			std::vector<double> length_to_parent_;

			int* depth_;

			/// Index of the root of the tree.
			int root_;

			/// Helper function for printing
			void dfs_print(std::ostream &os, int node) const;

		public:

			/// Construct a new empty tree.
			Tree() : root_(-1) {};

			// --- tree statistics... ----------------------------------
			/**
			* \brief Get the number of nodes in the tree
			*
			* Notice, this is the total number of nodes in the tree,
			* <em>not</em> the number of leaves in the the tree.
			*
			* \returns Number of leaves in the tree.
			*/
			int size() { return parent_.size(); }

			// --- predicates ------------------------------------------

			/// Predicate testing if \a node is a leaf.
			bool is_leaf(int node) const
			{ 
				assert(node > -1 && node < (int)child_.size());
				return ( child_[node].size() <= 0 );
			}

			/// Predicate testing if \a node is an inner node.
			bool is_inner(int node) const
			{ 
				assert(node > -1 && node < (int)child_.size());
				return ( child_[node].size() > 0 );
			}

			/// Predicate testing if \a node is the root.
			bool is_root(int node) const
			{ 
				return node == root_;
			}

			// --- node access -----------------------------------------

			int root() const { return root_; }

			/// Get the label associated with \a node.
			const std::string &label(int node) const
			{
				assert( node >=0 && node < (int)label_.size());
				return label_[node];
			}

			std::vector< std::vector<int> > get_edge() const
			{
				return child_;
			}

			std::vector<int> get_child( int node ) const
			{
				assert( node >= 0 && node < child_.size() );
				return child_[node];
			}

			int get_child(int node, int num) const
			{
				assert( node < child_.size() && num < child_[node].size() );
				return child_[node][num];
			}

			/// Get the parent of \a node.
			int parent(int node) const
			{
				assert(node < (int)parent_.size());
				return parent_[node];
			}

			double length_to_parent(int node) const
			{
				assert(node < (int)length_to_parent_.size());
				return length_to_parent_[node];
			}

			int leafNum();
			// --- construction ----------------------------------------

			/**
			* \brief Insert a new node in the tree.
			*
			* Adds a node to the tree.  If children are specified, their
			* parent index is updated (potentially deleting their
			* previous parent).
			*
			*
			* \returns		Index of the new node.
			*/
			int add_node( const std::string label );
			int add_node( const std::string label,  std::vector<int> tChild, std::vector<double> cLength );
			
			const std::string getMaxlenLabel() const;
			
			/**
			* \brief Sets the root index.
			*
			* The \a root index will be considered the root of the tree.
			* There is no guarantee that the sub-tree rooted in \a root
			* will contain the full structure; that must be guaranteed by
			* the tree building code.
			*
			* \param root 		Index to be used as root.
			*
			* \pre \a root must be an index containing an actual node,
			* i.e. one that was returned by add_node().
			*/
			void set_root(int root)
			{
				assert(root < (int)parent_.size());
				root_ = root;
			}


			// --- misc. -----------------------------------------------

			/**
			* \brief Prints the tree structure to a stream.
			*
			* Prints the tree, both in its internal representation and in
			* Newick format.
			*
			* If the tree is rooted in a leaf, that leaf will be printed
			* as a left-most leaf in the Newick representation.
			*/
			void print(std::ostream &os) const;
		};

		inline std::ostream &operator<<(std::ostream &os, const Tree &tree)
		{
			tree.print(os);
			return os;
		}
	}
}

#endif
