#include "newick.h"
#include "tree.h"
using namespace BiRC::treelib;

#include <boost/spirit/core.hpp>
#include <stack>
#include <iterator>
#include <istream>
#include <cassert>

namespace {
	struct newick_grammar : public boost::spirit::grammar<newick_grammar>
	{
		mutable std::auto_ptr<Tree> tree;
		mutable std::stack<int> node_stack;
		mutable std::stack<double> branch_len_stack;

		struct handle_branch_start
		{
			newick_grammar &self;
			handle_branch_start(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr, Itr) const
			{
				self.node_stack.push( -1 );
			}
		};

		struct handle_leaf
		{
			newick_grammar &self;
			handle_leaf(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr beg, Itr end) const
			{
				std::string label(beg,end);
				self.node_stack.push(self.tree->add_node(label));
			}
		};

		struct handle_inner_node
		{
			newick_grammar &self;
			handle_inner_node(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr, Itr) const
			{
				std::vector<int> leaves;
				std::vector<double> length;

				int tleaf;
				double tlen;
				
				tleaf = self.node_stack.top();
				self.node_stack.pop();
				
				while ( tleaf != -1 )
				{
					tlen = self.branch_len_stack.top();
					self.branch_len_stack.pop();

					leaves.push_back( tleaf );
					length.push_back( tlen );

					if ( self.node_stack.empty() )
					{
						break;
					}

					tleaf = self.node_stack.top();
					self.node_stack.pop();
				}
				
				self.node_stack.push(self.tree->add_node("", leaves, length ) );
			}
		};

		struct handle_weight
		{
			newick_grammar &self;
			handle_weight(newick_grammar &self) : self(self) {}

			void operator()(double d) const
			{
				self.branch_len_stack.push( d );
			}
		};

		struct handle_weight_default
		{
			newick_grammar &self;
			handle_weight_default(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr, Itr) const
			{
				self.branch_len_stack.push( 0.0 );
			}
		};

		struct handle_tree_start
		{
			newick_grammar &self;
			handle_tree_start(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr, Itr) const
			{
				self.tree->set_root(self.node_stack.top());
			}
		};

		struct handle_tree_end
		{
			newick_grammar &self;
			handle_tree_end(newick_grammar &self) : self(self) {}

			template <typename Itr>
			void operator()(Itr, Itr) const
			{
				std::vector<int> leaves;
				std::vector<double> length;

				int tleaf;
				double tlen;

				tleaf = self.node_stack.top();
				self.node_stack.pop();

				while ( tleaf != -1 )
				{
					tlen = self.branch_len_stack.top();
					self.branch_len_stack.pop();

					leaves.push_back( tleaf );
					length.push_back( tlen );

					if ( self.node_stack.empty() )
					{
						break;
					}
					tleaf = self.node_stack.top();
					self.node_stack.pop();
				}

				self.tree->set_root( self.tree->add_node( "", leaves, length ) );
			}
		};

		handle_leaf leaf_handler;
		handle_branch_start branch_start_handler;
		handle_inner_node inner_node_handler;
		handle_weight weight_handler;
		handle_weight_default weight_default_handler;
		handle_tree_start tree_start_handler;
		handle_tree_end tree_end_handler;

		newick_grammar()
			: tree(new Tree),
			leaf_handler(*this),
			branch_start_handler(*this),
			inner_node_handler(*this),
			weight_handler(*this),
			weight_default_handler(*this),
			tree_start_handler(*this),
			tree_end_handler(*this)
		{
		}

		template <typename scanner_t>
		struct definition
		{
			definition(newick_grammar const &self)
			{
				using namespace boost::spirit;
				
				leaf_r = (alnum_p >> *(~ch_p(':')))[self.leaf_handler];

				inner_node_r =
					('(' >> epsilon_p[self.branch_start_handler]) >> edge_r >> *(',' >> edge_r) >> ')' 
					>> epsilon_p[self.inner_node_handler]
				;

				edge_r =
					leaf_r >> ((':' >> (real_p[self.weight_handler]))
					| epsilon_p[self.weight_default_handler])
					|
					inner_node_r >> ((':' >> (real_p[self.weight_handler]))
					| epsilon_p[self.weight_default_handler])
					;

				// The split in start and end is because we need to
				// follow different rules depending on the number of
				// children of the root in the newick format, and we
				// cannot distinguish between 2 and 3 children from the
				// start of the parsing.
				tree_start_r =
					leaf_r >> ((':' >> (real_p[self.weight_handler]))
					| epsilon_p[self.weight_default_handler]) >> ';' >> epsilon_p[self.tree_start_handler]
				|
					('(' >> epsilon_p[self.branch_start_handler]) >> edge_r >> ',' >> edge_r >> tree_end_r
					;
				tree_end_r =
					epsilon_p >> *( ',' >> edge_r ) >>  ')' >> ((':' >> (real_p[self.weight_handler]))
					| epsilon_p[self.weight_default_handler]) >> ';' 
					>> epsilon_p[self.tree_end_handler]
				;
			}

			boost::spirit::rule<scanner_t> leaf_r, inner_node_r, edge_r;
			boost::spirit::rule<scanner_t> tree_start_r, tree_end_r;
			const boost::spirit::rule<scanner_t> &
				start() const { return tree_start_r; }
		};   
	};
}

namespace {
	template <typename Itr>
	std::auto_ptr<Tree>
		parse(Itr beg, Itr end)
	{
		newick_grammar grammar;
		boost::spirit::parse_info<Itr> info =
			boost::spirit::parse(beg, end, grammar, boost::spirit::space_p);

		if (info.full) return grammar.tree;
		else           return std::auto_ptr<Tree>(0);
	}
}


std::auto_ptr<Tree>
BiRC::treelib::parse_newick(const std::string &str)
{
	return ::parse(str.begin(), str.end());
}


std::auto_ptr<Tree>
BiRC::treelib::parse_newick(const char *beg, const char *end)
{
	return ::parse(beg, end);
}
