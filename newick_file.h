#ifndef NEWICK_FILE_H_INCLUDED
#define NEWICK_FILE_H_INCLUDED

#include <string>
#include <memory>

// FIXME: There is no errorhandling in these functions!!!

namespace BiRC { 
	namespace treelib {
		// forward decl.
		class Tree;

		/**
		* Parse a string in Newick format from a file.
		*
		* \param fname	Name of the file to parse.
		* \returns	A tree corresponding to the string, or 0 in case of errors.
		*/
		std::auto_ptr<Tree> parse_newick_file(const char *fname);

		/**
		* Parse a string in Newick format from a file.
		*
		* \param fname	Name of the file to parse.
		* \returns	A tree corresponding to the string, or 0 in case of errors.
		*/
		std::auto_ptr<Tree> parse_newick_file(const std::string &fname);
	}
}

#endif
