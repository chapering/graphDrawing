
#include "newick_file.h"
#include "newick.h"
#include "tree.h"

#include <fstream>
#include <iostream>
#include <string>
/*
#include <QFile>
#include <QString>
#include <QTextStream>
*/
using namespace BiRC::treelib;
using namespace std;

// FIXME: There is no error handling here!
std::auto_ptr<Tree> BiRC::treelib::parse_newick_file(const char *fname)
{
/*
	QFile curFile( fname );
	if ( curFile.open( QFile::ReadOnly | QFile::Text ))
	{
		QTextStream in( &curFile );
		string str = in.readAll().toStdString();
		curFile.close();

		std::auto_ptr<Tree> tree( parse_newick(str) );
		return tree;
	}
	else	return std::auto_ptr<Tree>(0);
	*/

    string str;
    ifstream curFile(fname, ios::in);

    curFile >> str;

    if(!str.empty())
    {
		/*
        cerr << "reading tree... " << fname <<  endl;
        cerr << "========> read ..... " << str << endl;
		*/

        std::auto_ptr<Tree> tree( parse_newick(str) );
        return tree;
    }
    else    return std::auto_ptr<Tree>(0);

    curFile.close();
}

std::auto_ptr<Tree> BiRC::treelib::parse_newick_file(const std::string &fname)
{
    return parse_newick_file(fname.c_str());
}
