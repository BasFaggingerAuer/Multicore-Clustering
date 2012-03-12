/*
Copyright 2012, Bas Fagginger Auer.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <exception>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>

#include "biograph.h"

using namespace clu;
using namespace std;

BioEdge::BioEdge() :
	x(0),
	y(0),
	score(0.0)
{

}

BioEdge::BioEdge(const int &_x, const int &_y, const double &_score) :
	x(_x),
	y(_y),
	score(_score)
{

}

BioEdge::~BioEdge()
{

}

BioGraph::BioGraph()
{

}

BioGraph::~BioGraph()
{

}

void BioGraph::clear()
{
	vertices.clear();
	edges.clear();
}

istream &BioGraph::readTAB2(istream &in, const string &experiment)
{
	//Reads a BioGRID TAB2 file from disk.
	int nrVertices = 0;
	map<int, int> identifiers;

	clear();

	while (in.good())
	{
		string line;

		getline(in, line);

		//Skip all comments and empty lines.
		while (!line.empty() && isspace(line[0])) line.erase(0, 1);

		if (line.length() < 2) continue;
		if (line[0] == '#') continue;
		
		//Find all words within this line, separated by tabs.
		vector<string> words;
		istringstream sline(line);
		string word;
		
		while (getline(sline, word, '\t')) words.push_back(word);
		
		if (words.size() != 24)
		{
			cerr << "Invalid number of entries!" << endl;
			throw exception();
		}
		
		//Only include desired interactions.
		if (words[11] != experiment) continue;
		
		const int u = atoi(words[3].c_str()), v = atoi(words[4].c_str());
		
		if (u == 0 || v == 0)
		{
			cerr << "Invalid identifiers!" << endl;
			throw exception();
		}
		
		//Have we already seen these identifiers?
		if (identifiers.find(u) == identifiers.end())
		{
			identifiers[u] = nrVertices++;
			vertices.push_back(words[7]);
		}
		
		if (identifiers.find(v) == identifiers.end())
		{
			identifiers[v] = nrVertices++;
			vertices.push_back(words[8]);
		}
		
		edges.push_back(BioEdge(identifiers[u], identifiers[v], atof(words[18].c_str())));
		//cout << vertices[identifiers[u]] << " -- " << vertices[identifiers[v]] << ": " << atof(words[18].c_str()) << endl;
	}
	
	assert(nrVertices == static_cast<int>(vertices.size()));
	
	if (nrVertices <= 0 || edges.empty() || identifiers.empty())
	{
		cerr << "Empty graph!" << endl;
		throw exception();
	}

#ifndef NDEBUG
	cerr << "Read a BioGRID TAB2 graph with " << vertices.size() << " vertices and " << edges.size() << " edges." << endl;
#endif
	
	return in;
}

