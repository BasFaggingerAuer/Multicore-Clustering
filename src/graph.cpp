/*
Copyright 2011, Bas Fagginger Auer.

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
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <cassert>

#include "graph.h"

using namespace clu;
using namespace std;

Graph::Graph() :
	nrVertices(0),
	nrEdges(0),
	neighbourRanges(),
	vertexWeights(),
	neighbours(),
	coordinates(),
	Omega(0)
{

}

Graph::~Graph()
{

}

bool Graph::empty() const
{
	if (nrVertices == 0 || nrEdges == 0 || neighbours.empty() || neighbourRanges.empty()) return true;

	return false;
}

void Graph::clear()
{
	nrVertices = 0;
	nrEdges = 0;
	neighbourRanges.clear();
	vertexWeights.clear();
	neighbours.clear();
	coordinates.clear();
}

istream &Graph::readMETIS(istream &in)
{
	//Reads a METIS graph file from disk.
	int flags = 0;
	int nrVertexWeights = 0;

	clear();

	//First read header.
	while (in.good())
	{
		string line;

		getline(in, line);

		//Skip all comments and empty lines.
		while (!line.empty() && isspace(line[0])) line.erase(0, 1);

		if (line.length() < 2) continue;
		if (line[0] == '%') continue;

		//Read header information.
		istringstream sline(line);

		if (!(sline >> nrVertices >> nrEdges))
		{
			cerr << "Unable to read METIS header!" << endl;
			throw exception();
		}

		//Read optional arguments.
		nrVertexWeights = 0;

		sline >> flags >> nrVertexWeights;

		break;
	}

	//Determine flags.
	const bool weightedEdges = ((flags % 10) != 0), weightedVertices = (((flags/10) % 10) != 0);
	//bool singleEdges = (((flags/100) % 10) != 0);

	if (nrVertexWeights != 0 && !weightedVertices)
	{
		cerr << "Invalid vertex weight specification!" << endl;
		throw exception();
	}

	if (weightedVertices && nrVertexWeights == 0) nrVertexWeights = 1;

	//Resize graph arrays.
	neighbourRanges.assign(nrVertices, make_int2(0, 0));
	neighbours.reserve(nrEdges);
	coordinates.assign(nrVertices, make_float2(0.0, 0.0));

	if (nrVertexWeights > 1) cerr << "Warning: we only record the first vertex weight instead of all " << nrVertexWeights << "!" << endl;
	if (weightedVertices || nrVertexWeights > 0) cerr << "Warning: for clustering all vertex weights are discarded!" << endl;

	vertexWeights.assign(nrVertices, 0);

	//Read vertex neighbours from stream.
	bool warnedSelf = false, warnedDuplicates = false, warnedMultiple = false;

	for (int i = 0; i < nrVertices; ++i)
	{
		if (!in.good())
		{
			cerr << "Read error!" << endl;
			throw exception();
		}

		string line;

		getline(in, line);
		
		istringstream sline(line);

		//Read vertex weights.
		if (weightedVertices)
		{
			int vwgt = 1;
			
			//Discard vertex weights when clustering.
			for (int j = 0; j < nrVertexWeights; ++j) sline >> vwgt;
		}

		//Read vertex neighbours.
		set<int> tmpNeighbours;
		int var;

		neighbourRanges[i].x = neighbours.size();

		while (sline >> var)
		{
			//Verify neighbour index.
			if (var < 1 || var > nrVertices)
			{
				cerr << "Invalid neighbour index (" << var << ", " << nrVertices << ")!" << endl;
				throw exception();
			}

			//Read edge weight if desired.
			int wgt = 1;

			if (weightedEdges) sline >> wgt;

			//Add new neighbour, self-edges are permitted, but not stored in the neighbour lists.
			//Verify that the neighbour has not occurred earlier.
			if (var - 1 != i)
			{
				if (tmpNeighbours.insert(var).second)
				{
					neighbours.push_back(make_int2(var - 1, wgt));
				}
				else if (!warnedDuplicates)
				{
					warnedDuplicates = true;
					cerr << "Warning: duplicate neighbours!" << endl;
				}
			}
			else
			{
				//If self-edges occur, store their weight in the vertex weight.
				if (vertexWeights[i] != 0 && !warnedMultiple)
				{
					warnedMultiple = true;
					cerr << "Warning: multiple self-edges!" << endl;
				}
				
				vertexWeights[i] += wgt;
				
				if (!warnedSelf)
				{
					warnedSelf = true;
					cerr << "This graph has self-edges." << endl;
				}
			}
		}
		
		neighbourRanges[i].y = neighbours.size();

#ifndef NDEBUG
		if (neighbourRanges[i].y > nrVertices + neighbourRanges[i].x)
		{
			cerr << "Too many vertex neighbours!" << endl;
			throw exception();
		}
#endif
	}
	
	setClusterWeights();

#ifndef NDEBUG
	cerr << "Read a METIS graph with " << nrVertices << " vertices and " << nrEdges << " edges." << endl;
#endif

	return in;
}

istream &Graph::readTAB2(istream &in, const string &experiment)
{
	//Reads a BioGRID TAB2 file from disk.
	map<int, int> identifiers;
	vector<int3> edges;

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
		if (identifiers.find(u) == identifiers.end()) identifiers[u] = nrVertices++;
		if (identifiers.find(v) == identifiers.end()) identifiers[v] = nrVertices++;
		
		edges.push_back(make_int3(identifiers[u], identifiers[v], static_cast<int>(ceil(fabs(1000.0f*atof(words[18].c_str()))))));
	}
	
	if (nrVertices <= 0 || edges.empty() || identifiers.empty())
	{
		cerr << "Empty graph!" << endl;
		throw exception();
	}
	
	//Convert edges to neighbour arrays.
	neighbourRanges.assign(nrVertices, make_int2(0, 0));
	vertexWeights.assign(nrVertices, 0);
	neighbours.assign(2*edges.size(), make_int2(0, 0));
	coordinates.assign(nrVertices, make_float2(0.0, 0.0));
	nrEdges = edges.size();
	
	for (vector<int3>::const_iterator i = edges.begin(); i != edges.end(); ++i)
	{
		neighbourRanges[i->x].y++;
		neighbourRanges[i->y].y++;
	}
	
	for (int i = 1; i < nrVertices; ++i) neighbourRanges[i].x += neighbourRanges[i - 1].x + neighbourRanges[i - 1].y;
	
	assert(neighbourRanges[nrVertices - 1].x + neighbourRanges[nrVertices - 1].y == 2*static_cast<int>(edges.size()));
	
	for (int i = 0; i < nrVertices; ++i) neighbourRanges[i].y = neighbourRanges[i].x;
	
	for (vector<int3>::const_iterator i = edges.begin(); i != edges.end(); ++i)
	{
		neighbours[neighbourRanges[i->x].y++] = make_int2(i->y, i->z);
		neighbours[neighbourRanges[i->y].y++] = make_int2(i->x, i->z);
	}
	
#ifndef NDEBUG
	//Verify ranges.
	assert(neighbourRanges[0].x == 0);
	
	for (int i = 0; i < nrVertices - 1; ++i) assert(neighbourRanges[i].y == neighbourRanges[i + 1].x);
	
	assert(neighbourRanges[nrVertices - 1].y == static_cast<int>(neighbours.size()));
#endif
	
	setClusterWeights();

#ifndef NDEBUG
	cerr << "Read a BioGRID TAB2 graph with " << nrVertices << " vertices and " << nrEdges << " edges." << endl;
#endif

	return in;
}

istream &Graph::readCoordinates(std::istream &in)
{
	if (nrVertices <= 0)
	{
		cerr << "Cannot read coordinates for an empty graph!" << endl;
		throw exception();
	}
	
	vector<double> xCoords, yCoords;
	
	xCoords.reserve(nrVertices);
	yCoords.reserve(nrVertices);
	
	for (int i = 0; i < nrVertices; ++i)
	{
		string line;
		
		if (!in.good())
		{
			cerr << "Faulty coordinate data stream!" << endl;
			throw exception();
		}
		
		getline(in, line);
		
		istringstream sline(line);
		
		double x = 0, y = 0;
		
		sline >> x >> y;
		
		xCoords.push_back(x);
		yCoords.push_back(y);
	}
	
	assert((int)xCoords.size() == nrVertices && (int)yCoords.size() == nrVertices);
	assert((int)coordinates.size() == nrVertices);
	
	//Normalise coordinates.
	const double xMin = *min_element(xCoords.begin(), xCoords.end());
	const double yMin = *min_element(yCoords.begin(), yCoords.end());
	const double xMax = *max_element(xCoords.begin(), xCoords.end());
	const double yMax = *max_element(yCoords.begin(), yCoords.end());
	
	for (int i = 0; i < nrVertices; ++i)
	{
		coordinates[i] = make_float2((xCoords[i] - xMin)/(xMax - xMin),
					(yCoords[i] - yMin)/(yMax - yMin));
	}
	
#ifndef NDEBUG
	cerr << "Read " << nrVertices << " coordinates." << endl;
#endif

	return in;
}

vector<int> Graph::random_shuffle()
{
	//Randomizes the vertex order of the graph.
	vector<int> permutation(nrVertices);
	
	for (int i = 0; i < nrVertices; ++i) permutation[i] = i;

	//Create a random permutation.
	std::random_shuffle(permutation.begin(), permutation.end());
	
	//Permute neighbour ranges and coordinates.
	if (true)
	{
		vector<int2> tmpRanges(nrVertices);
		
		for (int i = 0; i < nrVertices; ++i) tmpRanges[i] = neighbourRanges[permutation[i]];
		
		neighbourRanges = tmpRanges;
		
		vector<float2> tmpCoords(nrVertices);
		
		for (int i = 0; i < nrVertices; ++i) tmpCoords[i] = coordinates[permutation[i]];
		
		coordinates = tmpCoords;
	}

	//Permute vertex weights.
	if (true)
	{	
		vector<int> tmpWeights(nrVertices);
		
		for (int i = 0; i < nrVertices; ++i) tmpWeights[i] = vertexWeights[permutation[i]];
		
		vertexWeights = tmpWeights;
	}
	
	//Create inverse permutation.
	vector<int> invPermutation(nrVertices);
	
	for (int i = 0; i < nrVertices; ++i) invPermutation[permutation[i]] = i;
	
	//Apply inverse permutation to vertex indices and edges.
	for (vector<int2>::iterator i = neighbours.begin(); i != neighbours.end(); ++i) i->x = invPermutation[i->x];

	return permutation;
}

Graph Graph::createStarGraph(const int &nrPoints)
{
	Graph graph;
	
	graph.nrVertices = nrPoints + 1;
	graph.nrEdges = nrPoints;
	
	//Create edges from centre (vertex nrPoints) to each of the points.
	for (int i = 0; i < nrPoints; ++i)
	{
		graph.neighbourRanges.push_back(make_int2(i, i + 1));
		graph.neighbours.push_back(make_int2(nrPoints, 1));
	}
	
	//Add centre.
	graph.neighbourRanges.push_back(make_int2(nrPoints, 2*nrPoints));
	
	for (int i = 0; i < nrPoints; ++i)
	{
		graph.neighbours.push_back(make_int2(i, 1));
	}
	
	return graph;
}

Graph Graph::createTreeGraph(const int &, const int &)
{
	//TODO: Implement this.
	throw exception();
	
	return Graph();
}

Graph Graph::createGridGraph(const int &w, const int &h)
{
	Graph graph;
	
	graph.nrVertices = w*h;
	
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			//Presently only add internal edges.
			//TODO: Add border vertices.
			if (x > 0 && y > 0 && x < w - 1 && y < h - 1)
			{
				graph.neighbourRanges.push_back(make_int2(graph.neighbours.size(), graph.neighbours.size() + 4));
				
				graph.neighbours.push_back(make_int2((x - 1) + w*y, 1));
				graph.neighbours.push_back(make_int2((x + 1) + w*y, 1));
				graph.neighbours.push_back(make_int2(x + w*(y - 1), 1));
				graph.neighbours.push_back(make_int2(x + w*(y + 1), 1));
			}
		}
	}
	
	graph.nrEdges = graph.neighbours.size()/2;
	
	return graph;
}

Graph Graph::createLineGraph(const int &size)
{
	Graph graph;
	
	assert(size >= 1);

	graph.nrVertices = size;
	graph.nrEdges = size - 1;
	
	if (size > 1)
	{
		for (int x = 0; x < size; ++x)
		{
			if (x == 0)
			{
				graph.neighbourRanges.push_back(make_int2(graph.neighbours.size(), graph.neighbours.size() + 1));
				graph.neighbours.push_back(make_int2(x + 1, 1));
			}
			else if (x == size - 1)
			{
				graph.neighbourRanges.push_back(make_int2(graph.neighbours.size(), graph.neighbours.size() + 1));
				graph.neighbours.push_back(make_int2(x - 1, 1));
			}
			else
			{
				graph.neighbourRanges.push_back(make_int2(graph.neighbours.size(), graph.neighbours.size() + 2));
				graph.neighbours.push_back(make_int2(x - 1, 1));
				graph.neighbours.push_back(make_int2(x + 1, 1));
			}
		}
	}
	
	return graph;
}

void Graph::setClusterWeights()
{
	//Set vertex weights equal to the sum of the weights of the incoming edges.
	Omega = 0;
	
	for (int i = 0; i < nrVertices; ++i)
	{
		const int2 r = neighbourRanges[i];
		//Include self-edges.
		int w = 2*vertexWeights[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			w += neighbours[j].y;
			Omega += neighbours[j].y;
		}
		
		vertexWeights[i] = w;
	}
	
	assert((Omega & 1) == 0);
	Omega /= 2;

#ifndef NDEBUG
	cerr << "Set " << nrVertices << " clustering vertex weights." << endl;
#endif
}

