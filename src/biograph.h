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
#ifndef CLUSTER_BIO_GRAPH_H
#define CLUSTER_BIO_GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "graph.h"

namespace clu
{

class BioEdge
{
	public:
		BioEdge();
		BioEdge(const int &, const int &, const double &);
		~BioEdge();
		
		int x, y;
		double score;
};

class BioGraph
{
	public:
		BioGraph();
		~BioGraph();
		
		void clear();
		std::istream &readTAB2(std::istream &, const std::string &);

		template <typename T>
		void convert(Graph &g, const T &conv) const
		{
			//Convert edges to neighbour arrays.
			g.clear();
			
			if (vertices.empty() || edges.empty()) return;
			
			g.nrVertices = static_cast<int>(vertices.size());
			g.nrEdges = static_cast<int>(edges.size());
			g.neighbourRanges.assign(g.nrVertices, make_int2(0, 0));
			g.vertexWeights.assign(g.nrVertices, 0);
			g.neighbours.assign(2*g.nrEdges, make_int2(0, 0));
#ifndef LEAN
			g.coordinates.assign(g.nrVertices, make_float2(0.0, 0.0));
#endif
			
			for (std::vector<BioEdge>::const_iterator i = edges.begin(); i != edges.end(); ++i)
			{
				g.neighbourRanges[i->x].y++;
				g.neighbourRanges[i->y].y++;
			}
			
			for (int i = 1; i < g.nrVertices; ++i) g.neighbourRanges[i].x += g.neighbourRanges[i - 1].x + g.neighbourRanges[i - 1].y;
			
			assert(g.neighbourRanges[g.nrVertices - 1].x + g.neighbourRanges[g.nrVertices - 1].y == 2*g.nrEdges);
			
			for (int i = 0; i < g.nrVertices; ++i) g.neighbourRanges[i].y = g.neighbourRanges[i].x;
			
			//Score edges with the desired metric.
			for (std::vector<BioEdge>::const_iterator i = edges.begin(); i != edges.end(); ++i)
			{
				const int score = conv(i->score);
				
				g.neighbours[g.neighbourRanges[i->x].y++] = make_int2(i->y, score);
				g.neighbours[g.neighbourRanges[i->y].y++] = make_int2(i->x, score);
			}
			
#ifndef NDEBUG
			//Verify ranges.
			assert(g.neighbourRanges[0].x == 0);
			
			for (int i = 0; i < g.nrVertices - 1; ++i) assert(g.neighbourRanges[i].y == g.neighbourRanges[i + 1].x);
			
			assert(g.neighbourRanges[g.nrVertices - 1].y == static_cast<int>(g.neighbours.size()));
#endif
			
			g.setClusterWeights();
		};
		
		std::vector<std::string> vertices;
		std::vector<BioEdge> edges;
};

}

#endif

