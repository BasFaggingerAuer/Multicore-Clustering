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
#include <map>
#include <set>
#include <cassert>
#include <tbb/tbb.h>

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

//Calculate all integer edge weights in parallel.
template <typename T>
class BioGraphWeightFor
{
	public:
		BioGraphWeightFor(std::vector<int> &_v, const std::vector<BioEdge> &_w, const T &_conv) : v(_v), w(_w), conv(_conv) {};
		~BioGraphWeightFor() {};
		
		void operator () (const tbb::blocked_range<size_t> &r) const
		{
			for (size_t i = r.begin(); i != r.end(); ++i) v[i] = conv(w[i].score);
		};
		
	private:
		std::vector<int> &v;
		const std::vector<BioEdge> &w;
		const T &conv;
};

class BioGraph
{
	public:
		BioGraph();
		~BioGraph();
		
		void clear();
		std::istream &readTAB2(std::istream &, const std::string &);
		std::istream &readRAW(std::istream &, const double &);

		template <typename T>
		void convert(Graph &g, const T &conv) const
		{
			//Convert edges to neighbour arrays.
			g.clear();
			
			if (vertices.empty() || edges.empty()) return;
			
			g.nrVertices = static_cast<int>(vertices.size());
			g.neighbourRanges.assign(g.nrVertices, make_int2(0, 0));
			g.vertexWeights.assign(g.nrVertices, 0);
#ifndef LEAN
			g.coordinates.assign(g.nrVertices, make_float2(0.0, 0.0));
#endif
			
			//Calculate edge scores in parallel.
			std::vector<int> edgeScores(edges.size(), 0);
			BioGraphWeightFor<T> tmp(edgeScores, edges, conv);
			
			tbb::parallel_for(tbb::blocked_range<size_t>(0, edges.size()), tmp);
			

#ifndef NDEBUG
			std::cerr << "Scores lie between " << *min_element(edgeScores.begin(), edgeScores.end()) << " and " << *max_element(edgeScores.begin(), edgeScores.end()) << "." << std::endl;
#endif
			
			//Count number of edges.
			g.nrEdges = 0;
			
			for (std::vector<BioEdge>::const_iterator i = edges.begin(); i != edges.end(); ++i)
			{
				const int score = edgeScores[i - edges.begin()];
				
				assert(i->x >= 0 && i->x < g.nrVertices && i->y >= 0 && i->y < g.nrVertices);
				
				if (score > 0)
				{
					g.nrEdges++;
					g.neighbourRanges[i->x].y++;
					g.neighbourRanges[i->y].y++;
				}
			}
			
			if (g.nrEdges <= 0)
			{
				std::cerr << "WARNING: Empty graph with these edge weights!" << std::endl;
				g.clear();
				return;
			}
			
			g.neighbours.assign(2*g.nrEdges, make_int2(0, 0));
			
			for (int i = 1; i < g.nrVertices; ++i) g.neighbourRanges[i].x += g.neighbourRanges[i - 1].x + g.neighbourRanges[i - 1].y;
			
			assert(g.neighbourRanges[g.nrVertices - 1].x + g.neighbourRanges[g.nrVertices - 1].y == 2*g.nrEdges);
			
			for (int i = 0; i < g.nrVertices; ++i) g.neighbourRanges[i].y = g.neighbourRanges[i].x;
			
			//Score edges with the desired metric.
			for (std::vector<BioEdge>::const_iterator i = edges.begin(); i != edges.end(); ++i)
			{
				const int score = edgeScores[i - edges.begin()];
				
				if (score > 0)
				{
					g.neighbours[g.neighbourRanges[i->x].y++] = make_int2(i->y, score);
					g.neighbours[g.neighbourRanges[i->y].y++] = make_int2(i->x, score);
				}
			}
			
#ifndef NDEBUG
			//Verify ranges.
			assert(g.neighbourRanges[0].x == 0);
			
			for (int i = 0; i < g.nrVertices - 1; ++i) assert(g.neighbourRanges[i].y == g.neighbourRanges[i + 1].x);
			
			assert(g.neighbourRanges[g.nrVertices - 1].y == static_cast<int>(g.neighbours.size()));
#endif
			
			g.setClusterWeights();
			
#ifndef NDEBUG
			std::cerr << "Created a graph with " << g.nrVertices << " vertices and " << g.nrEdges << "/" << edges.size() << " (" << (100*g.nrEdges)/edges.size() << "%) edges from biological data." << std::endl;
#endif
		};
		
		std::vector<int> vertices;
		std::vector<BioEdge> edges;
};

class GeneOntology
{
	public:
		GeneOntology();
		~GeneOntology();
		
		void clear();
		std::istream &readGene2go(std::istream &);
		void clusterOntology(std::ostream &, const BioGraph &, const std::vector<int> &) const;
		
		std::map<int, std::set<std::string> > genes;
};

}

#endif

