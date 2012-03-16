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
#include <algorithm>
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
			cerr << "Invalid number of entries (" << words.size() << ")!" << endl;
			throw exception();
		}
		
		//Only include desired interactions.
		if (words[11] != experiment) continue;
		
		const int u = atoi(words[1].c_str()), v = atoi(words[2].c_str());
		
		if (u == 0 || v == 0)
		{
			cerr << "Invalid identifiers!" << endl;
			throw exception();
		}
		
		//Have we already seen these identifiers?
		if (identifiers.find(u) == identifiers.end())
		{
			identifiers[u] = nrVertices++;
			vertices.push_back(u);
		}
		
		if (identifiers.find(v) == identifiers.end())
		{
			identifiers[v] = nrVertices++;
			vertices.push_back(v);
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

	cerr << "Read a BioGRID TAB2 graph with " << vertices.size() << " vertices and " << edges.size() << " edges." << endl;
	
	return in;
}

istream &BioGraph::readRAW(istream &in, const double &minScore)
{
	//Reads a raw coexpression file from disk.
	//On each line: Entrez_Gene_ID1  Entrez_Gene_ID2  Correlation.
	int nrVertices = 0;
	map<int, int> identifiers;
#ifndef NDEBUG
	size_t count1 = 0, count2 = 0;
#endif

	clear();

	while (in.good())
	{
		string line;

		getline(in, line);

		//Skip all comments and empty lines.
		while (!line.empty() && isspace(line[0])) line.erase(0, 1);

		if (line.length() < 2) continue;
		if (line[0] == '#') continue;
		
		//Find all words within this line.
		size_t wordStart = 0, wordEnd = 0;
		vector<string> words;
		
		while ((wordEnd = line.find_first_of(" \t,;:|", wordStart + 1)) != string::npos)
		{
			if (wordEnd > wordStart + 1) words.push_back(line.substr(wordStart, wordEnd - wordStart));
			wordStart = wordEnd;
		}
		
		//Get the last word.
		words.push_back(line.substr(wordStart));
		
		if (words.size() != 3)
		{
			cerr << "Invalid number of entries (" << words.size() << ")!" << endl;
			
			throw exception();
		}
		
		const double score = atof(words[2].c_str());
		
		if (score >= minScore)
		{
			const int u = atoi(words[0].c_str()), v = atoi(words[1].c_str());
			
			if (u == 0 || v == 0)
			{
				cerr << "Invalid identifiers!" << endl;
				throw exception();
			}
			
			//Have we already seen these identifiers?
			if (identifiers.find(u) == identifiers.end())
			{
				identifiers[u] = nrVertices++;
				vertices.push_back(u);
			}
			
			if (identifiers.find(v) == identifiers.end())
			{
				identifiers[v] = nrVertices++;
				vertices.push_back(v);
			}
			
			edges.push_back(BioEdge(identifiers[u], identifiers[v], score));
			//cout << vertices[identifiers[u]] << " -- " << vertices[identifiers[v]] << ": " << score << endl;
		}
		
#ifndef NDEBUG
		if (vertices.size() > count1 + 100 || edges.size() > count2 + 1000)
		{
			count1 = vertices.size();
			count2 = edges.size();
			cerr << "\r" << count1 << ", " << count2;
			cerr.flush();
		}
#endif
	}
	
#ifndef NDEBUG
	cerr << endl;
#endif
	
	assert(nrVertices == static_cast<int>(vertices.size()));
	
	if (nrVertices <= 0 || edges.empty() || identifiers.empty())
	{
		cerr << "Empty graph!" << endl;
		throw exception();
	}

	cerr << "Read a RAW coexpression graph with " << vertices.size() << " vertices and " << edges.size() << " edges." << endl;
	
	return in;
}

GeneOntology::GeneOntology()
{

}

GeneOntology::~GeneOntology()
{

}

void GeneOntology::clear()
{
	genes.clear();
}

istream &GeneOntology::readGene2go(istream &in)
{
	//Read the Gene Ontology information from ftp://ftp.ncbi.nih.gov/gene/.
	//On each line: tx_id  Entrez_Gene_ID1  GO_ID  Evidence  Qualifier  GO_term  PubMed_Category.

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
		
		if (words.size() != 8)
		{
			cerr << "Invalid number of entries (" << words.size() << ")!" << endl;
			throw exception();
		}
		
		genes[atoi(words[1].c_str())].insert(words[5]);
	}
	
	cerr << "Read Gene Ontology data for " << genes.size() << " genes." << endl;
	
	return in;
}

void GeneOntology::clusterOntology(ostream &out, const BioGraph &graph, const vector<int> &cluster) const
{
	//Determines statistics about the given clustering.
	assert(graph.vertices.size() == cluster.size());
	
	const int nrClusters = *max_element(cluster.begin(), cluster.end()) + 1;
	
#ifndef NDEBUG
	//Check the total number of clusters.
	if (true)
	{
		set<int> tmp(cluster.begin(), cluster.end());
		
		assert(static_cast<int>(tmp.size()) == nrClusters);
	}
#endif
	
	//Find the total number of vertices in each cluster.
	vector<int> nrClusterVertices(nrClusters, 0);
	
	for (vector<int>::const_iterator i = cluster.begin(); i != cluster.end(); ++i) nrClusterVertices[*i]++;
	
	//Determine the number of singletons.
	int nrSingletons = 0;
	
	for (vector<int>::const_iterator i = nrClusterVertices.begin(); i != nrClusterVertices.end(); ++i)
	{
		assert(*i > 0);
		if (*i == 1) nrSingletons++;
	}
	
	//Determine cluster properties.
	vector<set<string> > clusterProperties(nrClusters, set<string>());
	
	for (int i = 0; i < static_cast<int>(graph.vertices.size()); ++i)
	{
		map<int, set<string> >::const_iterator j = genes.find(graph.vertices[i]);
		
		if (j != genes.end())
		{
			clusterProperties[cluster[i]].insert(j->second.begin(), j->second.end());
		}
		else
		{
			cerr << "Unable to find gene with Entrez ID " << graph.vertices[i] << "!" << endl;
			clusterProperties[cluster[i]].insert("UNKNOWN");
		}
	}
	
	//Write information.
	out << nrClusters << " clusters, of which " << nrSingletons << " (" << (100*nrSingletons)/nrClusters << "%) are singletons:" << endl;
	
	for (int i = 0; i < nrClusters; ++i)
	{
		if (nrClusterVertices[i] > 1)
		{
			out << i << ", " << clusterProperties[i].size() << "/" << nrClusterVertices[i] << " properties: ";
			
			for (set<string>::const_iterator j = clusterProperties[i].begin(); j != clusterProperties[i].end(); ++j) out << *j << ", ";
			
			out << endl;
		}
	}
	
	//Compare to average number of properties of each gene in the graph to determine the quality.
	double clusterAvg = 0.0;
	
	for (int i = 0; i < nrClusters; ++i)
	{
		clusterAvg += static_cast<double>(clusterProperties[i].size())/static_cast<double>(nrClusterVertices[i]);
	}
	
	clusterAvg /= static_cast<double>(nrClusters);
	
	long graphSum = 0;
	
	for (vector<int>::const_iterator i = graph.vertices.begin(); i != graph.vertices.end(); ++i)
	{
		map<int, set<string> >::const_iterator j = genes.find(*i);
		
		if (j != genes.end())
		{
			graphSum += j->second.size();
		}
		else
		{
			graphSum++;
		}
	}
	
	double graphAvg = static_cast<double>(graphSum)/static_cast<double>(graph.vertices.size());
	
	cerr << "Quality: " << graphAvg/clusterAvg << " for " << nrClusters << " clusters and " << nrSingletons << " singletons." << endl;
	out << "Average number of GO properties per gene in this graph: " << graphAvg << ", versus the average of the average number of GO properties per cluster: " << clusterAvg << ", indicates quality " << graphAvg/clusterAvg << "." << endl;
}

