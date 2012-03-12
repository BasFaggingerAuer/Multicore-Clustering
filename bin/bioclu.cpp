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
#include <iomanip>
#include <fstream>
#include <exception>

#include <graph.h>
#include <biograph.h>
#include <vis.h>
#include <cluster.h>
#include <clustertbb.h>
#include <clustercuda.h>

#include <cmath>

#include <SDL.h>

#include <boost/program_options.hpp>

using namespace clu;
using namespace std;

class pixel
{
	public:
		pixel(const unsigned char &_b, const unsigned char &_g, const unsigned char &_r) :
			r(_r), g(_g), b(_b)
		{

		};
		
		~pixel()
		{

		}
		
		inline pixel operator |= (const pixel &a)
		{
			r |= a.r;
			g |= a.g;
			b |= a.b;
			
			return *this;
		}
		
		unsigned char r, g, b;
};

class DrawerSDL : public Drawer
{
	public:
		DrawerSDL(SDL_Surface *, const int &);
		~DrawerSDL();
		
		void nextSlide();
		void prevSlide();
		
		void drawGraphMatrix(const Graph &);
		void drawGraphMatrixClustering(const Graph &, const std::vector<int> &);
		void drawGraphCoordinates(const Graph &);
		void drawGraphClustering(const Graph &, const std::vector<int> &);
		
	private:
		void addSlide();
		
		inline pixel fromHue(const int &c) const
		{
			const float h = 123.456f*(float)c + 654.123f;
			
			return pixel((unsigned char)(254.0f*min(max(0.0f, 2.0f*cosf(h)), 1.0f)), (unsigned char)(254.0f*min(max(0.0f, 2.0f*cosf(h + 2.0f*M_PI/3.0f)), 1.0f)), (unsigned char)(254.0f*min(max(0.0f, 2.0f*cosf(h + 4.0f*M_PI/3.0f)), 1.0f)));
		};
		
		SDL_Surface * const screen;
		const int size;
		pixel *buffer;
		const int pitch;
		const unsigned int clearColour;
		
		int curSlide;
		vector<SDL_Surface *> slides;
		
		SDL_Rect origRect, permRect;
};

class ScoreGaussian
{
	public:
		ScoreGaussian(const double &sigma) : factor(1.0/(sigma*sigma)) {};
		~ScoreGaussian() {};
		
		int operator () (const double &score) const {return 1 + static_cast<int>(ceil(2048.0*exp(factor*(fabs(score) - 1.0))));};
		
	private:
		const double factor;
};

int main(int argc, char **argv)
{
	int drawSize = 512;
	
	string experiment = "Negative Genetic";
	string fileName = "", shortFileName = "";
	string gnuplotFileName = "";

	//Parse command line options.
	try
	{
		boost::program_options::options_description desc("Options");
		
		desc.add_options()
		("help,h", "show this help message")
		("input-file", boost::program_options::value<string>(), "set graph input file")
		("experiment,e", boost::program_options::value<string>(), "set BioGRID experiment (e.g. Positive Genetic)");
		
		boost::program_options::positional_options_description pos;
		
		pos.add("input-file", -1);
		
		boost::program_options::variables_map varMap;
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(pos).run(), varMap);
		boost::program_options::notify(varMap);
		
		if (varMap.count("help"))
		{
			cerr << desc << endl;
			return 1;
		}
		
		if (varMap.count("input-file")) fileName = varMap["input-file"].as<string>();
		if (varMap.count("experiment")) experiment = varMap["experiment"].as<string>();
		
		if (fileName == "")
		{
			cerr << "You have to specify an input file!" << endl;
			throw exception();
		}
	}
	catch (exception &e)
	{
		cerr << "Invalid command line arguments!" << endl;
		return -1;
	}

	cerr << "Creating clustering of '" << fileName << "'..." << endl;

	//Read graph from disk.
	BioGraph bioGraph;

	try
	{
		//We are reading a BioGRID TAB2 file.
		ifstream file(fileName.c_str());

		bioGraph.readTAB2(file, experiment);
		file.close();
	
		shortFileName = fileName.substr(1 + fileName.find_last_of("/\\"));
	}
	catch (exception &e)
	{
		cerr << "An exception occured when reading " << fileName << " from disk!" << endl;
		return -1;
	}
	
	//Start SDL to display the results.
	if (SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		cerr << "Unable to initialise SDL: " << SDL_GetError() << endl;
		return -1;
	}
	
	SDL_Surface *screen = SDL_SetVideoMode(2*drawSize, drawSize, 24, SDL_SWSURFACE);
	
	if (!screen)
	{
		cerr << "Unable to create window: " << SDL_GetError() << endl;
		return -1;
	}
	
	//Set caption and clear screen.
	SDL_WM_SetCaption("Cluster", "Cluster");
	SDL_FillRect(screen, &screen->clip_rect, 0x00000000);
	SDL_Flip(screen);
	
	//Create graph drawer.
	DrawerSDL drawer(screen, drawSize);
	
	//Create clusterer.
	Cluster *cluster = new ClusterTBB();
	
	//Enter main loop.
	bool running = true;
	
	while (running)
	{
		//Take care of incoming events.
		int quality = -1;
		SDL_Event event;
		
		while (SDL_PollEvent(&event))
		{
			//Handle keypresses.
			if (event.type == SDL_KEYDOWN)
			{
				const SDLKey k = event.key.keysym.sym;
				
				if (k == SDLK_ESCAPE)
				{
					running = false;
				}
				else if (k == SDLK_LEFT)
				{
					drawer.prevSlide();
				}
				else if (k == SDLK_RIGHT)
				{
					drawer.nextSlide();
				}
				else if (k >= '1' && k <= '9')
				{
					quality = k - (int)('1');
				}
				else if (k == '0')
				{
					quality = 10;
				}
			}
			
			if (event.type == SDL_QUIT)
			{
				running = false;
			}
		}
		
		if (quality >= 0 && quality < 10)
		{
			//Start clustering.
			try
			{
				//Convert to ordinary graph.
				Graph graph;
				
				bioGraph.convert(graph, ScoreGaussian(0.05 + static_cast<double>(quality)/16.0));
				
				vector<int> component = cluster->cluster(graph, 1.0, 0);
	
				drawer.drawGraphMatrix(graph);
				drawer.drawGraphMatrixClustering(graph, component);
				
				cout << "Generated clustering with modularity " << Cluster::modularity(graph, component) << "." << endl;
			}
			catch (exception &e)
			{
				cout << "Unable to generate clustering!" << endl;
			}
			
			SDL_Flip(screen);
		}
		else if (quality == 10)
		{
			//Find best clustering.
			Graph graph;
			double bestModularity = -1.0;
			double bestSigma = 0.0;
			vector<int> bestClustering(bioGraph.vertices.size(), 0);
			
			for (double sigma = 0.05; sigma <= 0.50; sigma += 0.001)
			{
				bioGraph.convert(graph, ScoreGaussian(sigma));
				
				const vector<int> clustering = cluster->cluster(graph, 1.0, 0);
				const double modularity = Cluster::modularity(graph, clustering);
				
				if (modularity > bestModularity)
				{
					bestModularity = modularity;
					bestSigma = sigma;
					bestClustering = clustering;
					
					cout << "Found modularity " << bestModularity << " clustering at sigma " << bestSigma << "." << endl;
					drawer.drawGraphMatrixClustering(graph, bestClustering);
					SDL_Flip(screen);
				}
				
				SDL_Delay(25);
			}
			
			//Export related genes.
			const int nrClusters = *max_element(bestClustering.begin(), bestClustering.end()) + 1;
			vector<string> clusterMembers(nrClusters, "");
			
			cout << "For modularity " << bestModularity << " and sigma " << bestSigma << ", we found " << nrClusters << " clusters:" << endl;
			
			for (int i = 0; i < static_cast<int>(bioGraph.vertices.size()); ++i) clusterMembers[bestClustering[i]] += bioGraph.vertices[i] + ", ";
			for (int i = 0; i < nrClusters; ++i) cout << "Cluster " << i + 1 << ":" << endl << clusterMembers[i] << endl;
		}
		
		SDL_Delay(25);
	}
	
	//Free data.
	delete cluster;

	SDL_Quit();

	return 0;
}

DrawerSDL::DrawerSDL(SDL_Surface *_screen, const int &_size) :
	Drawer(),
	screen(_screen),
	size(_size),
	pitch(screen->w),
#ifdef SAVE_SLIDES
	clearColour(0xffffffff),
#else
	clearColour(0x00000000),
#endif
	curSlide(0)
{
	assert(screen);
	
	origRect.x = 0; origRect.y = 0; origRect.w = size; origRect.h = size;
	permRect.x = size; permRect.y = 0; permRect.w = size; permRect.h = size;
}

DrawerSDL::~DrawerSDL()
{
	//Free all slides.
	for (vector<SDL_Surface *>::iterator i = slides.begin(); i != slides.end(); ++i) SDL_FreeSurface(*i);
}

void DrawerSDL::nextSlide()
{
	if (slides.empty()) return;
	
	curSlide = (curSlide + 1) % slides.size();
	SDL_BlitSurface(slides[curSlide], 0, screen, 0);
	SDL_Flip(screen);
}

void DrawerSDL::prevSlide()
{
	if (slides.empty()) return;
	
	curSlide = (curSlide > 0 ? curSlide - 1 : slides.size() - 1);
	SDL_BlitSurface(slides[curSlide], 0, screen, 0);
	SDL_Flip(screen);
}

void DrawerSDL::addSlide()
{
	//Save the current screen as a slide and (optionally) write it to disk.
	//Make copy of screen.
	SDL_Surface *slide = SDL_DisplayFormat(screen);
	
	assert(slide);
	
	SDL_BlitSurface(screen, 0, slide, 0);
	slides.push_back(slide);
	
#ifdef SAVE_SLIDES
	char outFile[] = "map0000.bmp";
	
	sprintf(outFile, "map%04d.bmp", (int)slides.size());
	SDL_SaveBMP(screen, outFile);
#endif
}

void DrawerSDL::drawGraphCoordinates(const Graph &graph)
{
	
}

void DrawerSDL::drawGraphClustering(const Graph &graph, const vector<int> &cmp)
{
	
}

void DrawerSDL::drawGraphMatrix(const Graph &graph)
{
	//Clear part of the screen.
	SDL_FillRect(screen, &origRect, clearColour);
	
	SDL_LockSurface(screen);
	
	buffer = static_cast<pixel *>(screen->pixels);
	
	//Find minimum and maximum weight.
	int minWgt = INT_MAX, maxWgt = INT_MIN;
	
	for (vector<int2>::const_iterator i = graph.neighbours.begin(); i != graph.neighbours.end(); ++i)
	{
		minWgt = min(minWgt, i->y);
		maxWgt = max(maxWgt, i->y);
	}

#ifndef NDEBUG
	cerr << "Edge weights lie between " << minWgt << " and " << maxWgt << "." << endl;
#endif
	
	//No permutation.
	for (int i = 0; i < graph.nrVertices; ++i)
	{
		const int y = ((long)origRect.h*(long)i)/(long)graph.nrVertices;
		pixel *cp = &buffer[(y + origRect.y)*pitch + origRect.x];
		const int2 r = graph.neighbourRanges[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			const int2 n = graph.neighbours[j];
			pixel *dest = &cp[(((long)origRect.w*(long)n.x)/(long)graph.nrVertices)];
			
			*dest = pixel(max(dest->r, static_cast<unsigned char>((255*(n.y - minWgt))/maxWgt)), 0, 255);
		}
	}
	
	SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, origRect.x, origRect.y, origRect.w, origRect.h);
}

class SortByPart
{
	public:
		SortByPart(const vector<int> &_v) : v(_v) {};
		~SortByPart() {};
		
		bool operator () (const int &a, const int &b) const
		{
			return v[a] < v[b];
		};
		
	private:
		const vector<int> &v;
};

void DrawerSDL::drawGraphMatrixClustering(const Graph &graph, const vector<int> &cmp)
{
	assert((int)cmp.size() == graph.nrVertices);
	
	//Draw clustered graph matrix.
	SortByPart sorter(cmp);
	vector<int> pi(graph.nrVertices);
	
	for (int i = 0; i < graph.nrVertices; ++i) pi[i] = i;
	
	sort(pi.begin(), pi.end(), sorter);
	
	//Clear part of the screen.
	SDL_FillRect(screen, &permRect, clearColour);
	
	SDL_LockSurface(screen);
	
	buffer = static_cast<pixel *>(screen->pixels);
	
	//Find minimum and maximum weight.
	int minWgt = INT_MAX, maxWgt = INT_MIN;
	
	for (vector<int2>::const_iterator i = graph.neighbours.begin(); i != graph.neighbours.end(); ++i)
	{
		minWgt = min(minWgt, i->y);
		maxWgt = max(maxWgt, i->y);
	}

	//Draw permuted and coloured graph.
	vector<int> piInv(graph.nrVertices);
	
	for (int i = 0; i < graph.nrVertices; ++i) piInv[pi[i]] = i;
	
	for (int i = 0; i < graph.nrVertices; ++i)
	{
		const int y = ((long)permRect.h*(long)piInv[i])/(long)graph.nrVertices;
		pixel *cp = &buffer[(y + permRect.y)*pitch + permRect.x];
		const int2 r = graph.neighbourRanges[i];
		const int c0 = cmp[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			const int2 n = graph.neighbours[j];
			pixel *dest = &cp[(((long)permRect.w*(long)piInv[n.x])/(long)graph.nrVertices)];
			
			*dest = pixel(max(dest->r, static_cast<unsigned char>((255*(n.y - minWgt))/maxWgt)), (c0 == cmp[n.x] ? 128 : 0), 255);
		}
	}
	
	SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, permRect.x, permRect.y, permRect.w, permRect.h);
	
	addSlide();
}

