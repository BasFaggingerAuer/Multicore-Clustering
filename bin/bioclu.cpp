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
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#define SAVE_SLIDES

using namespace std;
using namespace clu;
using namespace boost;

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
		void drawGraphMatrixPermutation(const Graph &, const std::vector<int> &);
		void drawGraphMatrixClustering(const Graph &, const std::vector<int> &);
		void drawGraphCoordinates(const Graph &);
		void drawGraphClustering(const Graph &, const std::vector<int> &);
		
	private:
		void addSlide();
		
		SDL_Surface * const screen;
		const int size;
		const int pitch;
		const unsigned int clearColour;
		vector<int> heatMap;
		vector<pixel> hueMap;
		
		int curSlide;
		vector<SDL_Surface *> slides;
		
		SDL_Rect origRect, permRect, hueRect;
};

//Different available scoring functions for the graph edges.
class ScoreFunction
{
	public:
		ScoreFunction(const string &_name) : name(_name), t(0.5) {};
		virtual ~ScoreFunction() {};
		
		void parameter(const double &_t)
		{
			assert(_t >= 0.0 && _t <= 1.0);
			t = _t;
		};
		
		virtual int operator () (const double &) const = 0;
		
	public:
		const string name;
		
	protected:
		double t;
};

class ScoreThreshold : public ScoreFunction
{
	public:
		ScoreThreshold() : ScoreFunction("threshold") {};
		~ScoreThreshold() {};
		
		int operator () (const double &score) const
		{
			return (fabs(score) >= t ? 1 : 0);
		};
};

class ScorePower : public ScoreFunction
{
	public:
		ScorePower() : ScoreFunction("power") {};
		~ScorePower() {};
		
		int operator () (const double &score) const
		{
			return 1 + static_cast<int>(floor(256.0*pow(fabs(score), 4.0*t)));
		};
};

class ScoreGaussian : public ScoreFunction
{
	public:
		ScoreGaussian() : ScoreFunction("Gaussian") {};
		~ScoreGaussian() {};
		
		int operator () (const double &score) const
		{
			const double sigma = 0.55 - 0.5*t;
			return static_cast<int>(floor(256.0*exp((min(fabs(score), 1.0) - 1.0)/(sigma*sigma))));
		};
};

class ScoreScale : public ScoreFunction
{
	public:
		ScoreScale() : ScoreFunction("scale") {};
		~ScoreScale() {};
		
		int operator () (const double &score) const
		{
			return static_cast<int>(256.0*fabs(score));
		};
};

int main(int argc, char **argv)
{
	int drawSize = 1024;
	int scoreMode = 0;
	
	string experiment = "Negative Genetic";
	double minCorr = 0.0;
	string geneOntologyFileName = "";
	string fileName = "";
	string gnuplotFileName = "";

	//Parse command line options.
	try
	{
		program_options::options_description desc("Options");
		
		desc.add_options()
		("help,h", "show this help message")
		("input-file", program_options::value<string>(), "set graph input file")
		("experiment,e", program_options::value<string>(), "set BioGRID experiment (e.g. Positive Genetic)")
		("corr,c", program_options::value<double>(), "set minimum correlation")
		("ontology,o", program_options::value<string>(), "read Gene Ontology data")
		("score,s", program_options::value<int>(), "set score mode (0 = none, 1 = treshhold, 2 = power, 3 = Gaussian)");
		
		program_options::positional_options_description pos;
		
		pos.add("input-file", -1);
		
		program_options::variables_map varMap;
		program_options::store(program_options::command_line_parser(argc, argv).options(desc).positional(pos).run(), varMap);
		program_options::notify(varMap);
		
		if (varMap.count("help"))
		{
			cerr << desc << endl;
			return 1;
		}
		
		if (varMap.count("input-file")) fileName = varMap["input-file"].as<string>();
		if (varMap.count("ontology")) geneOntologyFileName = varMap["ontology"].as<string>();
		if (varMap.count("experiment")) experiment = varMap["experiment"].as<string>();
		if (varMap.count("corr")) minCorr = varMap["corr"].as<double>();
		if (varMap.count("score")) scoreMode = varMap["score"].as<int>();
		
		if (fileName == "")
		{
			cerr << "You have to specify an input file!" << endl;
			throw std::exception();
		}
	}
	catch (std::exception &e)
	{
		cerr << "Invalid command line arguments!" << endl;
		return -1;
	}

	cerr << "Creating clustering of '" << fileName << "'..." << endl;

	//Read graph from disk.
	BioGraph bioGraph;

	try
	{
		if (fileName.find(".tab2.txt") != string::npos)
		{
			//We are reading a BioGRID TAB2 file.
			ifstream file(fileName.c_str());

			bioGraph.readTAB2(file, experiment);
			file.close();
		}
		else if (fileName.find(".txt.gz") != string::npos)
		{
			//We are reading a compressed RAW file.
			ifstream file(fileName.c_str(), ios_base::binary);
			iostreams::filtering_istream inStream;

			inStream.push(iostreams::gzip_decompressor());
			inStream.push(file);
			bioGraph.readRAW(inStream, minCorr);
			file.close();
		}
		else
		{
			//We are reading an uncompressed RAW file.
			ifstream file(fileName.c_str());

			bioGraph.readRAW(file, minCorr);
			file.close();
		}
	}
	catch (std::exception &e)
	{
		cerr << "An std::exception occured when reading '" << fileName << "' from disk!" << endl;
		return -1;
	}
	
	//Read ontology if desired.
	GeneOntology ontology;
	
	if (!geneOntologyFileName.empty())
	{
		try
		{
			ifstream file(geneOntologyFileName.c_str(), ios_base::binary);
			iostreams::filtering_istream inStream;

			inStream.push(iostreams::gzip_decompressor());
			inStream.push(file);
			ontology.readGene2go(inStream);
			file.close();
		}
		catch (std::exception &e)
		{
			cerr << "An std::exception occured when reading '" << fileName << "' from disk!" << endl;
			return -1;
		}
	}
	
	//Start SDL to display the results.
	if (SDL_Init(SDL_INIT_VIDEO) != 0)
	{
		cerr << "Unable to initialise SDL: " << SDL_GetError() << endl;
		return -1;
	}
	
	SDL_Surface *screen = SDL_SetVideoMode(32 + 2*drawSize, drawSize, 24, SDL_SWSURFACE);
	
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
	ClusterTBB *cluster = new ClusterTBB();
	
	//Create score function.
	ScoreFunction *score;
	
	if (scoreMode == 1) score = new ScoreThreshold();
	else if (scoreMode == 2) score = new ScorePower();
	else if (scoreMode == 3) score = new ScoreGaussian();
	else score = new ScoreScale();
	
	cerr << "Using " << score->name << " score mode (" << scoreMode << "), experiment '" << experiment << "', and minimum correlation " << minCorr << " ..." << endl;
	
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
				
				score->parameter(static_cast<double>(quality)/9.0);
				bioGraph.convert(graph, *score);
				
				if (!graph.empty())
				{
					vector<int> permute(graph.nrVertices);
					vector<int> component = cluster->clusterPermute(permute, graph, 1.0, 0);
		
					drawer.drawGraphMatrix(graph);
					drawer.drawGraphMatrixPermutation(graph, permute);
					
					cerr << "Generated clustering with modularity " << Cluster::modularity(graph, component) << "." << endl;
				
					if (!ontology.genes.empty()) ontology.clusterOntology(cout, bioGraph, component);
				}
			}
			catch (std::exception &e)
			{
				cerr << "Unable to generate clustering!" << endl;
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
			
			for (double sigma = 0.0; sigma <= 1.0; sigma += 0.01)
			{
				score->parameter(sigma);
				bioGraph.convert(graph, *score);
				
				cerr << "\r" << sigma << "      ";
				
				//There should not be any isolated vertices.
				bool isolated = false;
				
				for (vector<int2>::const_iterator i = graph.neighbourRanges.begin(); i != graph.neighbourRanges.end() && !isolated; ++i) if (i->x >= i->y) isolated = true;
				
				if (!isolated && !graph.empty())
				{
					vector<int> permute(graph.nrVertices);
					const vector<int> clustering = cluster->clusterPermute(permute, graph, 1.0, 0);
					const double modularity = Cluster::modularity(graph, clustering);
					
					if (modularity > bestModularity)
					{
						bestModularity = modularity;
						bestSigma = sigma;
						bestClustering = clustering;
						
						cerr << "Found modularity " << bestModularity << " clustering at sigma " << bestSigma << "." << endl;
						drawer.drawGraphMatrix(graph);
						drawer.drawGraphMatrixPermutation(graph, permute);
						SDL_Flip(screen);
						SDL_Delay(10);
					}
				}
			}
			
			cerr << endl;
			
			//Export related genes.
			cout << "Found modularity " << bestModularity << " clustering at sigma " << bestSigma << " for " << score->name << " scoring method." << endl;
			if (!ontology.genes.empty()) ontology.clusterOntology(cout, bioGraph, bestClustering);
			
			cerr << "Done." << endl;
		}
		
		SDL_Delay(25);
	}
	
	//Free data.
	delete score;
	delete cluster;

	SDL_Quit();

	return 0;
}

DrawerSDL::DrawerSDL(SDL_Surface *_screen, const int &_size) :
	Drawer(),
	screen(_screen),
	size(_size),
	pitch(screen->w),
	clearColour(0x00000000),
	heatMap(size*size, 0),
	curSlide(0)
{
	assert(screen);
	
	origRect.x = 0; origRect.y = 0; origRect.w = size; origRect.h = size;
	permRect.x = size; permRect.y = 0; permRect.w = size; permRect.h = size;
	hueRect.x = 2*size; hueRect.y = 0; hueRect.w = 32; hueRect.h = size;
	
	hueMap.assign(1024, pixel(0, 0, 0));
	
	for (size_t i = 0; i < hueMap.size(); ++i)
	{
		double s = static_cast<double>(i)/static_cast<double>(hueMap.size());
		
		s = 1.0 - pow(1.0 - s, 4.0);
		
		hueMap[i] = pixel(static_cast<unsigned char>(127.0*(1.0 + cos(1.5*M_PI*s + 2.0*M_PI/3.0))),
				static_cast<unsigned char>(127.0*(1.0 + cos(1.5*M_PI*s + 4.0*M_PI/3.0))),
				static_cast<unsigned char>(127.0*(1.0 + cos(1.5*M_PI*s + 0.0*M_PI/3.0))));
	}
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
	char outFile[] = "clu0000.bmp";
	
	sprintf(outFile, "clu%04d.bmp", (int)slides.size());
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
	//Clear heatmap.
	heatMap.assign(size*size, 0);

	for (int i = 0; i < graph.nrVertices; ++i)
	{
		const int y = ((long)size*(long)i)/(long)graph.nrVertices;
		int *row = &heatMap[y*size];
		const int2 r = graph.neighbourRanges[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			const int2 n = graph.neighbours[j];
			
			row[(((long)size*(long)n.x)/(long)graph.nrVertices)] += n.y;
		}
	}
	
	//Determine minimum and maximum.
	const long maxHeat = 1 + *max_element(heatMap.begin(), heatMap.end());
	
	//Clear part of the screen.
	SDL_FillRect(screen, &origRect, clearColour);
	
	SDL_LockSurface(screen);
	
	pixel *destRow = static_cast<pixel *>(screen->pixels);
	const int *src = &heatMap[0];
	
	destRow = &destRow[origRect.y*pitch + origRect.x];
	
	//Draw heatmap.
	for (int i = 0; i < size; ++i)
	{
		pixel *dest = destRow;
		
		for (int j = 0; j < size; ++j)
		{
			*dest++ = hueMap[(static_cast<long>(*src++)*hueMap.size())/maxHeat];
		}
		
		destRow += pitch;
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

void DrawerSDL::drawGraphMatrixPermutation(const Graph &graph, const vector<int> &pi)
{
	assert((int)pi.size() == graph.nrVertices);
	
	//Draw permuted graph matrix.
	vector<int> piInv(graph.nrVertices);
	
	for (int i = 0; i < graph.nrVertices; ++i) piInv[pi[i]] = i;
	
	//Clear heatmap.
	heatMap.assign(size*size, 0);

	for (int i = 0; i < graph.nrVertices; ++i)
	{
		const int y = ((long)size*(long)piInv[i])/(long)graph.nrVertices;
		int *row = &heatMap[y*size];
		const int2 r = graph.neighbourRanges[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			const int2 n = graph.neighbours[j];
			
			row[(((long)size*(long)piInv[n.x])/(long)graph.nrVertices)] += n.y;
		}
	}
	
	//Determine minimum and maximum.
	const long minHeat = *min_element(heatMap.begin(), heatMap.end()), maxHeat = 1 + *max_element(heatMap.begin(), heatMap.end());
	
	//Clear part of the screen.
	SDL_FillRect(screen, &permRect, clearColour);
	
	SDL_LockSurface(screen);
	
	pixel *destRow = static_cast<pixel *>(screen->pixels);
	const int *src = &heatMap[0];
	
	destRow = &destRow[permRect.y*pitch + permRect.x];
	
	//Draw heatmap.
	for (int i = 0; i < size; ++i)
	{
		pixel *dest = destRow;
		
		for (int j = 0; j < size; ++j)
		{
			*dest++ = hueMap[(static_cast<long>(*src++ - minHeat)*hueMap.size())/(maxHeat - minHeat)];
			//*dest++ = hueMap[(static_cast<long>(j)*hueMap.size())/size];
		}
		
		destRow += pitch;
	}
	
	//Draw hues.
	SDL_FillRect(screen, &hueRect, clearColour);
	
	destRow = static_cast<pixel *>(screen->pixels);
	destRow = &destRow[hueRect.y*pitch + hueRect.x];
	
	for (int i = 0; i < size; ++i)
	{
		pixel *dest = destRow;
		const pixel hue = hueMap[(static_cast<long>(size - i - 1)*hueMap.size())/size];
		
		for (int j = 0; j < 32; ++j)
		{
			*dest++ = hue;
		}
		
		destRow += pitch;
	}
	
	SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, permRect.x, permRect.y, permRect.w, permRect.h);
	
	addSlide();
}

void DrawerSDL::drawGraphMatrixClustering(const Graph &graph, const vector<int> &cmp)
{
	assert((int)cmp.size() == graph.nrVertices);
	
	//Draw clustered graph matrix.
	SortByPart sorter(cmp);
	vector<int> pi(graph.nrVertices);
	vector<int> piInv(graph.nrVertices);
	
	for (int i = 0; i < graph.nrVertices; ++i) pi[i] = i;
	
	sort(pi.begin(), pi.end(), sorter);
	
	for (int i = 0; i < graph.nrVertices; ++i) piInv[pi[i]] = i;
	
	//Clear heatmap.
	heatMap.assign(size*size, 0);

	for (int i = 0; i < graph.nrVertices; ++i)
	{
		const int y = ((long)size*(long)piInv[i])/(long)graph.nrVertices;
		int *row = &heatMap[y*size];
		const int2 r = graph.neighbourRanges[i];
		
		for (int j = r.x; j < r.y; ++j)
		{
			const int2 n = graph.neighbours[j];
			
			row[(((long)size*(long)piInv[n.x])/(long)graph.nrVertices)] += n.y;
		}
	}
	
	//Determine minimum and maximum.
	const long maxHeat = 1 + *max_element(heatMap.begin(), heatMap.end());
	
	//Clear part of the screen.
	SDL_FillRect(screen, &permRect, clearColour);
	
	SDL_LockSurface(screen);
	
	pixel *destRow = static_cast<pixel *>(screen->pixels);
	const int *src = &heatMap[0];
	
	destRow = &destRow[permRect.y*pitch + permRect.x];
	
	//Draw heatmap.
	for (int i = 0; i < size; ++i)
	{
		pixel *dest = destRow;
		
		for (int j = 0; j < size; ++j)
		{
			*dest++ = hueMap[(static_cast<long>(*src++)*hueMap.size())/maxHeat];
			//*dest++ = hueMap[(static_cast<long>(j)*hueMap.size())/size];
		}
		
		destRow += pitch;
	}
	
	//Draw hues.
	SDL_FillRect(screen, &hueRect, clearColour);
	
	destRow = static_cast<pixel *>(screen->pixels);
	destRow = &destRow[hueRect.y*pitch + hueRect.x];
	
	for (int i = 0; i < size; ++i)
	{
		pixel *dest = destRow;
		const pixel hue = hueMap[(static_cast<long>(size - i - 1)*hueMap.size())/size];
		
		for (int j = 0; j < 32; ++j)
		{
			*dest++ = hue;
		}
		
		destRow += pitch;
	}
	
	SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, permRect.x, permRect.y, permRect.w, permRect.h);
	
	addSlide();
}

