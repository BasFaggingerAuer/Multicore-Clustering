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
#include <string>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
	//Small program to extract relevant coexpression data and normalise the correlation coefficients.
	if (argc != 2)
	{
		cerr << "Usage: zcat coexpression.txt.gz | " << argv[0] << " 0.4 | gzip -c > data.txt.gz." << endl;
		return 1;
	}
	
	const double minScore = atof(argv[1]);
	long nrLines = 0, nextNrLines = 0;
	
	while (cin.good())
	{
		string line;

		getline(cin, line);

		//Skip all comments and empty lines.
		while (!line.empty() && isspace(line[0])) line.erase(0, 1);

		if (line.length() < 2) continue;
		if (line[0] == '#') continue;
		
		//Convert all separators to spaces.
		for (string::iterator i = line.begin(); i != line.end(); ++i) if (*i == '\t' || *i == ',' || *i == ';' || *i == ':' || *i == '|') *i = ' ';
		
		//Extract all words.
		vector<string> words;
		istringstream sline(line);
		string word;
		
		while (sline >> word) words.push_back(word);
		
		if (words.size() != 3)
		{
			cerr << "Invalid number of entries (" << words.size() << ") for '" << line << "'!" << endl;
			return -1;
		}
		
		const double score = atof(words[2].c_str()) - 1.0;
		
		if (score >= minScore)
		{
			//Output data.
			cout << atol(words[0].c_str()) << " " << atol(words[1].c_str()) << " " << score << endl;
		}
		
		if (++nrLines > nextNrLines)
		{
			cerr << "\r" << nextNrLines;
			nextNrLines += 100000;
			cerr.flush();
		}
	}
	
	cerr << endl << "Done." << endl;
	
	return 0;
}

