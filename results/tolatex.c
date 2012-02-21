/*
Results to LaTeX table converters, by Bas Fagginger Auer.
*/
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>

#define MAX_LINE 8192

void doubleToLaTeX(char *d, const double v)
{
	char *dp = d;
	
	sprintf(d, "%#.3f", v);
	
	while (*dp != '\0')
	{
		if (*dp == '.') *dp = '&';
		dp++;
	}
	
	/*
	char tmp[1024];
	char *cp = tmp;
	
	if (v < 0.0)
	{
		*dp++ = '-';
		*dp++ = '\0';
		return;
	}
	
	sprintf(tmp, "%.1e", v);
	
	while (*cp != '\0' && *cp != 'e') *dp++ = *cp++;
	
	cp++;
	*dp++ = ' ';
	*dp++ = '\\';
	*dp++ = 'c';
	*dp++ = 'd';
	*dp++ = 'o';
	*dp++ = 't';
	*dp++ = ' ';
	*dp++ = '1';
	*dp++ = '0';
	*dp++ = '^';
	*dp++ = '{';
	*dp = 'z';
	
	while (*cp != '\0')
	{
		if (*cp != '0' && *cp != '+') *dp++ = *cp;
		cp++;
	}
	
	if (*dp == 'z') *dp++ = '0';
	
	*dp++ = '}';
	*dp++ = '\0';
	*/
}

int main(int argc, char **argv)
{
	char line[MAX_LINE];
	FILE *in1, *in2;
	FILE *out = stdout;
	
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s cuda.txt tbb.txt\n", argv[0]);
		return EXIT_FAILURE;
	}
	
	in1 = fopen(argv[1], "r");
	in2 = fopen(argv[2], "r");
	
	if (in1 == NULL || in2 == NULL)
	{
		fprintf(stderr, "Unable to open either '%s' or '%s'!\n", argv[1], argv[2]);
		return EXIT_FAILURE;
	}
	
	setlocale(LC_NUMERIC, "en_US.utf-8");
	
	while (!feof(in1) || !feof(in2))
	{
		char name[MAX_LINE];
		long nrVertices, nrEdges;
		
		char time1Str[MAX_LINE], time2Str[MAX_LINE];
		double mod1 = -1.0, mod1Dev = 0.0, time1 = -1.0, time1Dev = -1.0;
		double mod2 = -1.0, mod2Dev = 0.0, time2 = -1.0, time2Dev = -1.0;
		
		if (fgets(line, MAX_LINE, in1) != NULL)
		{
			sscanf(line, "%s %ld %ld %lf %lf %lf %lf", name, &nrVertices, &nrEdges, &mod1, &mod1Dev, &time1, &time1Dev);
		}
		
		if (fgets(line, MAX_LINE, in2) != NULL)
		{
			sscanf(line, "%s %ld %ld %lf %lf %lf %lf", name, &nrVertices, &nrEdges, &mod2, &mod2Dev, &time2, &time2Dev);
		}
		
		doubleToLaTeX(time1Str, time1);
		doubleToLaTeX(time2Str, time2);
		
		fprintf(out, "\\texttt{%s} & %'ld & %'ld & %.2f & %s & %.2f & %s \\\\\n", name, nrVertices, nrEdges, mod1, time1Str, mod2, time2Str);
	}
	
	fclose(in1);
	fclose(in2);
	
	return EXIT_SUCCESS;
}

