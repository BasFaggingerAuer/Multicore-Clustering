#!/bin/sh
#
# Script which generates a list of all graphs, together with their number of edges.
#
# Usage: ./stats.sh /dimacs/data/test /dimacs/data/test2
#

rm -f results/graphs.txt

for var in "$@"
do
	find $var -type f -printf '%s\t%p\n' | sort -n | cut -f 2 | grep .graph.bz2 | xargs -n 1 build/bin/graphstat >> results/graphs.txt
done
