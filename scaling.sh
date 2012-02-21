#!/bin/sh
#
# Script which measures performance scaling as the number of threads increases.
#
# Usage: ./scaling.sh /dimacs/data/random
#
nravg=4;

rm -f results/*.stxt

for f in `find $1 -type f -printf '%s\t%p\n' | sort -n | cut -f 2 | grep .graph.bz2`
do
	fshort=`echo $f | sed 's/.*\///g' | sed 's/\.bz2//g'`;
	echo $f $fshort;
	build/bin/clu -a $nravg -m 3 $f > results/$fshort.stxt;
done

