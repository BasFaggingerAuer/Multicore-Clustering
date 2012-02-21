#!/bin/sh
# 
# GPU clustering DIMACS benchmark script.
#
# Usage: ./dimacs.sh /data/benchmarkoutput /data/dimacsroot graphs-for-clustering.txt
#
tbbid=5
cudadevice=1
cudaid=6
prefix=$1
dimacsdata=$2

mkdir -p $prefix
mkdir -p $prefix/tbb
mkdir -p $prefix/cuda
rm -f $prefix/graphs.txt

cat $3 | xargs -n 1 build/bin/graphstat -p $dimacsdata >> $prefix/graphs.txt

rm -f $prefix/tbb/*.ptn
rm -f $prefix/tbb/*.eval
rm -f $prefix/cuda/*.ptn
rm -f $prefix/cuda/*.eval

cat $prefix/graphs.txt | sort -n | cut -f 2 | xargs -n 1 build/bin/dimacs -d -1 -i $tbbid -p $prefix/tbb
cat $prefix/graphs.txt | sort -n | cut -f 2 | xargs -n 1 build/bin/dimacs -d $cudadevice -i $cudaid -p $prefix/cuda
