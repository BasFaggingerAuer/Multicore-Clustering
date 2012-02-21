#!/bin/sh
# 
# GPU clustering benchmarking script.
# First run stats.sh to generate a graphs file.
#
# Usage: ./benchmark.sh
#
tbbfile="results/tbb.txt";
cudafile="results/cuda.txt";
device=1
nravg=16

rm -f $tbbfile
cat results/graphs.txt | sort -n | cut -f 2 | xargs -n 1 build/bin/clu -a $nravg -m 1 >> $tbbfile
rm -f $cudafile
cat results/graphs.txt | sort -n | cut -f 2 | xargs -n 1 build/bin/clu -a $nravg -d $device -m 2 >> $cudafile
