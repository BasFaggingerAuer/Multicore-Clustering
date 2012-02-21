set autoscale
unset label
unset grid
set log x
set format x "10^{%L}"
set log y
set format y "10^{%L}"
set xtic nomirror auto
set ytic mirror auto
set grid ytics noxtics
set title "Number of clusters"
set xlabel "Number of graph edges"
set ylabel "Number of clusters (-)"
set key top left
plot "cuda.txt" using 3:8 title 'CUDA' with points, \
"tbb.txt" using 3:8 title 'TBB' with points
