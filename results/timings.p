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
set title "Clustering time"
set xlabel "Number of graph edges"
set ylabel "Clustering time (s)"
set key top left
plot "cuda.txt" using 3:6:7 title 'CUDA' with yerrorbars, \
"tbb.txt" using 3:6:7 title 'TBB' with yerrorbars
