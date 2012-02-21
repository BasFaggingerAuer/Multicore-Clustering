set autoscale
set log y
set log x
unset label
unset grid
set xtic nomirror auto
set ytic nomirror auto
set title "Matching time scaling"
set xlabel "Number of CPU threads"
set ylabel "Relative matching time (%)"
set key right top
set xtics 2
set xr [1:8]
set ytics (10,20,30,40,50,60,70,80,90,100)
set grid ytics noxtics
set yr [20:100.1]
plot "rgg_n_2_20_s0.graph.stxt" using 1:(100/$1) title 'ideal scaling' with lines, \
"rgg_n_2_15_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_15\_s0' with linespoints, \
"rgg_n_2_16_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_16\_s0' with linespoints, \
"rgg_n_2_17_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_17\_s0' with linespoints, \
"rgg_n_2_18_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_18\_s0' with linespoints, \
"rgg_n_2_19_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_19\_s0' with linespoints, \
"rgg_n_2_20_s0.graph.stxt" using 1:(a=$3,0/0) every ::0::0 notitle, '' using 1:(100*$3/a) title 'rgg\_n\_2\_20\_s0' with linespoints
