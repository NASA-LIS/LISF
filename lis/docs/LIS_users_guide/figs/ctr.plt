#set terminal postscript color 
#set output 'ps.ps'

set xrange [-98:-97.0]
set yrange [27.5:28.5]

set xtics rotate
set xtics -100, .25, -90
set ytics 25, .25, 35

set mxtics
set mytics

#set xlabel "create\nsome\nspace"

unset key
unset border

set grid xtics ytics mxtics mytics

plot "center_sw_box.dat" using 1:2 with lines lt -1 lw 2, \
     "center_sw_box.dat" using 1:2 with points
