set terminal postscript color 
set output 'sw.ps'

set xrange [-98.0:-97.25]
set yrange [27.5:28.25]

set xtics rotate
set xtics -100, .25, -90
set ytics 25, .25, 35

set mxtics
set mytics

#set xlabel "create\nsome\nspace"

unset key
unset border

set grid xtics ytics mxtics mytics

set label 1 "South-west corner" at -97.72,27.68
set label 2 "South-west  center" at -97.59,27.82

set arrow from -97.72,27.68 to -97.75,27.75 head filled
set arrow from -97.59,27.82 to -97.625,27.875 head filled

plot "center_sw_box.dat" using 1:2 with lines lt -1 lw 2,\
     "-" using 1:2 with points lt 7
         -97.75 27.75
         -97.625 27.875
         e
