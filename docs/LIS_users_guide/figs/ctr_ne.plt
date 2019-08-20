set terminal postscript color 
set output 'ne.ps'

set xrange [-93.25:-92.5]
set yrange [31.5:32.25]

set xtics rotate
set xtics -100, .25, -90
set ytics 25, .25, 35

set mxtics
set mytics

#set xlabel "create\nsome\nspace"

unset key
unset border

set grid xtics ytics mxtics mytics

set label 1 "North-east corner" at -92.7,32.05
set label 2 "North-east  center" at -92.835,31.93

set arrow from -92.7,32.05 to -92.75,32.0 head filled
set arrow from -92.835,31.93 to -92.875,31.875 head filled

plot "center_ne_box.dat" using 1:2 with lines lt -1 lw 2,\
     "-" using 1:2 with points lt 7
         -92.75 32.0
         -92.875 31.875
         e
