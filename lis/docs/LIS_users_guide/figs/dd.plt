set terminal postscript color 
set output 'ps.ps'

set xrange [-98.5:-92.25]
set yrange [27:33]
set xtics rotate
set xtics -100, .25, -90
set ytics 25, .25, 35

#set xlabel "create\nsome\nspace"

unset key
unset border

set grid

# draw the original sub-domain
#set arrow from -97.6,27.9 to -92.9,27.9 nohead
#set arrow from -92.9,27.9 to -92.9,31.9 nohead
#set arrow from -92.9,31.9 to -97.6,31.9 nohead
#set arrow from -97.6,31.9 to -97.6,27.9 nohead

# draw the new snapped-onto domain
#set arrow from -97.75,27.75 to -92.75,27.75 nohead
#set arrow from -92.75,27.75 to -92.75,32.0  nohead
#set arrow from -92.75,32.0  to -97.75,32.0  nohead
#set arrow from -97.75,32.0  to -97.75,27.75 nohead

# draw arrows showing adjustment of the extreme points
set arrow from -97.6,27.9 to -97.75,27.75 head filled
set arrow from -92.9,31.9 to -92.75,32.0  head filled

plot "snap_grid.dat" using 1:2 with lines lt -1 lw 2, \
    "orig_box.dat" using 1:2 with lines lt 3 lw 2
