set term postscript eps enhanced color
set key under
set border 3 front 
set grid xtics nomxtics ytics nomytics
set tics out
set xtics nomirror rotate by -90 scale 0
set ytics border nomirror scale 0
#set xtics border nomirror 
 
# nb       style
# 0        nop
# 1        +
# 2        *
# 3        +*   
# 4        empty  square
# 5        filled square 
# 6        empty  circle
# 7        filled circle
# 8        empty  triangle
# 9        filled triangle

set style line 1  lt 1  lc  rgb "#123456"  lw 6 pt 0 ps 3

set output "data/source.eps"
set title  "Ricker wavelet"
set xlabel "time (iterations)"
set ylabel "amplitude"
set xrange[0:100] 
set nologscale  xy  
plot "data/source.dat"  using 1:2 title "" w lp ls 1
