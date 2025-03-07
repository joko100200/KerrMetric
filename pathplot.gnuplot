set ticslevel 0
set xrange[-4000:4000]
set yrange[-4000:4000]
set zrange[-4000:4000]
set size ratio -1
set parametric
set urange [0:2*pi]
set vrange [0:pi]
set isosamples 20,20
splot 200*cos(u)*sin(v), 200*sin(u)*sin(v), 200*cos(v),\
 'xyz.out' with lines
pause -1
