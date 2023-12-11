#!/usr/bin/gnuplot 

set term png enhanced size 600,300 

set size ratio -1
set contour
unset surface
set view map
unset key

set output "img/dzeta-1000.png"
set cbr [-200:350]
set cntrparam levels increment -200,10,350
sp "files/Q_-1000.txt" u 1:2:4:4 w l lt -1 palette  t '' 

set output "img/dzeta-4000.png"
set cbr [-800:1200]
set cntrparam levels increment -800,20,1200
sp "files/Q_-4000.txt" u 1:2:4:4 w l lt -1 palette  t '' 