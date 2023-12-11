#!/usr/bin/gnuplot 

set term png enhanced size 600,300 

set size ratio -1
set contour
unset surface
set view map
unset key

#-1000 psi

set cntrparam levels increment -55,0.2,-50
set cbr [-55:-50]
set o "img/psi-1000.png"
set title "Q=-1000, psi"
sp "files/Q_-1000.txt" u 1:2:3:3 w l lt -1 palette  t ''

#-4000 psi

set output "img/psi-4000.png"
set cbr [-218:-202]
set cntrparam levels increment -218,0.5,-202
sp "files/Q_-4000.txt" u 1:2:3:3 w l lt -1 palette  t '' 

#4000 psi

set output "img/psi4000.png"
set cbr [202:218]
set cntrparam levels increment 202,1,218
sp "files/Q_4000.txt" u 1:2:3:3 w l lt -1 palette  t '' 