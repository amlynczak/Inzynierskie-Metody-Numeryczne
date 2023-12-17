# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

reset
set term gif size 800,300 animate delay 10
set output "D0.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("files/0zad5_it=%i.txt",i)
  splot file u 1:2:3 with pm3d title sprintf("t=%i",i)
} 

reset
set term gif size 800,300 animate delay 10
set output "D01.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("files/0.1zad5_it=%i.txt",i)
  splot file u 1:2:3 w pm3d title sprintf("t=%i",i)
} 