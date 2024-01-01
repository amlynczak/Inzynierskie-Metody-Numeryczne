reset
set term gif size 300,300 animate delay 20
set output "img/T.gif"
set view map # widok z gory
set size ratio -1
set cbr [0:]

file = "files/T_100.txt"
splot file u 1:2:3 with pm3d 

file = "files/T_200.txt"
splot file u 1:2:3 with pm3d 

file = "files/T_500.txt"
splot file u 1:2:3 with pm3d 

file = "files/T_1000.txt"
splot file u 1:2:3 with pm3d 

file = "files/T_2000.txt"
splot file u 1:2:3 with pm3d 


reset
set term gif size 300,300 animate delay 20
set output "img/NT.gif"
set view map # widok z gory
set size ratio -1
set cbr [0:]

file = "files/NT_100.txt"
splot file u 1:2:3 with pm3d 

file = "files/NT_200.txt"
splot file u 1:2:3 with pm3d 

file = "files/NT_500.txt"
splot file u 1:2:3 with pm3d 

file = "files/NT_1000.txt"
splot file u 1:2:3 with pm3d 

file = "files/NT_2000.txt"
splot file u 1:2:3 with pm3d 