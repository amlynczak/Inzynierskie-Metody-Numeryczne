set terminal pngcairo enhanced font 'Verdana,12'
set pm3d map
set size ratio 1.0

set output 'img/T_iter_100.png'
set title 'it=100'
splot 'files/T_100.txt' with pm3d

set output 'img/T_iter_200.png'
set title 'it=200'
splot 'files/T_200.txt' with pm3d

set output 'img/T_iter_500.png'
set title 'it=500'
splot 'files/T_200.txt' with pm3d

set output 'img/T_iter_1000.png'
set title 'it=1000'
splot 'files/T_1000.txt' with pm3d

set output 'img/T_iter_2000.png'
set title 'it=2000'
splot 'files/T_2000.txt' with pm3d



set output 'img/NT_iter_100.png'
set title 'it=100'
splot 'files/NT_100.txt' with pm3d

set output 'img/NT_iter_200.png'
set title 'it=200'
splot 'files/NT_200.txt' with pm3d

set output 'img/NT_iter_500.png'
set title 'it=500'
splot 'files/NT_500.txt' with pm3d

set output 'img/NT_iter_1000.png'
set title 'it=1000'
splot 'files/NT_1000.txt' with pm3d

set output 'img/NT_iter_2000.png'
set title 'it=2000'
splot 'files/NT_2000.txt' with pm3d