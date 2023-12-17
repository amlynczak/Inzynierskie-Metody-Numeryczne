#!/bin/bash

gcc crank-nicolson.c -lm &
wait
./a.out &
wait
echo 'rysowanie wykresow, map oraz gifow...' 
python3 maps.py &
python3 plot_c.py &
python3 plot_xsr.py &
gnuplot plot5.plt &
wait
echo 'wygenerowane gify: D0.gif oraz D01,gif'