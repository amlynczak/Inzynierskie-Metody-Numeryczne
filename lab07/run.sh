#!/bin/bash

gcc navier_stokes.c -lm &
wait
./a.out &
wait
echo 'rysowanie wykresow...'
gnuplot psi.plot &
gnuplot dzeta.plot &
python3 maps.py &
wait
echo 'wykresy konturowe sa w folderze ./img'