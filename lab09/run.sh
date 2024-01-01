#!/bin/bash

gcc crank-nicolson.c -lm -lgslcblas -lgsl &
wait
./a.out &
wait
echo 'rysowanie map...'
gnuplot plots.plot &
gnuplot gifs.plot &
wait
echo 'wygenerowane obrazy znajduja sie w /img'