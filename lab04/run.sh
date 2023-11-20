#!/bin/bash

gcc poisson_global.c -lm &
wait
./a.out &
wait
python3 plot_glob.py &
gcc poisson_local.c -lm &
wait
./a.out &
wait
python3 plot_local.py &
wait
echo 'END'