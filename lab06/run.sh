#!/bin/bash

gcc rho_0.c -lm &
wait
./a.out &
wait
gcc rho_n_0.c -lm &
wait
./a.out &
wait
python3 plots.py &
wait
echo ' '