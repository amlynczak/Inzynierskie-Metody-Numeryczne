#!/bin/bash

gcc poisson_multigrid.c -lm &
wait
./a.out &
wait
python3 plots.py &
wait
echo ' ' 