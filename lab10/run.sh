#!/bin/bash

gcc verlet.c -lm &
wait
./a.out &
wait
python3 plots.py &
python3 maps.py &