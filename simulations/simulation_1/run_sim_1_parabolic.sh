#!/bin/bash

testing=$HOME/nonlinear-strpde
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
script=sim_1_parabolic.cpp

if [ -f script ]; then rm script; fi
g++ -o script $script -I$cpp -I$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s

for i in $(seq 0 2); do
    for j in $(seq 0 29); do
        ./script $i $j
    done
done
