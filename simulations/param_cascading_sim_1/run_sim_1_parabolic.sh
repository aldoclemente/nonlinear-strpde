#!/bin/bash

lbfgspp=/usr/include/LBFGSpp/include
eigen=/usr/include/eigen3

testing=$HOME/nonlinear-strpde
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
script=sim_1_parabolic.cpp
exe=$(basename $script .cpp)

if [ -f $exe ]; then rm script; fi
g++ -o $exe $script -I$cpp -I$core -I$eigen -I$lbfgspp -O2 -std=c++20 -march=native -s

for i in $(seq 0 2); do
    for j in $(seq 0 29); do
        ./$exe $i $j
    done
done
