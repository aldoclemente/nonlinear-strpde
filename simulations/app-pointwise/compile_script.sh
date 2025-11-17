#!/bin/bash

testing=$HOME/nonlinear-strpde
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
exe=$(basename -- $1 .cpp)

if [ -f $exe ]; then rm $exe; fi
g++ -o $exe $1 -I$cpp -I$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s
