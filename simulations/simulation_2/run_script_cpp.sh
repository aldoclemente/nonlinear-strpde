#!/bin/bash

testing=$HOME/nonlinear-strpde
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core

if [ -f script ]; then rm script; fi
g++ -o script $1 -I$cpp -I$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s
./script
