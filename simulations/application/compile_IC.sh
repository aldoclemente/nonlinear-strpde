#!/bin/bash

testing=$HOME/nonlinear-strpde
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
script=IC.cpp

if [ -f script ]; then rm script; fi
g++ -o script_IC $script -I$cpp -I$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s
chown 1000:1000 script_IC 