#!/bin/bash

cpp=$HOME/fdaPDE-cpp
dir=$(dirname $1)
script=$(basename $1)

cd $dir 
g++ -o script $script -I$cpp -I$cpp/fdaPDE/core/ -I/usr/include/eigen3 -O2 -std=c++20 -g -march=native -DFDAPDE_NO_DEBUG

