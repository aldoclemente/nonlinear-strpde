#!/bin/bash

cpp=$HOME/fdaPDE-cpp
#dir=$(dirname $1)
#script=$(basename $1)
if [ -d $1 ]; then
    dir=$1
elif [ -f $1 ]; then
    dir=$(dirname $1)
fi

cd $dir

# $2 == 0, 1, 2, 3 (obs)
# $3 == 0, ..., 29 (sim)
# $4 == lambda ---> se non fornito lambda = 1
./script $2 $3 $4
