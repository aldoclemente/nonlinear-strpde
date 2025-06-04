#!/bin/bash

for i in $(seq 0 20); do 
    for j in $(seq 0 29); do
        ./script 1 $j $i
    done
done

