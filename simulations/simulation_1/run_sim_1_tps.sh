#!/bin/bash

for i in $(seq 0 2); do
    for j in $(seq 0 29); do
    Rscript sim_1_tps.R $i $j
    done
done  
