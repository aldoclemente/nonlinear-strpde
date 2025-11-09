#!/bin/bash

locs=1000
date > log.txt
for i in $(seq 0 29); do
    echo $i 
    ./simulazione_ic_stimata $locs $i
done 
date >> log.txt
