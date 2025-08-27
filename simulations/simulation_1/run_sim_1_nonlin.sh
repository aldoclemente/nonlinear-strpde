#!/bin/bash

mkdir -p qsuboutput

for i in $(seq 0 2); do
    for j in $(seq 0 29); do
        qsub -o qsuboutput/output-$i-$j.txt -e qsuboutput/error-$i-$j.txt -N sim_1-$i-$j -v sigma=$i,sim=$j run_sim_1_nonlin.qsub.sh
    done
done
