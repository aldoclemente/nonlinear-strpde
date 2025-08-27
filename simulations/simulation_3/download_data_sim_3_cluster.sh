#!/bin/bash

declare -a sigma=("0.00" "0.05" "0.10")

for i in "${sigma[@]}"; do
    for j in $(seq 0 29); do
        scp clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_3/sigma_$i/$j/estimate_nonlinear.mtx $(pwd)/../data/simulation_3/sigma_$i/$j/estimate_nonlinear.mtx
        scp clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_3/sigma_$i/$j/rmse_nonlinear.mtx $(pwd)/../data/simulation_3/sigma_$i/$j/rmse_nonlinear.mtx
    done
done

