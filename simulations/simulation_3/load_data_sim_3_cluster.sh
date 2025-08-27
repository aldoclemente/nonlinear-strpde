#!/bin/bash

declare -a sigma=("0.00" "0.05" "0.10")

scp $(pwd)/../data/simulation_3/incidence_matrix.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_3/    

scp $(pwd)/script_3_nonlin clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_3/
scp $(pwd)/run_sim_3_nonlin.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_3/
scp $(pwd)/run_sim_3_nonlin.qsub.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_3/

for i in "${sigma[@]}"; do
    for j in $(seq 0 29); do
        scp -r $(pwd)/../data/simulation_3/sigma_$i/$j/obs.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_3/sigma_$i/$j/obs.mtx
    done
done

