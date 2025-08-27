#!/bin/bash

declare -a sigma=("0.00" "0.05" "0.10")

scp $(pwd)/../data/simulation_1/test_locs.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/
scp $(pwd)/../data/simulation_1/time_mesh.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/    

scp $(pwd)/script_1_nonlin clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_1/
scp $(pwd)/run_sim_1_nonlin.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_1/
scp $(pwd)/run_sim_1_nonlin.qsub.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/simulation_1/

for i in "${sigma[@]}"; do
    for j in $(seq 0 29); do
        scp -r $(pwd)/../data/simulation_1/sigma_$i/$j/obs.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/sigma_$i/$j/obs.mtx
        scp -r $(pwd)/../data/simulation_1/sigma_$i/$j/locs.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/sigma_$i/$j/locs.mtx
        scp -r $(pwd)/../data/simulation_1/sigma_$i/$j/exact.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/sigma_$i/$j/exact.mtx
        scp -r $(pwd)/../data/simulation_1/sigma_$i/$j/test_obs.mtx clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/simulation_1/sigma_$i/$j/test_obs.mtx
    done
done
