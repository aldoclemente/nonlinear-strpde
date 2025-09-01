#!/bin/bash

scp $(pwd)/../data/application//mesh/incidence_matrix.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/mesh/
scp $(pwd)/../data/application/mesh/nodes.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/mesh/
scp $(pwd)/../data/application/mesh/cells.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/mesh/
scp $(pwd)/../data/application/mesh/boundary.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/mesh/

scp $(pwd)/../data/application/obs0.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/
scp $(pwd)/../data/application/obs.csv clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/data/application/

scp $(pwd)/script_IC clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/
scp $(pwd)/script_amy-negative clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/

scp $(pwd)/run_app-amy-negative.qsub.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/
scp $(pwd)/run_app-amy-negative.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/
scp $(pwd)/run_IC.qsub.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/
scp $(pwd)/run_IC.sh clemente@kami.inside.mate.polimi.it:/u/clemente/jmva/simulations/application/