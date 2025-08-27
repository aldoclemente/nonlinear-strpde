#!/bin/bash

#PBS -S /bin/bash
#PBS -l select=1:ncpus=2:mpiprocs=1,walltime=06:00:00
#PBS -j oe

APPTAINER=/opt/mox/apptainer/bin/apptainer
IMG=$HOME/fdapde-docker_latest.sif
simpath=$HOME/jmva/simulations/simulation_1/
cd $simpath

$APPTAINER exec $IMG ./script_1_nonlin $sigma $sim
