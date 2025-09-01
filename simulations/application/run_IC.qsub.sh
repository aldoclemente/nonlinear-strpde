#!/bin/bash

#PBS -S /bin/bash
#PBS -l select=1:ncpus=96:mpiprocs=48,walltime=48:00:00
#PBS -j oe

APPTAINER=/opt/mox/apptainer/bin/apptainer
IMG=$HOME/fdapde-docker_latest.sif
simpath=$HOME/jmva/simulations/application/
cd $simpath

date
$APPTAINER exec $IMG ./script_IC
date