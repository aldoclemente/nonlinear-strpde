#!/bin/bash

#PBS -S /bin/bash
#PBS -l select=1:ncpus=2:mpiprocs=1,walltime=48:00:00
#PBS -l mem=256gb
#PBS -j oe
#PBS -N app-init

#APPTAINER=/opt/mox/apptainer/bin/apptainer
#IMG=$HOME/fdapde-docker_latest.sif
#simpath=$HOME/jmva/simulations/application/
#cd $simpath

#date
#$APPTAINER exec $IMG ./script_IC
#date

cd "${PBS_O_WORKDIR}"

# clean log_ic
rm -f log_ic.txt compile_ic.txt

date >> log_ic.txt
echo "cleaning log_ic files..." >> log_ic.txt

# init spack
echo "loading spack" >> log_ic.txt
if [ ! -f "${HOME}/spack.sh" ]; then
    echo "spack configuration file not found" >> log_ic.txt
    exit 1
fi
source "${HOME}/spack.sh" >> log_ic.txt 2>&1

echo "loading fdapde_devel environment" >> log_ic.txt
spack env activate fdapde_devel >> log_ic.txt 2>&1

# check packages
EIGEN_LOC="$(spack location -i eigen)"
if [ -z "$EIGEN_LOC" ]; then
    echo "Eigen dependency not found" >> log_ic.txt
    exit 1
fi
echo "Eigen library found at ${EIGEN_LOC}" >> log_ic.txt

# compile
cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
script=IC.cpp
exe=$(basename -- "$script" .cpp)

echo "compiling..." >> log_ic.txt
g++ -o $exe $script \
    -I$cpp \
    -I$core \
    -I../include \
    -I"$EIGEN_LOC/include/eigen3" \
    -O2 -std=c++20 >> compile_ic.txt 2>&1

# run
echo "launching executable..." >> log_ic.txt
./$exe $N >> log_ic.txt 2>&1
date >> log_ic.txt
