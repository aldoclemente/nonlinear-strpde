#!/bin/bash

core=$HOME/Desktop/fdaPDE-core
cpp=$HOME/Desktop/fdaPDE-cpp
nonlinear=$HOME/Desktop/nonlinear-strpde/simulations
docker run --rm -v $cpp:/root/fdaPDE-cpp -v $nonlinear:/root/nonlinear-strpde -ti aldoclemente/fdapde-docker /bin/bash
