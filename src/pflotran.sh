#!/bin/bash -l
nreaz=$1
ncore=$2
pflotran_exe=$3

name="1dthermal"
dir="pflotran/"

cd $dir

srun -n $ncore $pflotran_exe -pflotranin $name".in" -stochastic -num_realizations $nreaz -num_groups $nreaz -screen_output off ;


wait

cd ../
