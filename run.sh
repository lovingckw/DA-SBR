#!/bin/bash -l
#SBATCH -A m1800
#SBATCH -q regular
#SBATCH -N 3
#SBATCH -t 00:30:00
#SBATCH -L SCRATCH
#SBATCH -J esmda-flux
#SBATCH -C haswell

module load python/3.6-anaconda-4.4
python src/esmda.py;
wait
