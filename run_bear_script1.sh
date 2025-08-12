#!/bin/bash

#SBATCH --nodes 4
#SBATCH --ntasks-per-node 8
#SBATCH --time 0-00:05:00 # days-hours:minutes:seconds
#SBATCH --qos bbshort  # bbshort allows up to 10 minutes for jobs

set -e

module purge 
module load bluebear
module load bear-apps/2022b
module load Python/3.10.8-GCCcore-12.2.0
module load mpi4py/3.1.4-gompi-2022b

mpiexec -n ${SLURM_NTASKS} python hello.py