#!/bin/bash
#SBATCH --account=broderra-mrc-alt-telomere   # account name 
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks=4             
#SBATCH --cpus-per-task=2      # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs

#SBATCH --job-name=nf_aligner     
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --time=1:30:00           # time limit for the whole run, in the form of d-hh:mm:ss

module purge; module load bluebear
module load bear-apps/2022b
module load Nextflow/24.04.2

nextflow run nf_aligner_130625.nf