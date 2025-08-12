#!/bin/bash
# create a interactive node
module load slurm-interactive
fisbatch_screen --nodes=1-1 --time=3:0:0 --mem=36G --ntasks=16

# Activate your environment
CONDA_ENV_PATH="/rds/projects/b/broderra-mrc-alt-telomere/${USER}_conda_env"
mamba activate "${CONDA_ENV_PATH}"

# change directory
cd /rds/projects/b/broderra-mrc-alt-telomere/Anu/teloseq

# Run Snakemake with the number of cores assigned by Slurm
snakemake --use-conda --configfile config.yml -pr --cores 16  --ri  --detailed-summary --report

snakemake --use-conda --configfile config.yml --cores 16 --rerun-incomplete >output.txt 2>&1

snakemake --use-conda --configfile config.yml --cores 16  --ri  --detailed-summary --report

