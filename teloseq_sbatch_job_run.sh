#!/bin/bash
#SBATCH --job-name=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4  # match this to the number of cores you want to run the job with
#SBATCH --time=04:00:00  # this is arbitrarily set here to 4 hours
#SBATCH --mem=16G  # for each ntask, 4GB RAM is assigned. But I <think> snakemake doesn't make good use of ntasks so we'll specify the RAM. Think in blocks of 4G
#SBATCH --qos=bbdefault
#SBATCH --account=your_account_here
#SBATCH --mail-type=FAIL

set -e

# Load environment modules
module purge; module load bluebear
module load bear-apps/2022b
module load Miniforge3/24.1.2-0
module load bear-apps/2023a/live
module load dorado/0.9.0-foss-2023a-CUDA-12.1.1

# start interactive session
module load slurm-interactive
fisbatch_screen --nodes=1-1 --time=3:0:0 --mem=36G --ntasks=16

# Enable mamba/conda
eval "$(${EBROOTMINIFORGE3}/bin/conda shell.bash hook)"
source "${EBROOTMINIFORGE3}/etc/profile.d/mamba.sh"

# Activate your environment
CONDA_ENV_PATH="/rds/projects/b/broderra-mrc-alt-telomere/${USER}_conda_env"
mamba activate "${CONDA_ENV_PATH}"

# Navigate to the pipeline directory
cd /rds/projects/b/broderra-mrc-alt-telomere/Anu/teloseq

# Run Snakemake with the number of cores assigned by Slurm
snakemake --use-conda --configfile config.yml -pr --cores "$SLURM_CPUS_PER_TASK"

snakemake --use-conda --configfile config.yml -pr --cores 16  --ri  --detailed-summary
snakemake --use-conda --configfile config.yml -pr --cores 16  --ri >output.txt 2>&1

# to only get the debug report
snakemake --cores 16 --debug-dag >dag.report.txt 2>&1
# to create html report
snakemake --cores 16 --dryrun --report report.html

