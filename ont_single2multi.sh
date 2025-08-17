#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --qos=bbdefault
#SBATCH --nodes=1                    # number of nodes to allocate, default is 1
#SBATCH --ntasks=24                  # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1            # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --mem=200G                   # memory required per node, in the form of [num][M|G|T]
#SBATCH --job-name=ont-fast5-api  # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt       # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt        # error file name will contain job name + job ID
#SBATCH --time=1-10:00:00            # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############

set -e

module purge;module load bluebear
module load bear-apps/2022a
module load ont-fast5-api/4.1.1-foss-2022a
echo "${SLURM_JOB_ID}: Job ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX} in the array"
############# MY CODE #############

echo $PWD
single_to_multi_fast5 -i  /rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/2018_08_16_CAM_ONT_2085_1cycle_B_ligation \
                      -s   /rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/cam_ont_multiread \
                      -t 24 -n 100 


# option at commandline worked
# single_to_multi_fast5 -i ./2018_08_16_CAM_ONT_2085_1cycle_B_ligation/ -s ./cam_ont_multiread/ -n100 -t16                      


