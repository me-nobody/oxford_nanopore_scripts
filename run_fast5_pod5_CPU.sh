#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --qos=bbdefault
#SBATCH --nodes=1                    # number of nodes to allocate, default is 1
#SBATCH --ntasks=24                 # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1           # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --mem=200G                   # memory required per node, in the form of [num][M|G|T]
#SBATCH --job-name=pod5-file-format  # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt       # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt        # error file name will contain job name + job ID
#SBATCH --time=1-10:00:00            # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############

set -e

module purge
module load bear-apps/2023a
module load pod5-file-format/0.3.23-foss-2023a
echo "${SLURM_JOB_ID}: Job ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX} in the array"
############# MY CODE #############

echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/2018_08_16_CAM_ONT_2085_1cycle_B_ligation/"
for file in $(ls /rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/2018_08_16_CAM_ONT_2085_1cycle_B_ligation/);
do 
    pod5 convert fast5 /rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/2018_08_16_CAM_ONT_2085_1cycle_B_ligation/$file  --output ./cam_ont_output_pod5s/${file:1:-6}.pod5 
done 



############# END #############
echo "this means  pod5 file format ran"
