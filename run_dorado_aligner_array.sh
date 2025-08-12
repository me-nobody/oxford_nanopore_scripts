#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --array=1-2
#SBATCH --qos=bbshort
#SBATCH --ntasks=8                # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --job-name=dorado_aligner   # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --time=0:9:00           # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############
module purge
module load bear-apps/2023a/live
module load --ignore_cache dorado/0.9.0-foss-2023a-CUDA-12.1.1
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #####################
echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models" 
echo $(ls -lh $dir_path/EBT_U2OS_15112023/fast5/)

file_list=("$dir_path/EBT_U2OS_15112023/output_pod5s","$dir_path/Nanopore_271123/output_pod5s")

INPUT_FILENAME=${file_list[${SLURM_ARRAY_TASK_ID}]}  # Look-up using array index

echo "I am array index ${SLURM_ARRAY_TASK_ID} and am processing file: $file_list[${SLURM_ARRAY_TASK_ID}]"

# dorado aligner hac $dir_path/EBT_U2OS_15112023/output_pod5s/ \
#        --device cuda:0 \
#        --models-directory $model_path \
#        --output-dir $dir_path/EBT_U2OS_15112023/output

############# END #############
echo "this means  dorado basecaller ran"
