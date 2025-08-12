#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=18      # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node
#SBATCH --job-name=dorado_basecaller   # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --time=2:00:00           # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############
module purge
module load bear-apps/2023a/live
module load dorado/0.9.0-foss-2023a-CUDA-12.1.1
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
# export CUDA_VISIBLE_DEVICES=0
echo $(lspci | grep -i nvidia)
echo $(nvidia-smi -L)
echo "cuda visible devices $CUDA_VISIBLE_DEVICES"
echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models" 
bonito_model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/bonito_models/teloseq_bonito.0.6.2_R941C.fine_tuned" 
echo $(ls -lh $dir_path/EBT_U2OS_15112023/fast5/)
dorado basecaller $bonito_model_path $dir_path/EBT_U2OS_15112023/output_pod5s/ \
       --device cuda:0 \
       --models-directory $bonito_model_path > $dir_path/EBT_U2OS_15112023/output/EBT_bonito_model_dorado_basecalled.bam
       

# dorado basecaller sup $dir_path/Nanopore_271123/output_pod5s/ \
#        --device cuda:0 \
#        --models-directory $model_path > $dir_path/Nanopore_271123/output/positive_sup_basecalled.bam
        

############# END #############
echo "this means  dorado basecaller ran"
