#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks=18            # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1     # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --job-name=bnbc6.3   # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --time=1:0:00           # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############
module purge
module load bluebear
module load bear-apps/2023a
module load Bonito/0.8.1-foss-2023a-CUDA-12.1.1
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
# export CUDA_VISIBLE_DEVICES=0
echo $(lspci | grep -i nvidia)
echo $(nvidia-smi -L)
echo "cuda visible devices $CUDA_VISIBLE_DEVICES"
echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/bonito_models" 
reference_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes"
echo $(ls -lh $dir_path/EBT_U2OS_15112023/output_pod5s/)
# bonito basecaller $model_path/teloseq_bonito.0.6.2_R941C.fine_tuned $dir_path/EBT_U2OS_15112023/output_pod5s/ \
#        --reference $reference_path/T2Tconsortium_v2.0.fasta \
#        --max-reads 1000 > $dir_path/EBT_U2OS_15112023/output/EBT__bonito_basecalled_1000reads.bam      

bonito basecaller $model_path/teloseq_bonito.0.6.2_R941C.fine_tuned $dir_path/Nanopore_271123/output_pod5s/ > $dir_path/Nanopore_271123/output/positive_control_bonito_basecalled_no_aln.bam      

############# END #############
echo "this means  bonito basecaller_6.3 ran"
