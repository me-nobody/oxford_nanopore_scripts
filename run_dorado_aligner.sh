#!/bin/bash -l

############# SLURM SETTINGS ##############
#SBATCH --account=broderra-mrc-alt-telomere   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=18      # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node
#SBATCH --job-name=drd_algn     # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --time=6:00:00           # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
############# LOADING MODULES #############
module purge
module load bluebear
module load bear-apps/2023a/live
module load --ignore_cache dorado/0.9.0-foss-2023a-CUDA-12.1.1
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
echo $(lspci | grep -i nvidia)
echo $(nvidia-smi -L)
echo "cuda visible devices $CUDA_VISIBLE_DEVICES"
echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/karlseder_data"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models" 
reference_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes"
echo $(ls -lh $dir_path/)
# dorado aligner $reference_path/T2Tconsortium_v2.0.fasta $dir_path/EBT_U2OS_15112023/output/EBT_bonito_basecalled_no_aln.bam >  $dir_path/EBT_U2OS_15112023/output/EBT_bonito_basecalled_aligned.bam
       
dorado aligner $reference_path/T2Tconsortium_v2.0.fasta $dir_path/SRR28825768.fastq >  $dir_path/SRR28825768_aligned.bam


# dorado aligner $reference_path/T2Tconsortium_v2.0.fasta $dir_path/Nanopore_271123/output/positive_control_bonito_basecalled_no_aln.bam > $dir_path/Nanopore_271123/output/positive_control_bonito_basecalled_aligned.bam


############# END #############
echo "this means  dorado aligner ran"
