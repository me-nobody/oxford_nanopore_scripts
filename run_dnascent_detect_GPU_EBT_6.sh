#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=dnascent_GPU     # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt      # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt       # error file name will contain job name + job ID
#SBATCH --gres=gpu:a30:1            # GPU per node
#SBATCH --qos=bbgpu
#SBATCH --time=6:00:00              # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem-per-gpu=54G           # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1                   # 2 nodes, with each node having a GPU
#SBATCH --ntasks-per-gpu=1         # each GPU will handle 1 task. script needs multiparallel
#SBATCH --cpus-per-task=18          # BlueBEAR has 4 GPU in GPU node with 72 CPU per GPU

############# LOADING MODULES (optional) #############
module purge
module load bluebear
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############


echo $PWD
bam_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output_data"
container_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/containers" 
reference_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes"
index_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/dnascent"

echo $(wc -l $reference_path/T2Tconsortium_v2.0.fasta)

apptainer run --env CUDA_VISIBLE_DEVICES=all --nv $container_path/DNAscent.sif detect \
              -b $bam_path/EBT_sup_T2T_sorted.bam  \
              -r "$reference_path/T2Tconsortium_v2.0.fasta" \
              -i "$index_path/ebt_index.dnascent" \
              -o "$index_path/ebt_u2os_15112023.detect" \
              -t 18 --GPU 1
    
############# END #################################
echo "this means  DNAscent ran"

