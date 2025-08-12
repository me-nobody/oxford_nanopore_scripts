#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=dnascent_CPU    # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt     # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt      # error file name will contain job name + job ID
#SBATCH --ntasks=18                # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --qos=bbdefault
#SBATCH --time=6:00:00             # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem-per-cpu=4G           # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1                  # 2 nodes, with each node having a GPU
#SBATCH --cpus-per-task=2          # BlueBEAR has 4 GPU in GPU node with 72 CPU per GPU

############# LOADING MODULES (optional) #############
module purge
module load bluebear
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############


echo $PWD
bam_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output_data"
container_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/containers" 
reference_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes"
index_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/dnascent"

echo $(wc -l $reference_path/T2Tconsortium_v2.0.fasta)

apptainer run  $container_path/DNAscent.sif detect \
              -b $bam_path/positive_control_sup_T2T_aligned.bam  \
              -r "$reference_path/T2Tconsortium_v2.0.fasta" \
              -i "$index_path/posctl_index.dnascent" \
              -o "$index_path/posctl_271123_CPU.bam" \
              -t 18 -q 20 -l 1000 --GPU 0
    
############# END #################################
echo "this means  DNAscent ran"

