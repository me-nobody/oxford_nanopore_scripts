#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=dnascent_GPU    # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt     # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt      # error file name will contain job name + job ID
#SBATCH --ntasks=20                # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=18                # each task is a thread and make changes to script too
#SBARCH --mem=200G
#SBATCH --time=6:00:00             # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm

############# LOADING MODULES (optional) #############
module purge
module load bluebear
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
# RUN THIS SCRIPT FROM PWD
echo $PWD

echo $(wc -l ../../reference_genomes/S288C_reference_sequence_R9-1-1_19990210.fsa)

# -q and -l is not needed as wf-basecaller has separated the high quality data

apptainer run --nv  ../../containers/dnascent_3.1.2.sif detect \
              -b ../cam_ont_aligned_bam/cam.ont.complete.bam  \
              -r ../../reference_genomes/S288C_reference_sequence_R9-1-1_19990210.fsa \
              -i cam_ont_index.dnascent \
              -o cam.ont.complete.detect -t 18  --GPU 0
    
############# END #################################
echo "this means  DNAscent detect ran"

