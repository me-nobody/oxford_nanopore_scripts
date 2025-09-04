#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=forksense_GPU    # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt     # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt      # error file name will contain job name + job ID
#SBATCH --ntasks=18               # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:1
#SBARCH --mem=200G
#SBATCH --time=1:00:00             # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm



############# LOADING MODULES (optional) #############
module purge
module load bluebear
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############

echo $PWD

# RUN THIS SCRIPT WITHIN THE DNASCENT FOLDER

echo $(wc -l cam_ont.detect)

apptainer run  --nv ../../containers/dnascent_3.1.2.sif forkSense \
              --detect cam.ont.complete.detect  \
              --output cam.ont.complete.bed \
              --order "BrdU,EdU" \
              --markAnalogues \
              --markTerminations \
              --markForks \
              --makeSignatures \
              --threads 18
              
                  
############# END #################################
echo "this means  DNAscent Forksense ran"

