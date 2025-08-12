#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=forksense_CPU    # some descriptive job name of your choice
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
detect_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/dnascent"
container_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/containers" 



echo $(wc -l $detect_path/ebt_u2os_15112023_CPU.detect)

apptainer run  $container_path/DNAscent.sif forkSense \
              --detect $detect_path/ebt_u2os_15112023_CPU.detect  \
              --output "$detect_path/ebt_u2os_CPU_forksense.bed" \
              --order "EdU,BrdU" \
              --markAnalogues \
              --markTerminations \
              --markForks \
              --makeSignatures \
              --threads 18
              
                  
############# END #################################
echo "this means  DNAscent Forksense ran"

