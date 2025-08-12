#!/bin/bash -l

############# SLURM SETTINGS #############
############# LOADING MODULES #############
module purge
module load bear-apps/2023a/live
module load dorado/0.9.0-foss-2023a-CUDA-12.1.1
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
# export CUDA_VISIBLE_DEVICES=0

echo $PWD
dir_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models" 
reference_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes"
echo $(ls -lh $dir_path/)
dorado aligner $reference_path/T2Tconsortium_v2.0.fasta $dir_path/EBT_U2OS_15112023/output_pod5s/ \     
       --threads=48 \  
       --output-dir $dir_path/EBT_U2OS_15112023/output

dorado aligner $reference_path/T2Tconsortium_v2.0.fasta $dir_path/Nanopore_271123/output_pod5s/ \     
       --threads=48 \  
       --output-dir $dir_path/Nanopore_271123/output


############# END #############
echo "this means  dorado basecaller ran"
