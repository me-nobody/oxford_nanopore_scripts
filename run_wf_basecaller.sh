#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=broderra-mrc-alt-telomere
#SBATCH --job-name=wf_bc_GPU       # some descriptive job name of your choice
#SBATCH --output=%x-%j_out.txt     # output file name will contain job name + job ID
#SBATCH --error=%x-%j_err.txt      # error file name will contain job name + job ID
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=18                # each task is a thread and make changes to script too
#SBARCH --mem=200G
#SBATCH --time=6:00:00             # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
         
        

############# LOADING MODULES (optional) #############
module purge
module load bluebear
module load bear-apps/2022b/live
module load Nextflow/24.04.2
echo "Hello from $SLURM_JOB_NODELIST"
############# MY CODE #############
input_files="/rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/2018_08_16_CAM_ONT_2085_1cycle_B_ligation"
model_path="/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models"
reference="/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta"

# nextflow run epi2me-labs/wf-basecalling \
#     --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0' \
#     --dorado_ext 'fast5' \
#     --input $input_files \
#     --ref $reference \
#     --output_fmt 'bam' \
#     --basecaller_model_path "${model_path}/dna_r9.4.1_e8_sup@v3.6" \
#     --out_dir "/rds/projects/b/broderra-mrc-alt-telomere/Anu/wf-basecalling/data" \
#     -profile standard

# apptainer updated in the base.config
# script run within wf-basecalling as it has the base.config with apptainer updated
# quick multi-liner to convert FAST5 to POD5

nextflow run epi2me-labs/wf-basecalling \
    --dorado_ext 'fast5' \
    --input "../boemo_lab_data/cam_ont_multiread/" \
    --ref "../reference_genomes/T2Tconsortium_v2.0.fasta" --output_fmt 'bam' \
    --basecaller_model_path "../dorado_models/dna_r9.4.1_e8_sup@v3.3/" \
    --out_dir "./data" -profile apptainer -resume --cuda_device 'cuda:all'
