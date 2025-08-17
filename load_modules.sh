# project_path
# /rds/projects/b/broderra-mrc-alt-telomere/Anu



# pull the nextflow wf-teloseq workflow from Oxford Nanopore Docker hub
apptainer pull docker://ontresearch/wf-teloseq
# create a interactive session
module load slurm-interactive
fisbatch_screen --nodes=1-1 --ntasks=16 --time=3:0:0 --mem=36G

# ont-fast5-api
module purge;module load bluebear
module load bear-apps/2022a
module load ont-fast5-api/4.1.1-foss-2022a

#PYSAM
module load bear-apps/2024a
module load SAMtools/1.21-GCC-13.3.0
module load Pysam/0.22.1-GCC-13.3.0

# IGV
module load bear-apps/2024a
module load Mesa/24.1.3-GCCcore-13.3.0
module load IGV/2.19.5-Java-21


# Accessing pod5-file-format 0.3.23-foss-2023amosule
module load bear-apps/2023a
module load pod5-file-format/0.3.23-foss-2023a

# subsample a bamfile for development and testing
module load bear-apps/2022a/live
module load SAMtools/0.1.20-GCC-11.3.0
#
module load bear-apps/2023a/live
module load SAMtools/1.21-GCC-12.3.0

# find partially aligned reads in bam file
samtools view yourfile.bam | awk '$6 ~ /[SH]/ {print $1, $6}'


# samtools sort bam file
samtools sort -@16 file > sorted_file.bam

# samtools index bam file
samtools index -@16 -b -o sorted_file.bam.bai sorted_file.bam

# I have chosen EBT_U2OS_15112023 for development and testing
# select a small chromsome( chr 22 ) for testing the pipeline
samtools view -b -o EBT_HP_chr22.bam EBT_HP_aligned_130625.bam chr22

samtools stats <file.bam> |grep ^SN|cut -f 2-

# converting the bam file to FastQ as required by wf-teloseq. It actually requires only the basecalled unaligned data
samtools fastq --threads=8 -o EBT_HP_chr22.fq EBT_HP_chr22.bam

# extract reads with longer length
samtools view -h aligned.bam | awk 'length($10) > 30 || $1 ~ /^@/' | samtools view -bS - > filtered.bam

# samtools convert unaligned FASTQ to bam file
samtools import <input.fastq> -o <output.bam>

awk 'match($0, /L[[:alpha:]]T/) {
print RSTART, substr($0, RSTART, RLENGTH)}' file

samtools view <infile.bam> |awk 'match($10,/pattern/){print RSTART,substr($10,RSTART,RLENGTH)}'
# improved version
awk '{
   n = 0
   while (match($0, /L[[:alpha:]]T/)) {
      n += RSTART
      print n, substr($0, RSTART, RLENGTH)
      $0 = substr($0, RSTART + 1)
   }
}' file

# import nextflow module
module load bear-apps/2022b/live
module load Nextflow/24.04.2

# loading and running epi2me basecalling and aligning workflow
nextflow run epi2me-labs/wf-basecalling --help


# load DORADO
module load bear-apps/2023a/live
module load --ignore_cache dorado/0.9.0-foss-2023a-CUDA-12.1.1
# to run nextflow workflow we need a 1) script 2) config file
# nextflow run hello.nf -c custom.config
nextflow run main.nf -c nextflow_teloseq.config \
 --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output/EBT_HP_chr22.fq \
 --skipmapping -profile apptainer
# config file has skipmapping set to false and also path to reference file
nextflow run main.nf -c nextflow.config \
 --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output/EBT_sup_T2T_basecalled.fastq \
 --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta   -profile apptainer

nextflow run main.nf -c nextflow.config \
 --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output_data/EBT_sup_T2T_basecalled.fastq \
 --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta  --skipmapping -profile apptainer

nextflow run main.nf -c nextflow.config \
 --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output_data/positive_control_sup_T2T_basecalled.fastq \
 --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta  --skipmapping -profile apptainer


nextflow run main.nf -c nextflow.config --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output_data/EBT_bonito_basecalled.fastq \
  --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta  -profile apptainer

# epi2me workflow
nextflow run main.nf -c nextflow_epi2me.config --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/karlseder_data/SRR28825768.fastq \
  --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta  -profile apptainer -resume

# nanopore workflow
nextflow run main.nf -c nextflow.config --fastq /rds/projects/b/broderra-mrc-alt-telomere/Anu/karlseder_data/SRR28825768.fastq \
  --reference /rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/T2Tconsortium_v2.0.fasta  \
  --skipmapping -profile apptainer

# run bonito view 
module load bear-apps/2023a
module load Bonito/0.8.1-foss-2023a-CUDA-12.1.1
bonito view /rds/projects/b/broderra-mrc-alt-telomere/Anu/bonito_models/teloseq_bonito.0.6.2_R941C.fine_tuned

# Slurm script to view batch jobs
sacct --user=dasaz --format=jobname,avecpu,cputime,ncpus,ntasks,reqcpu,reqmem

# TARPON
--outdir_positive "/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/tarpon_outdir"
--indir_positive "/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output"
--indir_positive_f1 = "positive_control_sup_T2T_basecalled.fastq"
--indir_positive_f2 = "positive_control_hac_T2T_basecalled.fastq"

--outdir_EBT "/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/tarpon_outdir"
--indir_EBT "/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output"
--indir_EBT_f1 = "EBT_hac_T2T_basecalled.fastq"
--indir_EBT_f2 = "EBT_sup_T2T_basecalled.fastq"

# time for jobs to wait in a queue excluding interactive shell
sacct -X --starttime YYYY-MM-DD -o JobID,JobNAme,Planned |sort -u -k3 -n|grep -v FISBATCH

# arrange according to file type
ls -alXh

# wf-teloseq did not run properly with Bonito custom model basecalled data, even with lowered parameters. I could see alignments in IGV but these were
# mostly in locations other than Dorado basecalled data.

# pycoQC
module load pycoQC/2.5.0.21-foss-2019b-Python-3.7.4
# mutliqc
moodule load bear-apps/2019b/live
module load MultiQC/1.9-foss-2019b-Python-3.7.4
module load FastQC/0.11.9-Java-11
fastqc -t16 --noextract --dir /scratch/dasaz -o ../output_summary/EBT_sup_T2T_basecalled.fastq

module load bear-apps/2023a
module load NanoPlot/1.43.0-foss-2023a

# SRA toolkit
module load bear-apps/2023a
module load SRA-Toolkit/3.0.10-gompi-2023a


# remove symbolic links
find -L . -name . -o -type d -prune -o -type l -exec rm {} +

# search for a string in all files in a folder
grep -rni "string" *
