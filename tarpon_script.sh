#!/usr/bin/bash

dir_path='/rds/projects/b/broderra-mrc-alt-telomere/Anu'
# nextflow run $dir_path/TARPON/main.nf --input $dir_path/pilot_data/EBT_U2OS_15112023/output_data/EBT_sup_basecalled.bam \
#     --capture_probe_sequence TTTTCCTGTACTTCGTTCAGTTACGTATTGCT \
#     --sample_name EBT_U2OS_15112023 \
#     --outdir $dir_path/dump \
#     -profile singularity

nextflow run $dir_path/TARPON/main.nf --input $dir_path/pilot_data/EBT_U2OS_15112023/output_data/EBT_sup_basecalled.bam \
    --capture_probe_sequence TTTTCCTGTACTTCGTTCAGTTACGTATTGCT \
    --sample_name EBT_U2OS_15112023 \
    --outdir $dir_path/dump \
    -profile singularity
