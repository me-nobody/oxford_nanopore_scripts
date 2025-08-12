#!/usr/bin/env nextflow

process ALIGN { 
    beforeScript '''\
        module purge; module load bluebear
        module load bear-apps/2023a/live
        module load dorado/0.9.0-foss-2023a-CUDA-12.1.1
    '''.stripIndent()  

    
    input: 
    path reference_genome
    each path(basecalled_bam_file)
    path dir_path


    output: 
    path 'chunk_*' 

    script: 
    """
    dorado aligner  ${reference_genome} ${basecalled_bam_file} --output-dir ${dir_path}
    """
} 



workflow { 
         def dir_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data'
         reference = channel.fromPath("/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/human_pangenome_reference_v1.1.fasta",checkIfExists:true)
         bam_files = [file('/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output/EBT_basecalled_100625.bam'),
                      file('/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output/positive_control_basecalled_100625.bam') ]
         ALIGN(reference,bam_files,dir_path)
}