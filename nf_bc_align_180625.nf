#!/usr/bin/env nextflow


process BASECALL { 
    beforeScript '''\
        module purge; module load bluebear
        module load bear-apps/2023a/live
        module load dorado/0.9.0-foss-2023a-CUDA-12.1.1
    '''.stripIndent()  

    
    input: 
    path model_path
    each path(pod5_file)
    path dir_path

    output: 
    path 'chunk_*' 

    script: 
    """
    dir_path=''
    dorado basecaller sup ${pod5_file} --models-directory ${model_path} --output-dir ${dir_path}
    """
} 


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
    dir_path=''
    dorado aligner  ${reference_genome} ${basecalled_bam_file} --output-dir ${dir_path}
    """
} 


workflow { 
         def dir_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data'
         def model_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/dorado_models'
         def reference = channel.fromPath("/rds/projects/b/broderra-mrc-alt-telomere/Anu/reference_genomes/human_pangenome_reference_v1.1.fasta")
         pod5_file = [file('/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output_pod5s'),
                      file('/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output_pod5s') ]
         BASECALL(model_path,pod5_file,dir_path)
}