import pysam
import os 
import time

import logging
logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(
        filename="bam_analysis.log",
        encoding="utf-8",
        filemode="a",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
    )

dir_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/output_data'
bam_file = 'positive_control_T2T_sorted_120625.bam'
index_file = 'positive_control_T2T_sorted_120625.bam.bai'


# dir_path = "/rds/projects/b/broderra-mrc-alt-telomere/Anu/karlseder_data"
# bam_file = 'SRR28825768_sorted.bam'
# index_file = 'SRR28825768_sorted.bai'

# dir_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/'
# bam_file = 'SRR28825768_sorted.bam'
# index_file = 'SRR28825768_sorted.bai'

file_handler = logging.FileHandler("bam_analysis_"+bam_file[:-4]+".log", mode="a", encoding="utf-8")
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
logger.addHandler(file_handler)

if os.path.exists(dir_path):
    align_file = os.path.join(dir_path,bam_file)
    align_index_file = os.path.join(dir_path,index_file)

contig_list = ['chr'+str(i) for i in range(1,23)]
contig_list = contig_list + ['chrX','chrY','chrM'] # chrM is ignored for calculation
logger.info(f'{len(contig_list)} chromosomes in reference')
# file to read from
samfile = pysam.AlignmentFile(align_file,'rb',index_filename = align_index_file)
# header of samfile
# header = pysam.AlignmentHeader.copy(samfile)
# file to write to
telofile = pysam.AlignmentFile(os.path.join(dir_path,bam_file[:-4]+"_telomere.bam"), "wb", template=samfile)
logger.info(f'number of reads {samfile.get_index_statistics()}')

for contig in contig_list:    # iterate over chromosomes
    logger.info(f'                                                        ') # line gap
    if contig != 'chrM':
        upstream = downstream = 100000
        contig_end = samfile.get_reference_length(contig)
        kb100_downstream = contig_end-downstream
        logger.info(f'number of upstream reads {samfile.count(contig,0,100000)} for {contig}')
        logger.info(f'number of downstream reads {samfile.count(contig,kb100_downstream,contig_end)} for {contig}')
        upstream_reads = samfile.fetch(contig,0,upstream)
        ucount = 0
        for read in upstream_reads:       
            read_length = read.query_length # some read lenghts are returning zero
            if read_length >= 1000:
                logger.info(f'upstream read {read.query_name} length {read_length}')
                telofile.write(read)
                # print(read.to_string())
                # time.sleep(1)
                ucount += 1
        logger.info(f'{ucount} upstream reads  for {contig}')        

        downstream_reads = samfile.fetch(contig,kb100_downstream,contig_end)
        dcount = 0
        for read in downstream_reads:       
            read_length = read.query_length # some read lenghts are returning zero
            if read_length >= 1000:
                # logger.info(f'downstream read {read.query_name} length {read_length}')
                telofile.write(read)
                # print(read.to_string())
                # time.sleep(1)
                dcount +=1
        logger.info(f'{dcount} downstream reads  for {contig}')    
    elif contig == 'chrM':
        contig_end =  samfile.get_reference_length(contig)
        # logger.info(f'number of reads in mitochondria {samfile.count(contig)}')   
    else:
        # logger.info("read error")    
        print('pass')
    print(ucount)    
    print(dcount)    
# for contig in contig_list:
#     reads = samfile.fetch(contig,0,100000)
#     for read in reads:       
#         logger.info(f'chromosome {contig}')
#         logger.info(f'read {read.query_name} length {read.query_length}')
# close the file
samfile.close()
telofile.close()