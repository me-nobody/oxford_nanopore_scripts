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

dir_path = '/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/EBT_U2OS_15112023/output_data'
bam_file = 'EBT_telomere.bam'
index_file = 'EBT_telomere.bam.bai'

file_handler = logging.FileHandler("read_sequence_"+bam_file[:-4]+".log", mode="a", encoding="utf-8")
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler) # report to screen
logger.addHandler(file_handler)    # report to file

if os.path.exists(dir_path): # check if file path is true
    align_file = os.path.join(dir_path,bam_file)
    align_index_file = os.path.join(dir_path,index_file)

samfile = pysam.AlignmentFile(align_file,'rb',index_filename = align_index_file)

with open(os.path.join(dir_path,"seq_file.txt"),"a") as fh: # create output file path
    for read in samfile:
        name = ">"+read.query_name
        logger.info(f'name of the read {name}')
        seq = read.get_forward_sequence()
        seq = seq[1:200] # we just need the terminal end sequence
        logger.info(f'length of read {len(seq)}')
        fh.write(name+'\n')
        fh.write(seq)
        fh.write('\n')

# fh.seek(0) # to return the seek back to start        