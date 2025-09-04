# create a interactive session
module load slurm-interactive
fisbatch_tmux --nodes=1-1 --ntasks=16 --time=3:0:0 --mem=36G

# load ont-fast5-api
module load bear-apps/2022a
module load ont-fast5-api/4.1.1-foss-2022a

# path=/rds/projects/b/broderra-mrc-alt-telomere/Anu/boemo_lab_data/cam_ont_multiread

# script
from ont_fast5_api.fast5_interface import *
# Path to your FAST5 file

fast5_path = "../boemo_lab_data/cam_ont_multiread/batch_421.fast5"

# f5_file = get_fast5_file(fast5_path,'r')
# read_ids = f5_file.get_read_ids()   
# reads = f5_file.get_reads()   

with get_fast5_file(fast5_path,'r'):
    f5_file = get_fast5_file(fast5_path,'r')
    read_ids = f5_file.get_read_ids() 
    reads = f5_file.get_reads()   # this is a generator
    for read in reads:
        context_tags = read.get_context_tags()
        file_name = context_tags.get('user_filename_input','none')
        model = context_tags.get('local_bc_temp_model','none') # info on the model used
        kit = context_tags.get('sequencing_kit','none') # sequencing kit
        print(f'sample {file_name} processed with {kit} and read by {model}')   
    
    # # Iterate through all reads in the file (supports multi-read FAST5)
    # for read_id in f5.get_read_ids():
    #     read = f5.get_read(read_id)
        
    #     # Get context_tags and tracking_id (contains chemistry info)
    #     context_tags = read.get_context_tags()
    #     tracking_id = read.get_tracking_id()
        
    #     print(f"Read ID: {read_id}")
    #     print("Context tags:")
    #     for k, v in context_tags.items():
    #         print(f"  {k}: {v}")
        
    #     print("Tracking info:")
    #     for k, v in tracking_id.items():
    #         print(f"  {k}: {v}")
        
    #     # Access raw signal if needed
    #     signal = read.get_raw_data()
    #     print(f"Signal length: {len(signal)}")
