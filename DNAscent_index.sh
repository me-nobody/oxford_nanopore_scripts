# commands to check for container
apptainer run-help ./containers/DNAscent.sif
# create a directory dnascent within the data folder Anu/pilot_data/EBT_U2OS_15112023/dnascent
apptainer inspect ../../../containers/DNAscent.sif
# create an interactive environment
module load slurm-interactive
fisbatch_screen --nodes=1-1 --time=3:0:0 --mem=36G --ntasks=16
# run DNAscent index

apptainer run ../../../containers/DNAscent.sif index -f ../output_pod5s --output ebt_index.dnascent

   

