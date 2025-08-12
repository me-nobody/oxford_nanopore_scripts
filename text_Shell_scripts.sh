# script to create csv for column 7 in full_report_bam_stats.csv of teloseq

cat file.csv |cut -d "," -f 7|sed '2,${s/|/,/g
s/:/,/
s/_/,/
s/-/,/
s/\./,/}'>col7.txt

# remove the first line
sed '1d' modified_col7_bam_stats_file.txt > modified.col7.bam.stats.file.txt

# add the header
awk 'BEGIN {print "rname,arm,ref,chrom,parent,start,end"}' '{print $0}' modified.col7.bam.stats.file.txt

# delete col 7 from original file
cut -d "," -f7 --complement test.csv > newtest.txt
# OR
# get the number of fields
awk -F , 'END{print NF}' file
# delete the column
awk -F , '{print $0}' file|cut -d "," -f1-6,8-30 >output

# insert a line in text file
sed '1i\rname,arm,ref,chrom,parent,start,end' col7.txt >updatedcol7.txt

# add data from col7 txt file as column to newtest.txt file 
paste newtest.txt col7.txt > run_EBT_1_U2OS.txt

# rsync protocol
rsync -av --ignore-existing dasaz@bluebear.bham.ac.uk:/rds/projects/b/broderra-mrc-alt-telomere/Anu/pilot_data/Nanopore_271123/dnascent/*.txt* /home/anubratadas/Documents/ALT_LAB_BHAM/err_logs/dnascent_log_files/

# awk does not match
awk  '$0 !~/#/ {print $1}' leftForks_DNAscent_forkSense_stressSignatures.bed 
