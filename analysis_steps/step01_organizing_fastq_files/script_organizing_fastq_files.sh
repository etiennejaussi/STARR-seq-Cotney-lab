##########
# This script allows to organize the fastq files
# The fastq files come from this directory : /archive/labs/Cotney/DATA/STARR-Seq/ChIPSTARRseq_test/test_run_3/TYAN01*
##########

# create a directory for the unconcatenated fastq files
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/unconcatenated_fastq

# move the unconcatenated fastq files in the fastq/unconcatenated_fastq directory
mv /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/*L00* /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/unconcatenated_fastq

# create a directory for the concatenated fastq files
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq

# move the concatenated fastq files in the fastq/concatenated_fastq directory
mv /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/*L00* /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq

# create a texte file which contains only the names of the samples
ls *.fastq.gz | sed 's/.fastq.gz//g' | cut -d "_" -f 1 | uniq > sample_list.txt
