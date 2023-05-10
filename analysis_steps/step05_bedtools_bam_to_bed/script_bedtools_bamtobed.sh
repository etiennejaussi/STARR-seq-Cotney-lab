#!/bin/bash

#SBATCH --job-name=bamtobed
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=30G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o bamtobed_%j.out
#SBATCH -e bamtobed_%j.err

##########
# This script allows to transform .bam files to .bed files
##########

# create a directory for bed files
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed

# modules
module load BEDtools/2.29.0

# run bedtools
# bed files are already sorted by chro number (chr1, chr2, etc...)
bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid1.bed

bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid2.bed

bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid3.bed

bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA1.bed

bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA2.bed

bedtools bamtobed -i /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam \
	| cut -f1-3 \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA3.bed