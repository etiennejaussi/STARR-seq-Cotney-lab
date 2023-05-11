#!/bin/bash

#SBATCH --job-name=quantify_regulatory_regions
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=30G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o quantify_regulatory_regions_%j.out
#SBATCH -e quantify_regulatory_regions_%j.err

############################################
#This script allows to estimate the coverage
############################################

# create directory for coverage result
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/

# export the temp variable
export TMPDIR=/home/FCAM/ejaussi/tmp/

# modules
module load BEDtools/2.29.0

#####################################################################
# This part allows to keep peaks that intersect between dna1 and rna1
#####################################################################

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/dna1_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna1_peaks.bed

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA1_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_peaks.bed

# run bedtools intersect
bedtools intersect \
	-wa \
	-wb \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna1_peaks.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_peaks.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_dna1_intersect.bed

# create a file that contain the original entry in a for each overlap
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_dna1_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna1_intersect.bed

# create a file that contain the original entry in b for each overlap
cut -f4-6 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_dna1_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_intersect.bed

#####################################################################
# This part allows to keep peaks that intersect between dna2 and rna2
#####################################################################

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/dna2_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna2_peaks.bed

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA2_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_peaks.bed

# run bedtools intersect
bedtools intersect \
	-wa \
	-wb \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna2_peaks.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_peaks.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_dna2_intersect.bed

# create a file that contain the original entry in a for each overlap
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_dna2_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna2_intersect.bed

# create a file that contain the original entry in b for each overlap
cut -f4-6 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_dna2_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_intersect.bed

#####################################################################
# This part allows to keep peaks that intersect between dna3 and rna3
#####################################################################

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/dna3_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna3_peaks.bed

# keep the chrom name, start and end columns
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA3_peaks.narrowPeak \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_peaks.bed

# run bedtools intersect
bedtools intersect \
	-wa \
	-wb \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna3_peaks.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_peaks.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_dna3_intersect.bed

# create a file that contain the original entry in a for each overlap
cut -f1-3 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_dna3_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna3_intersect.bed

# create a file that contain the original entry in b for each overlap
cut -f4-6 /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_dna3_intersect.bed \
				> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_intersect.bed

################################################################
# This part allows to obtain a file that contain all peaks areas
################################################################

# combine the files
cat /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna1_intersect.bed \
	/home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_intersect.bed \
	/home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna2_intersect.bed \
	/home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_intersect.bed \
	/home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/dna3_intersect.bed \
	/home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_intersect.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_peaks.bed

# sort by chromosome name and position then merge peaks
sort -k1,1 -k2,2n /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_peaks.bed \
	| bedtools merge -i - \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed

######################################
# This part allows to compute coverage
######################################

# run bedtools coverage for DNA1
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid1.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/Isoplasmid1_cov.txt

# run bedtools coverage for RNA1
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA1.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA1_cov.txt

# run bedtools coverage for DNA2
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid2.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/Isoplasmid2_cov.txt

# run bedtools coverage for RNA2
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA2.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA2_cov.txt

# run bedtools coverage for DNA3
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/Isoplasmid3.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/Isoplasmid3_cov.txt

# run bedtools coverage for RNA3
bedtools coverage \
	-counts \
	-a /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA_dna_merged.bed \
	-b /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bed/mRNA3.bed \
	> /home/FCAM/ejaussi/CHIP_STARR_seq/results/coverage/mRNA3_cov.txt
