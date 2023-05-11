#!/bin/bash

#SBATCH --job-name=macs2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=30G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o macs2_%j.out
#SBATCH -e macs2_%j.err

##########
# This script allows to identify areas in the genome that have been enriched with aligned reads with macs2
##########

# create a directory for the macs output files
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/dna_as_control
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/rna_as_control


# modules
module load macs2/2.1.2

#####################
# Narrow peak
#####################

#####
# DNA as control
#####
# run macs2 for sample 1
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA1 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA1.log

# run macs2 for sample 2
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA2 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA2.log

# run macs2 for sample 3
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA3 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna_as_control/mRNA3.log

#####
# RNA as control
#####
# run macs2 for sample 1
macs2 callpeak \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA1 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/mRNA1.log

# run macs2 for sample 2
macs2 callpeak \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA2 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/mRNA2.log

# run macs2 for sample 3
macs2 callpeak \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA3 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/rna_as_control/mRNA3.log

#####
# DNA
#####

# run macs2 for sample 1
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	-g 2.7e9 \
	-n dna1 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna1.log

# run macs2 for sample 2
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	-g 2.7e9 \
	-n dna2 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna2.log

# run macs2 for sample 3
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	-g 2.7e9 \
	-n dna3 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/narrow_peak/dna3.log

#####################
# Broad peak
#####################

#####
# RNA
#####

# run macs2 for sample 1
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	--broad-cutoff \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA1 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/mRNA1.log

# run macs2 for sample 2
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	--broad-cutoff \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA2 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/mRNA2.log

# run macs2 for sample 3
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam \
	-c /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	--broad-cutoff \
	--pvalue 5e-3 \
	-g 2.7e9 \
	-n mRNA3 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/mRNA3.log

#####
# DNA
#####

# run macs2 for sample 1
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	-g 2.7e9 \
	-n dna1 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/dna1.log

# run macs2 for sample 2
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	-g 2.7e9 \
	-n dna2 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/dna2.log

# run macs2 for sample 3
macs2 callpeak \
	-t /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam \
	--extsize 300 \
	--keep-dup all \
	--nomodel \
	--broad \
	-g 2.7e9 \
	-n dna3 \
	--outdir /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/ \
	2> /home/FCAM/ejaussi/CHIP_STARR_seq/results/macs2/broad_peak/dna3.log