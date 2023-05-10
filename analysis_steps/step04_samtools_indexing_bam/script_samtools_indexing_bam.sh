#!/bin/bash

#SBATCH --job-name=indexbam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=32G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o indexbam_%j.out
#SBATCH -e indexbam_%j.err

##########
# This script allows to index bam files
##########

# modules
module load samtools/1.9

# run samtools sort
## -b : create a BAI index
## -@ : nomber of threads

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam

samtools index -b -@ 6 \
	/home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam