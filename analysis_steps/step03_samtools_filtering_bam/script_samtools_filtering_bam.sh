#!/bin/bash

#SBATCH --job-name=filterbam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=32G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o filterbam_%j.out
#SBATCH -e filterbam_%j.err

##########
# This script allows to filter bam files
# www.biostars.org/p/95929/ : info for the filtering step on bowtie2 output
##########

# modules
module load samtools/1.9

# run samtools view for filtering
## -f 0x2 : get only "properly paired" alignments
## "XS:i:" : tag reads that have another "valid" mapping
## grep -v : remove reads that have another "valid" mapping

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_filtered.bam

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_filtered.bam

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_filtered.bam

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_filtered.bam

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_filtered.bam

samtools view -h -f 0x2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_sorted.bam \
	| grep -v "XS:i:" \
	| samtools view -h -S -b > /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_filtered.bam
