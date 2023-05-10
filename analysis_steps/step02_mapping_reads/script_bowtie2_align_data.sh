#!/bin/bash

#SBATCH --job-name=bowtie2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=64G
#SBATCH --mail-user=jaussi@uchc.edu
#SBATCH -o bowtie2_%j.out
#SBATCH -e bowtie2_%j.err

##########
# This script allows the mapping of the reads on the hg19 reference genome,
# convert output sam files to bam files and sort bam files
##########

# create a directory for the output sam files
mkdir /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam

# modules
module load bowtie2/2.2.9
module load samtools/1.9

# run bowtie2 for each sample with a for loop
## -p : nomber of threads
## -x : path to reference genome
## -1 : fastq files corresponding to the reads 1
## -2 : fastq files corresponding to the reads 2

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid1_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid1_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid1_sorted.bam

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid2_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid2_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid2_sorted.bam

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid3_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/Isoplasmid3_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/Isoplasmid3_sorted.bam

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA1_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA1_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA1_sorted.bam

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA2_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA2_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA2_sorted.bam

bowtie2 -p 12 \
	-x /isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_bowtie2/Human_genome \
	-1 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA3_R1.fastq.gz \
	-2 /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/concatenated_fastq/mRNA3_R2.fastq.gz \
	| samtools view -b -S -h - \
	| samtools sort -O bam -m 500M -@ 6 - \
> /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/bam/mRNA3_sorted.bam

#### script with for loop
# cat /home/FCAM/ejaussi/CHIP_STARR_seq/raw_data/fastq/sample_list.txt | while read sample
# do
# bowtie2 script (Replace name of sample by ${sample})
# done
####
