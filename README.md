# Introduction

The data comes from an experiment performed in Cotney lab 3 years ago. I have retrieved the fastq files in the following directory : `/archive/labs/Cotney/DATA/STARR-Seq/ChIPSTARRseq_test/test_run_3/TYAN01*`.

# Description of the analysis steps

## Step 01 :
* Organization of fastq files : For this, I used the script named : `script_organizing_fastq_files`

## Step 02 :
* Mapping the reads. For this, I used the script named : `script_bowtie2_align_data`

## Step 03 :
* Converting SAM files to BAM files. For this, I used the script named : `script_samtools_sam_to_bam`

## Step 04 :
* Sorting BAM files by genomic coordinates. For this, I used the script named : `script_samtools_sorting_bam_files`

## Step 05 :
* Filtering BAM files to remove unmapped reads, PCR or optical duplicate reads and low quality reads. For this, I used the script named : `script_samtools_filtering_bam_files`
