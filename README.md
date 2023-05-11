# Introduction

* The data comes from an experiment performed in Cotney lab 4 years ago. I have retrieved the fastq files in the following directory : `/archive/labs/Cotney/DATA/STARR-Seq/ChIPSTARRseq_test/test_run_3/TYAN01*`.
* Link to the hg19 reference genome : https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips. I downloaded the `hg19.fa.gz` file
* Link to the list of the chrommosome size : https://gist.github.com/nimezhu/d97a302816ea751e25571773d075b35e
* Link to Chip-seq data that I used for the step10 : https://cotneylab.cam.uchc.edu/~jcotney/CRANIOFACIAL_HUB/hg19/

# Description of the analysis steps

## Step01 : Organization of fastq files
* For this, I used the script named : `script_organizing_fastq_files`

## Step02 : Mapping the reads
* For this, I used the script named : `script_bowtie2_align_data`

## Step03 : Filtering BAM files
* For this, I used the script named : `script_samtools_filtering_bam`

## Step04 : Indexing BAM files
* For this, I used the script named : `script_samtools_indexing_bam`

## Step05 : Transform bam files to bed files
* For this, I used the script named : `script_bedtools_bamtobed`

## Step06 : Peak calling
* For this, I used the script named : `script_macs2_peak_calling_rna-ctr`

## Step07 : Obtain coverage
* For this, I used the script named : `script_bedtools_obtain_coverage`
