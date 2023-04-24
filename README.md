# Introduction

The data comes from an experiment performed in Cotney lab 3 years ago. I have retrieved the fastq files in the following directory : `/archive/labs/Cotney/DATA/STARR-Seq/ChIPSTARRseq_test/test_run_3/TYAN01*`.

Path to the GTF annotation file : `https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz`.
To decompress the GTF file, I used the following commande line : `gunzip Homo_sapiens.GRCh37.87.chr.gtf.gz`  
To remove Mitochondrial chromosomes, I used the following commande line :  
`cat Homo_sapiens.GRCh37.87.chr.gtf | awk '{print "chr"$0}' | grep -v chrG | grep -v chrH | grep -v chrMT > Homo_sapiens.GRCh37.87.nuclear_chr.gtf`

# Description of the analysis steps

## Step01 :
* Organization of fastq files : For this, I used the script named : `script_organizing_fastq_files`

## Step02 :
* Mapping the reads. For this, I used the script named : `script_bowtie2_align_data`

## Step03 :
* Converting SAM files to BAM files. For this, I used the script named : `script_samtools_sam_to_bam`

## Step04 :
* Sorting BAM files by genomic coordinates. For this, I used the script named : `script_samtools_sorting_bam`

## Step05 :
* Indexing BAM files. For this, I used the script named : `script_samtools_indexing_bam`

## Step06 :
* Filtering BAM files to remove unmapped reads, PCR or optical duplicate reads and low quality reads. For this, I used the script named : `script_samtools_filtering_bam`

## Step07 :
* Peak calling. For this, I used the script named : `script_macs2_peak_calling`

## Step08 :
* Generating bigwig files. For this, I used the script named : `script_bamcoverage_generate_bigwig`
