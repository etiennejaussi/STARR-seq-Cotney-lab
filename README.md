# Introduction

The data comes from an experiment performed in Cotney lab 3 years ago. I have retrieved the fastq files in the following directory : `/archive/labs/Cotney/DATA/STARR-Seq/ChIPSTARRseq_test/test_run_3/TYAN01*`.

Path to the GTF annotation file : `https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz`.
To decompress the GTF file, I used the following commande line : `gunzip Homo_sapiens.GRCh37.87.chr.gtf.gz`  
To remove Mitochondrial chromosomes, I used the following commande line : `cat Homo_sapiens.GRCh37.87.chr.gtf | awk '{print "chr"$0}' | grep -v chrG | grep -v chrH | grep -v chrMT > Homo_sapiens.GRCh37.87.nuclear_chr.gtf`

reference genome : https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

chrommosome size : https://gist.github.com/nimezhu/d97a302816ea751e25571773d075b35e

# Description of the analysis steps

## Step01 : Organization of fastq files
* For this, I used the script named :

## Step02 : Mapping the reads
* For this, I used the script named :

## Step03 : Filtering BAM files
* For this, I used the script named :

## Step04 : Indexing BAM files
* For this, I used the script named :

## Step05 : Extract fragments
* For this, I used the script named :

## Step06 : Peak calling
* For this, I used the script named :
