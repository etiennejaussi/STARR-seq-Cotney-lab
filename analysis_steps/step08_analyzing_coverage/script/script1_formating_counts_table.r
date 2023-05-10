#######
# This script allows to obtain the fomrating coverage counts table.
# We keep only row for which there were data for at least two RNA replicates
#######

# load library
library(tidyverse)

# load raw data
Isoplasmid1_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/Isoplasmid1_cov.txt", header=FALSE)
Isoplasmid2_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/Isoplasmid2_cov.txt", header=FALSE)
Isoplasmid3_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/Isoplasmid3_cov.txt", header=FALSE)
mRNA1_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/mRNA1_cov.txt", header=FALSE)
mRNA2_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/mRNA2_cov.txt", header=FALSE)
mRNA3_cov <- read.table("C:/Users/etien/Desktop/Stage_Uconn_M1/experiments/CHIP_STARR_seq/analysis_steps/step08_analyzing_coverage/data/mRNA3_cov.txt", header=FALSE)

# from dataframe to tibble
Isoplasmid1_cov <- as_tibble(Isoplasmid1_cov)
Isoplasmid2_cov <- as_tibble(Isoplasmid2_cov)
Isoplasmid3_cov <- as_tibble(Isoplasmid3_cov)
mRNA1_cov <- as_tibble(mRNA1_cov)
mRNA2_cov <- as_tibble(mRNA2_cov)
mRNA3_cov <- as_tibble(mRNA3_cov)

# join all counts
peak_coverage <- Isoplasmid1_cov %>% full_join(Isoplasmid2_cov, by=c("V1","V2","V3")) %>%
  full_join(Isoplasmid3_cov, by=c("V1","V2","V3")) %>% 
  full_join(mRNA1_cov, by=c("V1","V2","V3")) %>% 
  full_join(mRNA2_cov, by=c("V1","V2","V3")) %>%
  full_join(mRNA3_cov, by=c("V1","V2","V3"))
print(peak_coverage, width = Inf)

# rename columns
colnames(peak_coverage) <- c("chro","start","end","dna1_cov","dna2_cov","dna3_cov","rna1_cov","rna2_cov","rna3_cov")
print(peak_coverage, width = Inf)

# keep rna data for which there were data for at least two replicates
peak_coverage <- peak_coverage[which(peak_coverage$rna1_cov != 0 & peak_coverage$rna2_cov != 0), ]
peak_coverage <- peak_coverage[which(peak_coverage$rna1_cov != 0 & peak_coverage$rna3_cov != 0), ]
peak_coverage <- peak_coverage[which(peak_coverage$rna2_cov != 0 & peak_coverage$rna2_cov != 0), ]
peak_coverage

# create a table with the position information
cov_pos_info <- peak_coverage %>% select(chro:end)
cov_pos_info

# create a table with the coverage counts
cov_counts <- peak_coverage %>% select(dna1_cov:rna3_cov)
cov_counts

# save peak_coverage object
save(peak_coverage, file = "output/output1_formating_counts_table/peak_coverage.Rdata")

