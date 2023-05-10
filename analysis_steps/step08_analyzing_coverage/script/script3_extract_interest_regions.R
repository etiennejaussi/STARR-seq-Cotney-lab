#######
# This script allows to know the nomber of false postive results and filter the data
# to have the the most promising regions
#######

# load library
library(tidyverse)

# load data
load("output/2_differential_expression/all_cov_and_deseq_results.csv")

# count the nomber of peaks with an adjusted pvalue lower than 0.05
sum(cov_deseq$padj < 0.05, na.rm=TRUE)

# plot the distribution of the pvalues
hist(cov_deseq$pvalue)
abline(h=5200, col="red", lwd=3)

# nomber of pvalue above the red line
(length(which(cov_deseq$pvalue <= 0.05)) - 5200)

# nomber of false positives modified p-values (FDR = 0.1)
(length(which(cov_deseq$pvalue <= 0.05)) - 5200) * 0.1

# nomber of true positives modified p-values (FDR = 0.1)
(length(which(cov_deseq$pvalue <= 0.05)) - 5200) * 0.9

# order by log2FC
print(cov_deseq %>% arrange(desc(log2FoldChange)), width = Inf)

# order by adjusted pvalues
print(cov_deseq %>% arrange(padj), width = Inf)

# create a table that contain all regions of interest
interest_regions_all <- cov_deseq %>% select(chro, start, end)

# filter to obtain a table that contain data with an adjusted p-value < 0.05
cov_deseq_sig <- cov_deseq %>% filter(padj < 0.05)
print(cov_deseq_sig, width = Inf)

# order by adjusted pvalues
print(cov_deseq_sig %>% arrange(padj), width = Inf)

# create a table that contain significant regions of interest
interest_regions_sig <- cov_deseq_sig %>% select(chro, start, end)

# save the the tables that contain all and significant regions of interest in a text format
write.table(interest_regions_all, file="output/3_extract_interest_regions/all_interest_regions.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(interest_regions_sig, file="output/3_extract_interest_regions/sig_interest_regions.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

