#######
# This script allows to obtain the differential expression information thanks to Deseq2
# and combine the results with coverage counts
#######

# load library
library("DESeq2")
library("tidyverse")

# load data
load("output/output1_formating_counts_table/peak_coverage.Rdata")
peak_coverage

# ceate an individual name for each rows
peak_rank <- paste("rank", 1:124228, sep = "")

# add this individual name to the peak_coverage object
peak_coverage <- peak_coverage %>% mutate(peak_rank = peak_rank, .before = chro)
peak_coverage

# keep only the peak_rank and coverage counts columns
countData <- peak_coverage %>% select(-chro, -start, -end)
countData

# put the peak_rank column to raw names
countData <- countData %>% column_to_rownames(var = "peak_rank")
str(countData)

# transform to matrix
countData <- as.matrix(countData)
str(countData)

# create a table with sample information
coldata <- data.frame(
  sample = c( "dna1_cov", "dna2_cov", "dna3_cov", "rna1_cov", "rna2_cov", "rna3_cov" ),
  condition = c( "dna", "dna",  "dna", "rna", "rna", "rna" ),
  row.names = "sample" )

# tranform the column condition to factor
coldata$condition <- as.factor(coldata$condition)

# construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)
dds

# set control condition as reference
dds$condition <- relevel(dds$condition, ref = "dna")

# compute differential expression
dds <- DESeq(dds)

# get gene expression table
res <- results(dds, contrast = c("condition", "rna", "dna"), name="condition_rna_vs_dna")
res

# MA plot
plotMA(res, ylim = c(-10, 10))

# transform row names to culumn and change to tibble
res <- as.data.frame(res) %>% rownames_to_column(var = "peak_rank") %>% as_tibble()
res

# combine peak_coverage table with the DEseq results
cov_deseq <- peak_coverage %>% full_join(res, by = "peak_rank")
print(cov_deseq, width = Inf)

# save the result table
write.csv(cov_deseq, file = "output/output2_differential_expression/all_cov_and_deseq_results.csv")
