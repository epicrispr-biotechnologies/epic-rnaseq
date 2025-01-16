# deseq2_covariate_lrt.r
#
# This script uses DESeq2 to perform model selection
# from deseq2 models using different numbers of covariates

# Author: Tyler Borrman
# Date: 01/15/2025
#

# Load libraries
library(DESeq2)
library(tximport)
library(readr)

# Load data
metadata_path <- "/home/tylerborrman/epic-rnaseq_data/JIRA_BDS-1062_NHP_Heart_EPI321_vs_Ctrl/Input_Data/nhp_heart_rnaseq_epi321_vs_control_covariate_sex_samplesheet.csv"
metadata <- read_csv(metadata_path)
metadata <- metadata[c("sample", "condition", "sex")]
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample

# Generate DESeq2 objects
salmon_path <- "/home/tylerborrman/epic-rnaseq_data/JIRA_BDS-1062_NHP_Heart_EPI321_vs_Ctrl/salmon"

files <- file.path(
  salmon_path,
  metadata$sample,
  "quant.sf"
)
names(files) <- metadata$sample
tx2gene <- read_tsv(
  paste(salmon_path, "tx2gene.tsv", sep="/"),
  col_names = c("tx_id", "gene_id", "gene_name")
)

tx2gene <- tx2gene[c("tx_id", "gene_name")]
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# DESeq2 model with condition as only covariate
dds_unfiltered <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = ~ condition
)

# DESeq2 model with condition and sex as covariates
dds_unfiltered_sex <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = ~ sex + condition
)

# QC:filter
# Here we perform pre-filtering to keep only genes 
# that have a count of at least 10 reads summed across all samples.

# DESeq2 model with condition as only covariate
keep <- rowSums(counts(dds_unfiltered)) >= 10
dds <- dds_unfiltered[keep,]

# DESeq2 model with condition and sex as covariates
keep_sex <- rowSums(counts(dds_unfiltered_sex)) >= 10
dds_sex <- dds_unfiltered_sex[keep_sex,]

# Run DESeq2
dds <- DESeq(dds)
dds_sex <- DESeq(dds_sex)

res <- results(dds, contrast = c("condition","EPI321","CONTROL"), alpha = 0.01)
res_sex <- results(dds_sex, contrast = c("condition","EPI321","CONTROL"), alpha = 0.01)

dds_lrt <- DESeq(dds_sex, test ="LRT", reduced = ~ condition)

res_lrt <- results(dds_lrt)

print("Total genes analyzed:")
nrow(res_lrt) # 16,834
print("Number of genes with a significant LRT p-value: '~ sex + condition' vs '~ condition'")
sum(res_lrt$padj < 0.05, na.rm = TRUE) # 5,126 (30.5%)


# Let's add a random covariate to the model and see how it affects the results
metadata$random_covariate <- rep(c("A", "B"), 6)

# DESeq2 model with condition and random covariate
dds_unfiltered_rc <- DESeqDataSetFromTximport(
  txi,
  colData = metadata,
  design = ~ random_covariate + condition
)
keep_rc <- rowSums(counts(dds_unfiltered_rc)) >= 10
dds_rc <- dds_unfiltered_rc[keep_rc,]
dds_rc <- DESeq(dds_rc)
dds_lrt_rc <- DESeq(dds_rc, test ="LRT", reduced = ~ condition)
res_lrt_rc <- results(dds_lrt_rc)

print("Total genes analyzed:")
nrow(res_lrt_rc) # 16,834
print("Number of genes with a significant LRT p-value: '~ random_covariate + condition' vs '~ condition'")
sum(res_lrt_rc$padj < 0.05, na.rm = TRUE) # 0 (0%)
