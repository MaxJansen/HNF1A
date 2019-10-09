### This script selects cells/columns from the raw expression file 
### It uses a list of cells that passed thresholds I check locally
### Use R3.4.3 from /apps/well on Elder

### The thresholds are:
# Genes1tpm remove bottom 0.03 quantile 
# remove over 80% mapping to top 100
# remove more than 90% Unmapped
# At least 10% Assigned
# TrimReads > quantile 0.02

setwd("/well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_scRNA_Max")

sc_full_table <- read.table("sc_raw_geneexpr.txt", header = T)
filtered_sc <-
  read.csv(
    "/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/scripts/make_filtered_table/filtered_sc_names.txt"
  )

# Use the names of cells that passed the thresholds, then select. Dont' forget the GeneName col!
filtered_sc <- filtered_sc[,2]
filtered_names <- colnames(sc_full_table)[colnames(sc_full_table) %in% filtered_sc]
filtered_table <- sc_full_table[,c("GeneName",filtered_names)]
