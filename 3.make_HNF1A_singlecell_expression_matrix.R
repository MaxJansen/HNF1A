### Make the raw-single-cell expression matrix ###
### Runs with R3.4.3 on /apps/well
# This script will make a large table of expression of all single cells (columns)
# With a gene per row
# Although some cells have been sampled twice, the .tsv files used to make this table
# contain exactly one cell per column, creating exactly 1920 unique columns

# Load libraries
library(dplyr)

# Go to correct directories
setwd("/well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_scRNA_Max")

# Read all .tsv files fpr the 4 batches
Batch1 <-
  read.table(
    "Batch1/23.07.2019.Batch1.gene.counts.tsv",
    header = T,
    sep = "\t",
    row.names = 1
  )

Batch2 <-
  read.table(
    "Batch2/23.07.2019.Batch2.gene.counts.tsv",
    header = T,
    sep = "\t",
    row.names = 1
  )

Batch3 <-
  read.table(
    "Batch3/23.07.2019.Batch3.gene.counts.tsv",
    header = T,
    sep = "\t",
    row.names = 1
  )

Batch4 <-
  read.table(
    "Batch4/23.07.2019.Batch4.gene.counts.tsv",
    header = T,
    sep = "\t",
    row.names = 1
  )


# Remove empty 'GeneName' column in all batches except Batch1 (you only need one in the final table)
Batch2 <- subset(Batch2, select = -GeneName)
Batch3 <- subset(Batch3, select = -GeneName)
Batch4 <- subset(Batch4, select = -GeneName)

# colbind these df's
l <- list(Batch1, Batch2, Batch3, Batch4)
all_raw_geneexpr <- do.call(cbind, unname(l))
#all_raw_geneexpr <- bind_cols(Batch1, Batch2, Batch3, Batch4)
# This is the raw gene expression matrix

# Select singlecell columns: use colindex
sc_cols_index <- grep("singlecell", colnames(all_raw_geneexpr))
# Use index to select, don't forget the GeneName col!
sc_raw_geneexpr <- all_raw_geneexpr[,c(1,sc_cols_index)]

# Dimension check before writing:
#> dim(all_raw_geneexpr)
#[1] 63677  1921
#> dim(sc_raw_geneexpr)
#[1] 63677  1825
# Note: The actual number of cells are 1920 and 1824 (remember GeneName column!)

# Write tables:
write.table(all_raw_geneexpr, file = "all_raw_geneexpr.txt", row.names = TRUE, col.names = TRUE)
write.table(sc_raw_geneexpr, file = "sc_raw_geneexpr.txt", row.names = TRUE, col.names = TRUE)
