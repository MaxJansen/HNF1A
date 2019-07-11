###    Main description           ###
###    This script makes a large table of single cell sequencing read counts and associated numbers.
###    The final table should contain sample names, read file names, and a varierty of read counts.
###    There is an example sent to Max by Agata called: "31.05.2018.HNF1A_SS2_scRNA.map_qc_stats.txt"
###    End of main description    ###

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

###    Part 1. Raw reads    ###
###    This is one backbone of the table. It will start with sequence names,
###    and raw read counts. ###

# Import reads from both directions, change working directory accordingly
setwd("~/Oxford 2.0/HNF1A/Tables/raw_read_count")
raw_count_files <- list.files(path = " .", pattern = "Batch")
for (i in 1:length(raw_count_files))
  assign(raw_count_files[i],
         read.csv(raw_count_files[i], header = FALSE))

# Already Checked if all paired Batch files have equal number of lines using wc -l in bash
# Check if fastq _1 and _2 reads have equal number of raw read counts

# Function to compare
compare_reads1_2 <- function(Batchread1, Batchread2) {
  for (i in seq(2, dim(Batchread1)[1], 2)) {
    if (Batchread1[i, ] != Batchread2[i, ]) {
      print("inconsistent reads")
    }
  }
}

# Loop through read directions 1 and 2 to see if they have equal raw read counts
# Do a bit of cleaning to get Batch numbers without read direction
Batches <- ls()
Batches <- Batches[grep(".txt", Batches)]
Batches <- Batches[grep("Batch", Batches)]
Batches_1 <- Batches[grep("_1.txt", Batches)]
Batches_2 <- Batches[grep("_2.txt", Batches)]
Batches <- gsub("_1.txt", "", Batches)
Batches <- gsub("_2.txt", "", Batches)
Batches <- unique(Batches)

# Now you have unique Batch numbers, loop through paired Batch read directions and compare
for (i in Batches_1) {
  Batchreadname1 <- get(i)
  Batchreadname2 <- get(gsub("_1.txt", "_2.txt", i))
  print(i)
  print(gsub("_1.txt", "_2.txt", i))
  compare_reads1_2(Batchreadname1, Batchreadname2)
}

# Equal read numbers for both directions in all batches, so:
# Keep numerical values from reads 1 only, move values from these rows to new column
# Divide these lines counts by 4 for actual raw_reads column number

# This function makes clean tables and performs calculation
rawread_2_column <- function(Batch_info, Batch_name) {
  read_names <- Batch_info[seq(1, nrow(Batch_info), 2), ]
  read_names <- gsub("./", "", read_names)
  read_names <- gsub("_1.fastq.gz", "", read_names)
  raw_reads <- Batch_info[seq(2, nrow(Batch_info), 2), ]
  Batch_new <- data.frame(read_names, raw_reads)
  colnames(Batch_new) <- c("ReadName", "RawReads")
  Batch_new$RawReads <- as.numeric(as.character(Batch_new$RawReads))
  Batch_new$RawReads <-
    lapply(Batch_new$RawReads, function(x)
      x / 4)
  Batch_new$Batchname <- Batch_name
  return(Batch_new)
}

All_batches <- data.frame()
for (i in Batches) {
  i <- rawread_2_column(get(paste(i, "_1.txt", sep = "")), i)
  All_batches <- rbind(All_batches, i)
}

#Save All_batches.
All_batches <- as.data.frame(All_batches)
write.csv(All_batches, file = "All_batches_raw_reads.txt", row.names = TRUE)

###    Part 1. Finished.    ###

###    Part 2. Get RG_info.    ###
###    This is the other backbone of the table. RG_info files allow linking between
###    sequence name, sample name, and all other read count columns.

# Import RG_info tables
setwd("~/Oxford 2.0/HNF1A/Tables")
RG_info1 <- read.table("RG_info1")
RG_info2 <- read.table("RG_info2a")
RG_info2second <- read.table( ("RG_info2b"))
RG_info3 <- read.table("RG_info3")
RG_info4 <- read.table("RG_info4")

RG_info1$Batch <- as.factor("Batch1")
RG_info2$Batch <- as.factor("Batch2")
RG_info2second$Batch <- as.factor("Batch2second")
RG_info3$Batch <- as.factor("Batch3")
RG_info4$Batch <- as.factor("Batch4")

# Make a list of df's to edit all at once later
my.list <-
  list(RG_info1, RG_info2, RG_info2second, RG_info3, RG_info4)

# slim down the tables, keep name, date and batch. Bind df's and reformat column names.
select_column <- function(y) {
  z <- subset(y, select = c(V1, V3, V7, Batch))
  return(z)
}

my.list <- lapply(my.list, select_column)
RG_infoAll <- do.call(rbind, my.list)

# Tested if all values in RG_infoAll$V1 are unique, they are, so can be mapped directly to All_batches
# Link read filename to sample name
Mastertable <-
  merge(RG_infoAll, All_batches, by.x = "V1", by.y = "ReadName")
Mastertable <- Mastertable[, -4]
Mastertable <- Mastertable[, c(2, 5, 1, 3, 4)]
colnames(Mastertable) <-
  c("Name", "BatchName", "ReadNumber", "Date", "RawReads")
Mastertable$Name <- gsub('SM="', "", Mastertable$Name)
Mastertable$Name <- gsub('"', "", Mastertable$Name)

#Test to see if unique within batches. Answer: No.
#They are unique between batches, see quick_check_unique.R script for proof!

###    Part 2. Finished.    ###


###    Part 3. Add quantitative info from summary tables.    ###
###    This is some of the required quantitative data. (TPM and Trimmed not included)

# Import the quantitative data separately. This is a function that does that:
create_df_list <- function(directory) {
  setwd(directory)
  temp <- list.files(pattern = "*.summary")
  batch_quant_list <- lapply(setNames(temp, make.names(
    gsub("*.gene.counts.summary", "", temp)
  )),
  function(x)
    read.csv(x, stringsAsFactors = FALSE, header = F)[-1, ])
  return(batch_quant_list)
}

# Apply function to directories containing Batches
Batch1and2_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts/gene_counts_merged")
Batch2second_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts/gene_counts_Batch2/")
Batch3_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts/gene_counts_Batch3/")
Batch4_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts/gene_counts_Batch4/")

# The 4 lists of dataframes must be formatted for merging.
# Function splits columns and transposes df's:
format_quant_list <- function(Batch_list) {
  batch_qdf <- t(as.data.frame(Batch_list))
  new_batch_qdf <- as.data.frame(rownames(batch_qdf))
  for (i in 1:ncol(batch_qdf)) {
    new_cols <- str_split_fixed(batch_qdf[, i], "\t", 2)
    new_batch_qdf <- cbind(new_batch_qdf, new_cols)
  }
  rownames(new_batch_qdf) <- new_batch_qdf[, 1]
  new_batch_qdf <- new_batch_qdf[, -1]
  new_colnames <- new_batch_qdf[1, seq(1, ncol(new_batch_qdf), 2)]
  new_colnames <- lapply(new_colnames, as.character)
  new_batch_qdf <- new_batch_qdf[, -seq(1, ncol(new_batch_qdf), 2)]
  colnames(new_batch_qdf) <- new_colnames
  new_batch_qdf <- as.data.frame(new_batch_qdf)
  new_batch_qdf[] <-
    lapply(new_batch_qdf, function(x)
      as.numeric(as.character(x)))
  return(new_batch_qdf)
}

#Do this for all batch lists of df's:
Batch1and2_qdf <- format_quant_list(Batch1and2_df_list)
Batch2second_qdf <- format_quant_list(Batch2second_df_list)
Batch3_qdf <- format_quant_list(Batch3_df_list)
Batch4_qdf <- format_quant_list(Batch4_df_list)

# !!! Warning. Some notes on Batch2 and duplicates:
# Batch2 and Batch2second are duplicates, sum these values and then use (eventually) unique 'Name' column,
# not the Readnumber column. To do this, convert the batch named Batch2second to Batch2, then aggregate by
# Name + BatchName, to keep the BatchName in the Fulltable.
# You can do this, because in later steps Batch1and2_qdf will have summed values for the old and new
# Batch2 summary tables.
# Batch3 and Batch4 also contain duplicates. Do the same aggregation. These will have the same batchNames
# for duplicates though.

#Rename Batch2second to Batch2 in Mastertable:
Mastertable$BatchName[Mastertable$BatchName == "Batch2second"] <-
  "Batch2"

#First aggregate RawReads(numeric) for duplicated(Name) rows, then merge.
Mastertable$Name <- gsub("-", ".", Mastertable$Name)
Mastertable$RawReads <- as.numeric(Mastertable$RawReads)
Mastertable_agg <-
  aggregate(RawReads ~ Name + BatchName, data = Mastertable, FUN = sum)

#Merge for each separate batch.
Fullbatch1 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch1", ],
    Batch1and2_qdf,
    by.x = "Name",
    by.y = 0,
    all.x = TRUE
  )
Fullbatch2 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch2", ],
    Batch2second_qdf,
    by.x = "Name",
    by.y = 0,
    all.x = TRUE
  )
Fullbatch3 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch3", ],
    Batch3_qdf,
    by.x = "Name",
    by.y = 0,
    all.x = TRUE
  )
Fullbatch4 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch4", ],
    Batch4_qdf,
    by.x = "Name",
    by.y = 0,
    all.x = TRUE
  )

#Make Fulltable, check it, and then clean it.
Fulltable <- rbind(Fullbatch1, Fullbatch2, Fullbatch3, Fullbatch4)

#Check : Is RawReads about the same as rowsum?
Fulltable$RawSum <- rowSums(Fulltable[, 4:14])
ratio_df <-
  as.data.frame(cbind(Fulltable$RawReads, Fulltable$RawSum))
colnames(ratio_df) <- c("RawReads", "RawSum")
ratio_df$RawRatio <- ratio_df[, 1] / ratio_df[, 2]

#Find extreme rows (somewhat arbitrary numbers, but indicative):
ratio_df[(ratio_df$RawRatio > 1.01), ]
ratio_df[(ratio_df$RawRatio < 1), ]

#Summary check passed. Ratios near 1.
#Remove Irrelevant columns
Fulltable$Unassigned_Ambiguity <- NULL
Fulltable$Unassigned_MappingQuality <- NULL
Fulltable$Unassigned_FragmentLength <- NULL
Fulltable$Unassigned_Chimera <- NULL
Fulltable$Unassigned_Secondary <- NULL
Fulltable$Unassigned_Nonjunction <- NULL
Fulltable$Unassigned_Duplicate <- NULL

#Create 'pc' columns
Fulltable$Raw_reads_pc <- Fulltable[, 3] / Fulltable[, 8]
Fulltable$Assigned_pc <- Fulltable[, 4] / Fulltable[, 8]
Fulltable$MultiMapping_pc <- Fulltable[, 5] / Fulltable[, 8]
Fulltable$NoFeatures_pc <- Fulltable[, 6] / Fulltable[, 8]
Fulltable$Unmapped_pc <- Fulltable[, 7] / Fulltable[, 8]

###    Part 3. Finished.    ###

###    Part 4. Add quantitative info from genes TPM and gene.counts.    ###
###    For the gene.counts: Use a separate script on the cluster: single_cell_gene_count.R ###

# Gene TPM. The files are too large to import locally, use 2.count_genes_elder.R script on cluster
# '*gene.counts.tsv' files are also too large, use 2.calc_reads_pc.R to condense it to the required format

#Import results locally:
setwd("~/Oxford 2.0/HNF1A/Tables/table2/")
raw_count_files <- list.files(path = ".", pattern = ".txt")
for (i in 1:length(raw_count_files))
  assign(raw_count_files[i],
         read.csv(raw_count_files[i], header =
                    FALSE))

#Make a list of df's to edit all at once later
tpm.list <-
  list(tpm_batch1_and_2.txt,
       tpm_batch2.txt,
       tpm_batch3.txt,
       tpm_batch4.txt)
tpm_df <- do.call(rbind, tpm.list)
tpm_df <-
  as.data.frame(unlist(str_split_fixed(tpm_df[, 1], "\t", 2)))
colnames(tpm_df) <- c("Name", "Genes1tpm")
tpm_df$Genes1tpm <- as.numeric(as.character(tpm_df$Genes1tpm))
tpm_df <- aggregate(Genes1tpm ~ Name, data = tpm_df, FUN = sum)
tpm_df$Name <- as.character(tpm_df$Name)
gene_counts_list <-
  list(
    gene_counts_batch1_and_2.txt,
    gene_counts_batch2.txt,
    gene_counts_batch3.txt,
    gene_counts_batch4.txt
  )
gene_counts_df <- do.call(rbind, gene_counts_list)
gene_counts_df <-
  as.data.frame(unlist(str_split_fixed(gene_counts_df[, 1], "\t", 3)[, 2:3]))
colnames(gene_counts_df) <- c("Name", "top100_reads_pc")
gene_counts_df$top100_reads_pc <-
  as.numeric(as.character(gene_counts_df$top100_reads_pc))
gene_counts_df <-
  aggregate(top100_reads_pc ~ Name, data = gene_counts_df, FUN = mean)

Fullerbatch <-
  merge(Fulltable,
        tpm_df,
        by.x = "Name",
        by.y = "Name",
        all.x = TRUE)
Fullertable <-
  merge(
    Fullerbatch,
    gene_counts_df,
    by.x = "Name",
    by.y = "Name",
    all.x = TRUE
  )

#Make a subset for top100_pc comparison and a subset for 1tpm comparison
Fullertable_ss_select <-
  Fullertable[grep("singlecell", Fullertable$Name), ]
p <-
  ggplot(Fullertable_ss_select,
         aes(BatchName, Genes1tpm, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fullertable_ss_select,
         aes(BatchName, top100_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()


#Some quick plots [Save for later reference, but move eventually]
p <- ggplot(Fulltable, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, Assigned, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable,
         aes(BatchName, Unassigned_MultiMapping, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable,
         aes(BatchName, Unassigned_NoFeatures, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable,
         aes(BatchName, Unassigned_Unmapped, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, RawSum, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, Assigned_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, MultiMapping_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, NoFeatures_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, Unmapped_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, Raw_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

#This is a test
