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

All_batches$RawReads <- as.numeric(All_batches$RawReads)


# Do this for trimmed reads as well
# Import reads from both directions, change working directory accordingly
setwd("~/Oxford 2.0/HNF1A/Tables/trimmed_read_count/")
trim_count_files <- list.files(path = " .", pattern = "Batch")
for (i in 1:length(trim_count_files))
  assign(trim_count_files[i],
         read.csv(trim_count_files[i], header = FALSE))

# Loop through read directions 1 and 2 to see if they have equal raw read counts
# Do a bit of cleaning to get Batch numbers without read direction
Batches <- ls()
Batches <- Batches[grep(".txt", Batches)]
Batches <- Batches[grep("trim", Batches)]
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
# Divide these lines counts by 4 for actual trim_reads column number

# This function makes clean tables and performs calculation
trimread_2_column <- function(Batch_info, Batch_name) {
  read_names <- Batch_info[seq(1, nrow(Batch_info), 2), ]
  read_names <- gsub("./", "", read_names)
  read_names <- gsub("_1.fastq.gz", "", read_names)
  trim_reads <- Batch_info[seq(2, nrow(Batch_info), 2), ]
  Batch_new <- data.frame(read_names, trim_reads)
  colnames(Batch_new) <- c("ReadName", "TrimReads")
  Batch_new$TrimReads <- as.numeric(as.character(Batch_new$TrimReads))
  Batch_new$TrimReads <-
    lapply(Batch_new$TrimReads, function(x)
      x / 4)
  Batch_new$Batchname <- Batch_name
  return(Batch_new)
}

All_trim <- data.frame()
for (i in Batches) {
  i <- trimread_2_column(get(paste(i, "_1.txt", sep = "")), i)
  All_trim <- rbind(All_trim, i)
}

#Merge these tables
All_rawtrim <- 
  merge(All_batches, All_trim, by.x = "ReadName", by.y = "ReadName")
All_rawtrim <- All_rawtrim[, c(1, 2, 4, 3)]
colnames(All_rawtrim)[4] <- "Batchname"
#Save All_batches and All_trim
setwd("~/Oxford 2.0/HNF1A/Tables")
All_batches <- as.data.frame(All_batches)
write.csv(apply(All_batches, 2, as.character),
          file = "All_batches_raw_reads.txt",
          row.names = TRUE)

All_trim <- as.data.frame(All_trim)
write.csv(apply(All_trim, 2, as.character),
          file = "All_batches_trim_reads.txt",
          row.names = TRUE)

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
  merge(RG_infoAll, All_rawtrim, by.x = "V1", by.y = "ReadName")
Mastertable <- Mastertable[, -7]
Mastertable <- Mastertable[, c(2, 4, 1, 3, 5, 6)]
colnames(Mastertable) <-
  c("Name", "BatchName", "ReadNumber", "Date", "RawReads", "TrimReads")
Mastertable$Name <- gsub('SM="', "", Mastertable$Name)
Mastertable$Name <- gsub('"', "", Mastertable$Name)

# Test to see if unique within batches. Answer: No.
# They are unique between batches, see Legacy/quick_check_unique.R script for proof!

# A table of RawReads and Trimmed reads before Batch2 and Batch2second are aggregated:
Mastertable$TrimReads <- as.numeric(Mastertable$TrimReads)
Mastertable$TrimRaw_pc <- Mastertable[, 6] / Mastertable[, 5]


p <- ggplot(Mastertable, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() 

p <-
  ggplot(Mastertable, aes(BatchName, TrimReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Mastertable, aes(BatchName, TrimRaw_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()
###    Part 2. Finished.    ###


###    Part 3. Add quantitative info from summary tables.    ###
###    This is some of the required quantitative data. (TPM and not included)

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
Batch1_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts3/gene_counts_Batch1/")
Batch2_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts3/gene_counts_Batch2/")
Batch3_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts3/gene_counts_Batch3/")
Batch4_df_list <-
  create_df_list("~/Oxford 2.0/HNF1A/Tables/merged_counts3/gene_counts_Batch4/")

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
Batch1_qdf <- format_quant_list(Batch1_df_list)
Batch2_qdf <- format_quant_list(Batch2_df_list)
Batch3_qdf <- format_quant_list(Batch3_df_list)
Batch4_qdf <- format_quant_list(Batch4_df_list)

# !!! Warning. Some notes on Batch2 and duplicates:
# Batch2 and Batch2second are duplicates, sum their Raw and Trimmed values and then use (eventually) unique 'Name'
# column, not the Readnumber column. To do this, convert the batch named Batch2second to Batch2, then aggregate by
# Name + BatchName, to keep the BatchName in the Fulltable.
# You can do this, because Batch2_qdf will have summed values for the old and new
# Batch2 summary tables.
# Batch3 and Batch4 also contain duplicates. Do the same aggregation. These will have the same batchNames
# for duplicates though.

#Rename Batch2second to Batch2 in Mastertable:
Mastertable$BatchName[Mastertable$BatchName == "Batch2second"] <-
  "Batch2"

# First aggregate RawReads(numeric) for duplicated(Name) rows, merge later.
Mastertable$Name <- gsub("-", ".", Mastertable$Name)
Mastertable$RawReads <- as.numeric(Mastertable$RawReads)
Mastertable_agg <-
  aggregate(cbind(RawReads, TrimReads) ~ Name + BatchName, data = Mastertable, FUN = sum)

# Plot aggregated RawReads:
p <-
  ggplot(Mastertable_agg, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Reads") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)  

# Select controls and SC to see batch differences
Master_agg_control_select <-
  Mastertable_agg[grep("singlecell", Mastertable_agg$Name, invert = TRUE ), ]

p <- ggplot(Master_agg_control_select, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Reads") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)  

Master_agg_ss_select <-
  Mastertable_agg[grep("singlecell", Mastertable_agg$Name), ]

p <- ggplot(Master_agg_ss_select, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Reads") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)  



#Merge for each separate batch.
Fullbatch1 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch1", ],
    Batch1_qdf,
    by.x = "Name",
    by.y = 0,
    all.x = TRUE
  )
Fullbatch2 <-
  merge(
    Mastertable_agg[Mastertable_agg$BatchName == "Batch2", ],
    Batch2_qdf,
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

#Check : Is TrimReads equal to rowsum?
Fulltable$RawSum <- rowSums(Fulltable[, 5:15])
ratio_df <-
  as.data.frame(cbind(Fulltable$TrimReads, Fulltable$RawSum))
colnames(ratio_df) <- c("TrimReads", "RawSum")
ratio_df$TrimRatio <- ratio_df[, 2] / ratio_df[, 1]

#Find extreme rows (somewhat arbitrary numbers, but indicative):
ratio_df[(ratio_df$RawRatio > 1.01), ]
ratio_df[(ratio_df$RawRatio < 1), ]

#Summary check passed. Ratios near 1.
#Remove Irrelevant columns
Fulltable$Unassigned_MappingQuality <- NULL
Fulltable$Unassigned_FragmentLength <- NULL
Fulltable$Unassigned_Secondary <- NULL
Fulltable$Unassigned_Nonjunction <- NULL
Fulltable$Unassigned_Duplicate <- NULL

#Create 'pc' columns
Fulltable$Raw_sum_pc <- Fulltable[, 11] / Fulltable[, 4]
Fulltable$Assigned_pc <- Fulltable[, 5] / Fulltable[, 4]
Fulltable$Ambiguity_pc <- Fulltable[, 6] / Fulltable[, 4]
Fulltable$Multimapping_pc <- Fulltable[, 7] / Fulltable[, 4]
Fulltable$NoFeatures_pc <- Fulltable[, 8] / Fulltable[, 4]
Fulltable$Unmapped_pc <- Fulltable[, 9] / Fulltable[, 4]
Fulltable$Trim_pc <- Fulltable[, 4] / Fulltable[, 3]

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

### Note 18-07-2019: Redo merging from here ###
### I renamed tpm_batch1_and2 to batch1
# Make a list of df's to edit all at once later, start with gene read transcripts per million
tpm.list <-
  list(tpm_batch1.txt,
       tpm_batch2.txt,
       tpm_batch3.txt,
       tpm_batch4.txt)
tpm_df <- do.call(rbind, tpm.list)
tpm_df <- as.data.frame(unlist(str_split_fixed(tpm_df[, 1], "\t", 2)))
colnames(tpm_df) <- c("Name", "Genes1tpm")
tpm_df$Genes1tpm <- as.numeric(as.character(tpm_df$Genes1tpm))
tpm_df$Name <- as.character(tpm_df$Name)

# Gene_counts are library complexity 
gene_counts_list <-
  list(
    gene_counts_batch1.txt,
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

# Gene_1read are the genes with at least one read per cell
gene_1read_list <-
  list(
    gene_1read_batch1.txt,
    gene_1read_batch2.txt,
    gene_1read_batch3.txt,
    gene_1read_batch4.txt
  )
gene_1read_df <- do.call(rbind, gene_1read_list)
gene_1read_df <- 
  as.data.frame(unlist(str_split_fixed(gene_1read_df[, 1], "\t", 3)[, 2:3]))
colnames(gene_1read_df) <- c("Name", "gene_1read")
gene_1read_df$gene_1read <-
  as.numeric(as.character(gene_1read_df$gene_1read))

### Note 19-07-2019 add readsum ###
setwd("~/Oxford 2.0/HNF1A/Tables/sumcheck3/")
readsum_files <- list.files(path = ".", pattern = ".txt")
for (i in 1:length(readsum_files))
  assign(readsum_files[i],
         read.csv(readsum_files[i], header =
                    FALSE))
readsum.list <-
  list(Batch1.txt,
       Batch2.txt,
       Batch3.txt,
       Batch4.txt)
readsum_df <- do.call(rbind, readsum.list)
readsum_df <- as.data.frame(unlist(str_split_fixed(readsum_df[, 1], "\t", 2)))
colnames(readsum_df) <- c("Name", "ReadSum")
readsum_df$ReadSum <- as.numeric(as.character(readsum_df$ReadSum))
readsum_df$Name <- as.character(readsum_df$Name)

### End of adding readsum ###

Fullerbatch <- merge(
  Fulltable,
  tpm_df,
  by.x = "Name",
  by.y = "Name",
  all.x = TRUE
)

Intertable <- merge(
  Fullerbatch,
  gene_1read_df,
  by.x = "Name",
  by.y = "Name",
  all.x = TRUE
)

Intertable2 <- merge(
  Intertable,
  gene_counts_df,
  by.x = "Name",
  by.y = "Name",
  all.x = TRUE
)

Fullertable <- merge(
  Intertable2,
  readsum_df,
  by.x = "Name",
  by.y = "Name",
  all.x = TRUE
)

#Test: readsum/assigned reads ratio:
Fullertable$assign_ratio <- Fullertable[, 22]/Fullertable[, 5]

p <- ggplot(Fullertable,
            aes(BatchName, assign_ratio, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "ColSum/Assigned Reads ratio") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)  



# Make columns for better plotting:
Fullertable$cellLine <- Fullertable$Name
Fullertable$cellLine <- strsplit(Fullertable$cellLine, '.X')
Fullertable$Stage <- lapply(Fullertable$cellLine, `[[`, 2)
Fullertable$Stage <- unlist(Fullertable$Stage)

Fullertable$Stage <- as.character(Fullertable$Stage)
Fullertable$Stage <- strsplit(Fullertable$Stage, "[.]")
Fullertable$Differentiation <- lapply(Fullertable$Stage, `[[`, 1)
Fullertable$cellLine <- lapply(Fullertable$cellLine, `[[`, 1)

#Rename accordingly, the dash in some SFC012N19218 lines, 
#made it seem as if there was another unique cell line.
SFC012_N19218_row <- Fullertable$cellLine == "SFC012N19218"
Fullertable$cellLine[SFC012_N19218_row] <- "SFC012.N19218"

#Keep the actual 'Stage' part from the Stage column, also keep the 'sampleType'
Fullertable$Stage <- lapply(Fullertable$Stage, `[[`, 2)
Fullertable$Stage <- as.character(Fullertable$Stage)
Fullertable$Stage <- strsplit(Fullertable$Stage, '_')
Fullertable$sampleType <- lapply(Fullertable$Stage, tail, n=1)
Fullertable$Stage <- lapply(Fullertable$Stage, `[[`, 1)

Fullertable$Sample  <- paste(Fullertable$cellLine, Fullertable$Differentiation, sep = "_X")
Fullertable$Sample  <- paste(Fullertable$Sample, Fullertable$Stage, sep = "_")
Fullertable$Sample <- paste(Fullertable$Sample, Fullertable$sampleType, sep = "_")
Fullertable$sampleType <- unlist(Fullertable$sampleType)
Fullertable$Differentiation <- unlist(Fullertable$Differentiation)
Fullertable$Stage <- unlist(Fullertable$Stage)
Fullertable$cellLine <- unlist(Fullertable$cellLine)

# Add correction column:
Pro291fsinsC_row <- Fullertable$cellLine %in% c("SFC012.0420CLN","SFC012.0420CL4")
Fullertable$Correction[Pro291fsinsC_row] <- "Pro291fsinsC/+" 
Corrected_row <- Fullertable$cellLine %in% c("SFC012.C14222", "SFC012.N19218")
Fullertable$Correction[Corrected_row] <- "Corrected/+"
SB_row	<- Fullertable$cellLine %in% "SBAD3.1"
Fullertable$Correction[SB_row] <- "+/+(WT)"
# Stop here and save a few large tables:
# 1) Fullertable:
setwd("~/Oxford 2.0/HNF1A/Tables/")
write.table(Fullertable, file = "Fullertable.txt", row.names = FALSE)

#Select single cells and controls
Fullertable_ss_select <-
  Fullertable[grep("singlecell", Fullertable$sampleType), ]
Fullertable_control_select <-
  Fullertable[grep("singlecell", Fullertable$sampleType, invert = TRUE), ]
Fullertable_control_select$TypeBatch <- paste(Fullertable_control_select$sampleType,
                                              Fullertable_control_select$BatchName, sep = ' ')

### Boxplots per Batch ###

# Boxplots for assigned:
p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Assigned, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Assigned_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Assigned %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)  

p <- ggplot(Fullertable_control_select,
            aes(BatchType, Assigned, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9
),
axis.title.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, Assigned_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank()) + labs(x = "Batch", y = "Assigned %") +  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

# Boxplots for unassigned ambiguity:
p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Unassigned_Ambiguity , fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Ambiguity_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Unassigned Ambiguity %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)  

p <- ggplot(Fullertable_control_select,
            aes(BatchType, Unassigned_Ambiguity, fill =  BatchType))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9
),
axis.title.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, Ambiguity_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Unassigned Ambiguity %") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9
  ),
  axis.title.x = element_blank()
)

# Boxplots for No Features:
p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Unassigned_NoFeatures , fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, NoFeatures_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "No Features %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, Unassigned_NoFeatures, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "No Features %") + theme(
  axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9
),
axis.title.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, NoFeatures_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "No Features %") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9
),
axis.title.x = element_blank())

# Boxplots for Unmapped:
p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Unassigned_Unmapped , fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Unmapped_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Unmapped_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Unmapped %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(BatchType, Unassigned_Unmapped, fill =  BatchType))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9
),
axis.title.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, Unmapped_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Unmapped %") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())

# gene count data: 1TPM, gene_1read and %top100
p <- ggplot(Fullertable_ss_select,
            aes(BatchName, Genes1tpm, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Genes 1TPM") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, Genes1tpm, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Genes1TPM") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9),
  axis.title.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, gene_1read, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Genes 1 read") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, gene_1read, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Genes 1 read") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9),
  axis.title.x = element_blank())


p <- ggplot(Fullertable_ss_select,
            aes(BatchName, top100_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Top 100 reads %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(BatchType, top100_reads_pc, fill =  BatchType))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(TypeBatch, top100_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Top 100 reads %") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9),
  axis.title.x = element_blank())
### End of per batch boxplots ###

### By correction ### 

# Trimmed reads
p <- ggplot(Fullertable_ss_select,
            aes(Correction, TrimReads, fill =  Correction))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable_ss_select,
            aes(BatchName, TrimReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(x = "Batch", y = "Top 100 reads %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(Correction, TrimReads, fill =  Correction))
p + geom_boxplot() + theme_minimal()

# Genes > 1TPM
p <- ggplot(Fullertable_ss_select,
            aes(Correction, Genes1tpm, fill =  Correction))
p + geom_boxplot() + theme_minimal() + labs(x = "Correction", y = "Genes 1TPM") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(Correction, Genes1tpm, fill =  Correction))
p + geom_boxplot() + theme_minimal() + labs(x = "Correction", y = "Genes1TPM") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9),
  axis.title.x = element_blank())

# Top 100
p <- ggplot(Fullertable_ss_select,
            aes(Correction, top100_reads_pc, fill =  Correction))
p + geom_boxplot() + theme_minimal() + labs(x = "Correction", y = "Top 100 reads %") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20), legend.position = "none"
)

p <- ggplot(Fullertable_control_select,
            aes(Correction, top100_reads_pc, fill =  Correction))
p + geom_boxplot() + theme_minimal() + labs(x = "Correction", y = "Top 100 reads %") + theme(
  axis.title.y = element_text(size = 18),
  legend.position = "none",
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9),
  axis.title.x = element_blank())

### End of by correction ###

### Per sample ###
# Make some changes for boxplot aesthetics:
# Clip off '_single cell' from name. 

Fullertable_ss_select$Sample <- gsub("_", " ", Fullertable_ss_select$Sample)
Fullertable_ss_select$Sample <- gsub(" singlecell", "", Fullertable_ss_select$Sample)
test <- str_split_fixed(Fullertable_ss_select$Sample, " ", 3)[, c(2,3,1)]
cols <- cbind(test[, 1] , test[, 2], test[, 3])
Fullertable_ss_select$Sample <- apply(cols , 1 , paste , collapse = " " )

Fullertable_ss_select$Sample <- gsub("SFC012.0420CLN", "MODY1", Fullertable_ss_select$Sample)
Fullertable_ss_select$Sample <- gsub("SFC012.0420CL4", "MODY2", Fullertable_ss_select$Sample)
Fullertable_ss_select$Sample <- gsub("SFC012.C14222", "GE1", Fullertable_ss_select$Sample)
Fullertable_ss_select$Sample <- gsub("SFC012.N19218", "GE2", Fullertable_ss_select$Sample)

plot(Fullertable$Genes1tpm, Fullertable$top100_reads_pc)
plot(Fullertable$Genes1tpm, Fullertable$Assigned)
plot(Fullertable$Genes1tpm, Fullertable$gene_1read)

p <- ggplot(Fullertable_ss_select,
            aes(Sample, TrimReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(y = "Trimmed Reads") + theme(
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(
  angle = 90,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Sample, Genes1tpm, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(y = "Genes > 1 TPM") + theme(
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(
  angle = 90,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Sample, top100_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(y = "Top 100 Reads %") + theme(
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(
  angle = 90,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Sample, Assigned_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal() + labs(y = "Top 100 Reads %") + theme(
  axis.title.y = element_text(size = 18),
  axis.text.x = element_text(
  angle = 90,
  hjust = 1,
  vjust = 0.9),
  axis.title.x = element_blank())




#Boxplots for single cell stages and single cell diffs
p <- ggplot(Fullertable_ss_select,
            aes(Stage, RawReads, fill =  Stage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Stage, top100_reads_pc, fill =  Stage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Stage, Genes1tpm, fill =  Stage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Differentiation, RawReads, fill =  Differentiation))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Differentiation, top100_reads_pc, fill =  Differentiation))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_ss_select,
            aes(Differentiation, Genes1tpm, fill =  Differentiation))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

#Boxplots for control stages, select control first
Fullertable_control_select$typeStage <- paste(Fullertable_control_select$sampleType, 
                                              Fullertable_control_select$Stage, sep = '_')
p <- ggplot(Fullertable_control_select,
            aes(typeStage, RawReads, fill =  typeStage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(typeStage, top100_reads_pc, fill =  typeStage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(typeStage, Genes1tpm, fill =  typeStage))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

# Do the same you did for stages with differentiations
Fullertable_control_select$typeDiff <- paste(Fullertable_control_select$sampleType, 
                                              Fullertable_control_select$Differentiation, sep = '_')

p <- ggplot(Fullertable_control_select,
            aes(typeDiff, RawReads, fill =  typeDiff))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(typeDiff, top100_reads_pc, fill =  typeDiff))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

p <- ggplot(Fullertable_control_select,
            aes(typeDiff, Genes1tpm, fill =  typeDiff))
p + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_blank())

#Boxplots for sample types
p <- ggplot(Fullertable,
            aes(sampleType, Genes1tpm, fill =  sampleType))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable,
            aes(sampleType, Assigned_pc, fill =  sampleType))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable,
            aes(sampleType, top100_reads_pc, fill =  sampleType))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fullertable,
            aes(sampleType, RawReads, fill =  sampleType))
p + geom_boxplot() + theme_minimal()


#Some quick plots [Save for later reference, but move eventually]
p <- ggplot(Fulltable, aes(BatchName, RawReads, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, Assigned, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable,
         aes(BatchName, Unassigned_MultiMapping, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable,
         aes(BatchName, Unassigned_NoFeatures, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable,
         aes(BatchName, Unassigned_Unmapped, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, RawSum, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, Assigned_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, MultiMapping_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <-
  ggplot(Fulltable, aes(BatchName, NoFeatures_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, Unmapped_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()

p <- ggplot(Fulltable, aes(BatchName, Raw_reads_pc, fill =  BatchName))
p + geom_boxplot() + theme_minimal()


