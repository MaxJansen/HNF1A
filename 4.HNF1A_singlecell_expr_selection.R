### Run locally in R3.6.3
### This script filters out bad single cells ###
# You start with 1824 cells

# Load libraries
library(dplyr)

# Go to directory
setwd("~/Oxford 2.0/HNF1A/Tables/")

# Get the tables
Fullertable <- read.table("Fullertable.txt", header = T)

# Some checks
dim(Fullertable)
colnames(Fullertable)
rownames(Fullertable)

# Remove non-single cell
SC_table <- Fullertable[grep("singlecell", Fullertable$Name),]


# Start selecting, remove bottom 5% of genes1TPM
TPMselect <- SC_table[SC_table$Genes1tpm > quantile(SC_table$Genes1tpm, 0.03), ]
top100select <- TPMselect[TPMselect$top100_reads_pc < 0.80, ]
unmappedselect <- top100select[top100select$Unmapped_pc < 0.90, ]
assignedpcselect <- unmappedselect[unmappedselect$Assigned_pc > 0.1, ]
Rawreadselect <- assignedpcselect[assignedpcselect$TrimReads > quantile(assignedpcselect$TrimReads, 0.02), ]

#What percentage have we removed? We still have 93.5%
nrow(Rawreadselect)/nrow(SC_table) *100

Rawreadselect$Name
unique(Rawreadselect$cellLine)

Selected_list <- Rawreadselect$Name
write.table(Rawreadselect, file = "sc_filtered_metadata.txt")
write.csv(Selected_list, file = "filtered_sc_names.txt")
