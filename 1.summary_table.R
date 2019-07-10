### This is a script to make a summary table of the RG_info files ###
##
library(tidyr)
library(dplyr)

setwd("~/Oxford 2.0/HNF1A/Tables")

#Read tables and add Batch names accordingly
RG_info1 <- read.table("RG_info1")
RG_info2 <- read.table("RG_info2a")
RG_info3 <- read.table("RG_info3")
RG_info4 <- read.table("RG_info4")

RG_info1$Batch <- as.factor("Batch1")
RG_info2$Batch <- as.factor("Batch2")
RG_info3$Batch <- as.factor("Batch3")
RG_info4$Batch <- as.factor("Batch4")

#Make a list of df's to edit all at once later
my.list <- list(RG_info1,RG_info2,RG_info3,RG_info4)

#slim down the tables, keep name, date and batch. Bind df's and reformat column names.
select_column <- function(y) {
  z <- subset(y, select=c(V3,V7,Batch))  
  return(z)} 

my.list <- lapply(my.list, select_column)
RG_infoAll <- do.call(rbind, my.list)
RG_infoAll$V3 <- gsub('SM=', '', RG_infoAll$V3)
RG_infoAll$V3 <- gsub('"', '', RG_infoAll$V3)
RG_infoAll$V7 <- gsub('DT=', '', RG_infoAll$V7)
RG_infoAll$V7 <- gsub('"', '', RG_infoAll$V7)
colnames(RG_infoAll)[colnames(RG_infoAll)=="V3"] <- "sampleName" 
colnames(RG_infoAll)[colnames(RG_infoAll)=="V7"] <- "Date" 
RG_infoAll$Date <- as.factor(RG_infoAll$Date)

#Copy first column for later steps
RG_infoAll$cellLine <- RG_infoAll$sampleName
RG_infoAll$cellLine <- strsplit(RG_infoAll$cellLine, '-X')
RG_infoAll$Stage <- lapply(RG_infoAll$cellLine, `[[`, 2)
RG_infoAll$Stage <- unlist(RG_infoAll$Stage)

RG_infoAll$Stage <- strsplit(RG_infoAll$Stage, '-')
RG_infoAll$Differentiation <- lapply(RG_infoAll$Stage, `[[`, 1)

RG_infoAll$cellLine <- lapply(RG_infoAll$cellLine, `[[`, 1)

#Check how many unique cell lines there are
unique(RG_infoAll$cellLine)

#Rename accordingly, the dash in some SFC012N19218 lines, 
#made it seem as if there was another unique cell line.
SFC012_N19218_row <- RG_infoAll$cellLine == "SFC012N19218" 
RG_infoAll$cellLine[SFC012_N19218_row] <- "SFC012-N19218"
  
#Make a MODY column, based on the cell lines.
Pro291fsinsC_row <- RG_infoAll$cellLine %in% c("SFC012-0420CLN","SFC012-0420CL4")
RG_infoAll$MODY[Pro291fsinsC_row] <- "Pro291fsinsC/+" 
Corrected_row <- RG_infoAll$cellLine %in% c("SFC012-C14222", "SFC012-N19218")
RG_infoAll$MODY[Corrected_row] <- "Corrected/+"
SB_row	<- RG_infoAll$cellLine %in% "SBAD3-1"
RG_infoAll$MODY[SB_row] <- "+/+(WT)"

#Keep the actual 'Stage' part from the Stage column, also keep the 'sampleType'
RG_infoAll$Stage <- lapply(RG_infoAll$Stage, `[[`, 2)
RG_infoAll$Stage <- as.character(RG_infoAll$Stage)
RG_infoAll$Stage <- strsplit(RG_infoAll$Stage, '_')
RG_infoAll$sampleType <- lapply(RG_infoAll$Stage, tail, n=1)
RG_infoAll$Stage <- lapply(RG_infoAll$Stage, `[[`, 1)

#Convert to factors, useful for final summary
RG_infoAll$cellLine <- as.factor(as.character(RG_infoAll$cellLine))
RG_infoAll$Stage <- as.factor(as.character(RG_infoAll$Stage))
RG_infoAll$Differentiation <- as.factor(as.numeric(RG_infoAll$Differentiation))
RG_infoAll$MODY <- as.factor(RG_infoAll$MODY)
RG_infoAll$sampleType <- as.factor(as.character(RG_infoAll$sampleType))

#Summarise
summary_table <- summary(RG_infoAll)
show_me <- table(RG_infoAll)

#Save RG_infoAll table and summary
write.table(RG_infoAll, file = "complete_sample_table.txt", row.names = FALSE)
write.table(Summary_table, file = "summary_table.txt", row.names = FALSE)



### Some stuff in directory, using bash ####

#Find Batches for summaries by Agata. Bash scripts by Agata:
#cd /well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_SS2_scRNA_remap
#cat */*/RG_info | awk '{print $3}' | sort -u | wc -l
#cat */*/RG_info | awk '{print $3}' | sort -u | perl -i -pe 's/SM=//g' | perl -i -pe 's/"//g' | awk -F "_" '{print $1}' | sort | uniq -c

#Make dataframe from output of last script:
counts <- c(96, 192, 96, 96, 192, 192, 96, 192, 96, 96, 192, 192, 96, 96)
samples <- c('SBAD3-1-X2-PE','SBAD3-1-X2-ST7', 'SFC012-0420CL4-X1-ST7', 'SFC012-0420CL4-X2-PE', 
             'SFC012-0420CL4-X2-ST7', 'SFC012-0420CLN-X1-ST7', 'SFC012-0420CLN-X2-PE', 'SFC012-0420CLN-X2-ST7', 
             'SFC012-C14222-X1-ST7', 'SFC012-C14222-X2-PE', 'SFC012-C14222-X2-ST7', 'SFC012-N19218-X2-ST7', 
             'SFC012N19218-X1-ST7', 'SFC012N19218-X2-PE')

overview <- data.frame(samples, counts)

#Use this df to see to which Batches 
Batch_collection <- c()
for (i in overview$samples) {
   Batch_collection <- c(Batch_collection, c(i, unique(RG_infoAll$Batch[grep(i, RG_infoAll$sampleName)])))
 }

#Batch occurrences
#"SBAD3-1-X2-PE"         "2"                     
#"SBAD3-1-X2-ST7"        "3"                     
#"SFC012-0420CL4-X1-ST7" "1"                     
#"SFC012-0420CL4-X2-PE"  "2"                     
#"SFC012-0420CL4-X2-ST7" "3"                    "4"                     
#"SFC012-0420CLN-X1-ST7" "1"                     "2"                     
#"SFC012-0420CLN-X2-PE" "2"                     
#"SFC012-0420CLN-X2-ST7" "3"                     "4"                     
#"SFC012-C14222-X1-ST7" "1"                     
#"SFC012-C14222-X2-PE"   "2"                     
#"SFC012-C14222-X2-ST7"  "3"                    "4"                     
#"SFC012-N19218-X2-ST7"  "3"                     "4"                     
#"SFC012N19218-X1-ST7"  "1"                     
#"SFC012N19218-X2-PE"    "2"


### End of Bash stuff in directory ###
RG_infoAll$essentials <- paste(RG_infoAll$Batch, RG_infoAll$cellLine, RG_infoAll$Stage, 
                               RG_infoAll$Differentiation, RG_infoAll$MODY, RG_infoAll$sampleType)
desired_count_table <- aggregate(RG_infoAll$sampleName, by=list(essentials=RG_infoAll$essentials), FUN = length)


watch_this <- desired_count_table %>% separate(essentials, c("Batch", "cellLine", "Stage", "Differentiation", "Correction", "sampleType"), sep = " ")
write.table(watch_this, file = "desired_count_table.txt", row.names = FALSE)

class(desired_count_table$essentials)
