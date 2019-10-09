#!/usr/bin/Rscript

### This script splits single cell tables in PE and ST7 stage tables

stagesplit <- function(inFile) {
  library(tools)
  # Get input and output names right
  stopifnot(is.character(inFile))
  fileName <- file_path_sans_ext(inFile)
  PEName <- paste0(fileName, "_", "PE.txt")
  ST7Name <- paste0(fileName, "_", "ST7.txt")
  
  
  # Reading the input file
  mainTable <- read.table(inFile, header = TRUE)
  geneNames <- mainTable[1]
  
  # Perform the split and add geneName column:
  PETable <- mainTable[grep("PE", names(mainTable))]
  PETable <- cbind(geneNames, PETable)
  ST7Table <- mainTable[grep("ST7", names(mainTable))]
  ST7Table <- cbind(geneNames, ST7Table)
  
  
  # Write the two tables
  write.table(PETable,file = PEName)
  write.table(ST7Table, file = ST7Name)
  
  cat("... done splitting single cell table into two tables: one for ST7 and PE stage")
}
args = commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Need input file name, inFile",call. = TRUE)
} else {
  inFile=args[1]

}

stagesplit(inFile)