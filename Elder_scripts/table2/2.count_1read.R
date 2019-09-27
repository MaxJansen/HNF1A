# Import the 4 gene count tables (look into what they're actually called)
# There might be a more elegant way of doing this, but just run this script for all separately,
# Then move the results to local pc for Rstudio merging in 2.large_table.R

#Don't just list numbers, also add column names for merging later
library("data.table")

calc_1read <- function(inFile, outFile) {
  gene_counts = data.frame(fread(inFile,header = TRUE), row.names = TRUE)
  genes_1read = as.numeric(colSums(gene_counts[, 2:ncol(gene_counts)] >= 1))
  final_count = as.data.frame(cbind(colnames(gene_counts)[2:ncol(gene_counts)], genes_1read))
  final_count$genes_1read = as.numeric(as.character(final_count$genes_1read))
  write.table(final_count, file = outFile, sep = "\t", col.names = F)
}

#parse command line argument and run this function

args = commandArgs(trailingOnly = TRUE)
if (length(args)!=2) {
  stop("Need input and output file name, e-val = NA ===> e-val = 0\n",call. = TRUE)
} else {
  inFile=args[1]
  outFile=args[2]
}

calc_1read(inFile,outFile)
