# Import the 4 tpm tables
# There might be a more elegant way of doing this, but just run this script for all separately,
# Then move the results to local pc for Rstudio merging in 2.large_table.R

#Input is a "xxx.gene.tpm.tsv" file
#Output is a table of cell sample names and the number of transcripts per million
#Automate this to do it for multiple files using the .sh script on /well/mccarthy/users/maxlouis

calc_1tpm <- function(inFile, outFile) {
  tpm_counts=read.table(inFile,h=T,row.names=1)
  genes_1tpm=colSums(tpm_counts[,2:ncol(tpm_counts)]>=1)
  write.table(genes_1tpm, file = outFile, sep = "\t", col.names = F)
}

#parse command line argument and run this function

args = commandArgs(trailingOnly = TRUE)
if (length(args)!=2) {
  stop("Need input and output file name, e-val = NA ===> e-val = 0\n",call. = TRUE)
} else {
  inFile=args[1]
  outFile=args[2]
}

calc_1tpm(inFile,outFile)
