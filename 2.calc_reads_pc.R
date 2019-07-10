# Import the 4 tpm tables
# There might be a more elegant way of doing this, but just run this script for all separately,
# Then move the results to local pc for Rstudio merging in 2.large_table.R

#Input is a "xxx.gene.counts.tsv" file
#Output is a table of cell sample names and the percentage of reads mapping to top 100 genes
#Automate this to do it for multiple files using the .sh script on /well/mccarthy/users/maxlouis

calc_reads_pc <- function(inFile, outFile) {
  gene_counts=read.table(inFile,h=T,row.names=1)
  reads_pc=rep(0,ncol(gene_counts)-1)
  for (i in 2:ncol(gene_counts)){
    tmp=gene_counts[,i]
    names(tmp) = rownames(gene_counts)
    top_100 = names(sort(tmp,decreasing=T))[1:100]
    reads_pc[i-1] = sum(gene_counts[top_100,i])/sum(tmp)
    
  }
  final_count <- as.data.frame(cbind(colnames(gene_counts)[2:ncol(gene_counts)],reads_pc))
  final_count$reads_pc <- as.numeric(as.character(final_count$reads_pc))
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

calc_reads_pc(inFile,outFile)
