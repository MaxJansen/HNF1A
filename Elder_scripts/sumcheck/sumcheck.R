# Quick script to get colsums (sum of all gene reads) per cell.
# Compare these sums to assigned reads values from .summary-files

setwd("/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/sumcheck")
sourcedir <- "/well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_SS2_scRNA_remap/"
batchSummer <- function(batchfile, outfile) {
	batch_df <- read.table(paste0(sourcedir, batchfile), header = TRUE, row.names = 1, sep = "\t")
	sumlist <- colSums(batch_df[, 2:ncol(batch_df)])
	write.table(sumlist, file = outfile, sep = "\t", col.names = F)
	return(sumlist)
}

batchSummer("Batch1/17.07.2019.Batch1.gene.counts.tsv", "Batch1.txt")
batchSummer("Batch2/19.12.2018.Batch2.gene.counts.tsv", "Batch2.txt")
batchSummer("Batch3/17.04.2019.Batch3.gene.counts.tsv", "Batch3.txt")
batchSummer("Batch4/17.04.2019.Batch4.gene.counts.tsv", "Batch4.txt")

