summary_file=list.files(path=".", pattern=".gene.count.qc.summary.tsv$")
summary=read.table(summary_file,h=T,row.names=1)

gene_counts_file=list.files(path=".", pattern=".gene.counts.tsv$")
gene_counts=read.table(gene_counts_file,h=T,row.names=1)

tpm_counts_file=list.files(path=".", pattern=".gene.tpm.tsv$")
tpm_counts=read.table(tpm_counts_file,h=T,row.names=1)

base=gsub(".gene.counts.tsv","",gene_counts_file)

### get % of reads that map to top 100 genes
reads_pc=rep(0,ncol(gene_counts)-1)
for (i in 2:ncol(gene_counts)){
	tmp=gene_counts[,i]
	names(tmp) = rownames(gene_counts)
	top_100 = names(sort(tmp,decreasing=T))[1:100]
	reads_pc[i-1] = sum(gene_counts[top_100,i])/sum(tmp)
	
}

## pc reads in MT genes:
library("annotables")
ens_genes=sapply(strsplit(rownames(gene_counts),split="\\."), function(x) x[1])
MT_genes = which(ens_genes %in% grch37$ensgene[which(grch37$chr=="MT")])
MT_reads_pc=apply(gene_counts[,2:ncol(gene_counts)],2, function(x) sum(x[MT_genes])/sum(x))

## pc reads in ribosomal genes: biotype=rRNA
rRNA_genes = which(ens_genes %in% grch37$ensgene[which(grch37$biotype=="rRNA")])
rRNA_reads_pc=apply(gene_counts[,2:ncol(gene_counts)],2, function(x) sum(x[rRNA_genes])/sum(x))

## pc reads in protein_coding genes: biotype=protein_coding
prot_coding_genes = which(ens_genes %in% grch37$ensgene[which(grch37$biotype=="protein_coding")])
prot_coding_reads_pc=apply(gene_counts[,2:ncol(gene_counts)],2, function(x) sum(x[prot_coding_genes])/sum(x))

stats=data.frame(name=colnames(summary), category=sapply(strsplit(colnames(summary), split="_"), function(x) x[3]), plate=sapply(strsplit(colnames(summary), split="_"), function(x) x[1]), well=sapply(strsplit(colnames(summary), split="_"), function(x) x[2]), raw_reads=colSums(summary), reads_assigned=as.numeric(summary[11,]), reads_assigned_pc= as.numeric(summary[11,])/colSums(summary), reads_unmapped=as.numeric(summary[6,]),  reads_unmapped_pc= as.numeric(summary[6,])/colSums(summary), reads_nofeatures=as.numeric(summary[7,]), reads_nofeatures_pc=as.numeric(summary[7,])/colSums(summary), reads_multimap=as.numeric(summary[10,]), reads_multimap_pc= as.numeric(summary[10,])/colSums(summary) , top100g_reads_pc = reads_pc, genes_1read=as.numeric(colSums(gene_counts[,2:ncol(gene_counts)]>=1)), genes_1tpm=as.numeric(colSums(tpm_counts[,2:ncol(tpm_counts)]>=1)), MT_reads_pc=MT_reads_pc, rRNA_reads_pc=rRNA_reads_pc, prot_coding_reads_pc=prot_coding_reads_pc, INS_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="INS",2:ncol(tpm_counts)]), GCG_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="GCG",2:ncol(tpm_counts)]), NKX6_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="NKX6-1",2:ncol(tpm_counts)]), PDX1_cells = as.numeric(tpm_counts[tpm_counts$GeneName=="PDX1",2:ncol(tpm_counts)]))

### for SB_Ad3 - test:
stats=data.frame(name=colnames(summary), sample=sapply(strsplit(colnames(summary), split="_"), function(x) x[2]), 
category=sapply(strsplit(colnames(summary), split="_"), function(x) x[3]), plate=sapply(strsplit(colnames(summary), split="_"), function(x) x[4]), 
raw_reads=colSums(summary), reads_assigned=as.numeric(summary[11,]), reads_assigned_pc= as.numeric(summary[11,])/colSums(summary), reads_unmapped=as.numeric(summary[6,]),  reads_unmapped_pc= as.numeric(summary[6,])/colSums(summary), reads_nofeatures=as.numeric(summary[7,]), reads_nofeatures_pc=as.numeric(summary[7,])/colSums(summary), reads_multimap=as.numeric(summary[10,]), reads_multimap_pc= as.numeric(summary[10,])/colSums(summary) , top100g_reads_pc = reads_pc, genes_1read=as.numeric(colSums(gene_counts[,2:ncol(gene_counts)]>=1)), genes_1tpm=as.numeric(colSums(tpm_counts[,2:ncol(tpm_counts)]>=1)), MT_reads_pc=MT_reads_pc, rRNA_reads_pc=rRNA_reads_pc, prot_coding_reads_pc=prot_coding_reads_pc, INS_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="INS",2:ncol(tpm_counts)]), GCG_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="GCG",2:ncol(tpm_counts)]), NKX6_cells=as.numeric(tpm_counts[tpm_counts$GeneName=="NKX6-1",2:ncol(tpm_counts)]), PDX1_cells = as.numeric(tpm_counts[tpm_counts$GeneName=="PDX1",2:ncol(tpm_counts)]))

write.table(stats, paste0(base,".map_qc_stats.txt"),sep="\t",quote=F)


#### make some boxplots of the data:

# use ggplot2 colors:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = length(levels(stats$category))
cols = gg_color_hue(n)


pdf(paste0(base,".qc_plots.single_cell.pdf"), width=10)
for(i in 5:19){
	boxplot(stats[,i] ~ stats$category, col=cols, main=names(stats)[i], cex.lab=2, cex.axis=2)
}
dev.off()


##########  the same - but split by sample for the single-cells ############
stats$group=factor(paste(stats$category, stats$plate, sep="_"))
n = length(levels(stats$group))

# n = length(levels(stats$plate))
cols = gg_color_hue(n)

pdf(paste0(base,".qc_plots.by_plate.single_cell.pdf"), width=10)
par(las=2, mar=c(13,6,4,1))
for(i in 5:19){
#	boxplot(stats[which(stats$category=="singlecell"),i] ~ stats$plate[which(stats$category=="singlecell")], col=cols, main=names(stats)[i], cex.lab=2, cex.axis=1)
	boxplot(stats[which(stats$sample=="sc"),i] ~ stats$group[which(stats$sample=="sc")], col=cols, main=names(stats)[i], cex.lab=2, cex.axis=1)

}
dev.off()

###
 rg1=read.table("180410_K00150_0313_BHTV35BBXX/RG_info")
 rg1$sm=gsub("-",".",gsub("SM=","",gsub("\"","",rg1$V3)))
 rg2=read.table("180524_K00198_0314_BHV55CBBXX/RG_info")
 rg2$sm=gsub("-",".",gsub("SM=","",gsub("\"","",rg2$V3)))
 
 stats_0420CLN = stats[which(stats$plate=="SFC012.0420CLN.X1.ST7"),]
 stats_0420CLN$batch=rep(1,nrow(stats_0420CLN))
 stats_0420CLN[which(stats_0420CLN$name %in% rg2$sm),]$batch = 2
 
 
pdf(paste0(base,".qc_plots.by_seq_batch.single_cell.pdf"))
for(i in 5:19){
	boxplot(stats_0420CLN[,i] ~ stats_0420CLN$batch, col=cols, main=names(stats)[i], cex.lab=2, cex.axis=1)
}
dev.off()


######### boxplots - based on RNASeqMetrics 

qc_metrics=list.files(path="merged/QC/", pattern=".RNA_Metrics",recursive=T)

names=scan(paste0("merged/QC/",qc_metrics[1]), '', skip = 6, nlines = 1, sep = '\t')
values <- read.table(paste0("merged/QC/",qc_metrics[1]), skip = 7, nrows = 1,  sep = '\t')
values[,28]=gsub(".RNA_Metrics","",qc_metrics[1])
for (i in 2:length(qc_metrics)){
	tmp <- read.table(paste0("merged/QC/",qc_metrics[i]), skip = 7, nrows = 1, sep = '\t')
	values=rbind(values, tmp)
	values[i,28]=gsub(".RNA_Metrics","",qc_metrics[i])
}
names(values)=names
values$category=sapply(strsplit(colnames(summary), split="_"), function(x) x[2])
values$condition=sapply(strsplit(colnames(summary), split="_"), function(x) x[3])
values$plate=sapply(strsplit(colnames(summary), split="_"), function(x) x[4])
values$group=paste(values$condition,values$plate,sep="_")

pdf(paste0(base,".qc_plots.RNA_Metrics.pdf"))
par(las=2, mar=c(8,4,4,1))
for(i in 17:22){
	boxplot(values[which(values$category=="sc"),i] ~ values[which(values$category=="sc"),]$group, col=cols, main=names(values)[i], cex.lab=2, cex.axis=1)
}
for(i in 5:19){
	boxplot(stats[which(stats$sample=="sc"),i] ~ stats$group[which(stats$sample=="sc")], col=cols, main=names(stats)[i], cex.lab=2, cex.axis=1)
}
par(mfrow=c(2,4))
for (g in levels(stats$group)){
	plot(stats[stats$group==g,]$top100g_reads_pc, stats[stats$group==g,]$genes_1tpm, col=as.numeric(stats[stats$group==g,]$sample), pch=16, main=g, xlab="pc_reads_100top_genes",ylab="genes at 1TPM")
}
for (g in levels(stats$group)){
	plot(stats[stats$group==g,]$top100g_reads_pc, stats[stats$group==g,]$genes_1tpm, col=as.numeric(stats[stats$group==g,]$sample), pch=16, main=g, 
	xlab="pc_reads_100top_genes",ylab="genes at 1TPM", ylim=c(2000,6000))
}

dev.off()

### gene body coverage plots:
pdf(paste0(base,".gene_body_coverage_plots.RNA_Metrics.pdf"), width=18, height=12)
par(mfrow=c(4,6))
groups=unique(values$group)
for (g in groups){
	for (i in 1:length(qc_metrics)){
		sample=gsub(".RNA_Metrics","",qc_metrics[i])
		if(grepl(g,gsub(".RNA_Metrics","",qc_metrics[i]))){
			tmp <- read.table(paste0("merged/QC/",qc_metrics[i]), skip = 10, h=T, sep = '\t')
			plot(tmp[,1],tmp[,2],type="l", main=sample)
		}
	}
}
dev.off()

__END__
## htsjdk.samtools.metrics.StringHeader
# CollectRnaSeqMetrics REF_FLAT=/well/htseq/Genomes/annotation/human_g1k_v37.ERCC.refFlat STRAND_SPECIFICITY=NONE INPUT=B03_sc_cond1_NXT1.bam OUTPUT=../QC/B03_sc_cond1_NXT1.RNA_Metrics    MINIMUM_LENGTH=500 RRNA_FRAGMENT_PERCENTAGE=0.8 METRIC_ACCUMULATION_LEVEL=[ALL_READS] ASSUME_SORTED=true STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
## htsjdk.samtools.metrics.StringHeader
# Started on: Fri Nov 09 14:09:14 GMT 2018

## METRICS CLASS        picard.analysis.RnaSeqMetrics
PF_BASES        PF_ALIGNED_BASES        RIBOSOMAL_BASES CODING_BASES    UTR_BASES       INTRONIC_BASES  INTERGENIC_BASES        IGNORED_READS   CORRECT_STRAND_READS    INCORRECT_STRAND_READS  NUM_R1_TRANSCRIPT_STRAND_READS  NUM_R2_TRANSCRIPT_STRAND_READS  NUM_UNEXPLAINED_READS   PCT_R1_TRANSCRIPT_STRAND_READS  PCT_R2_TRANSCRIPT_STRAND_READS  PCT_RIBOSOMAL_BASES     PCT_CODING_BASES        PCT_UTR_BASES   PCT_INTRONIC_BASES      PCT_INTERGENIC_BASES    PCT_MRNA_BASES  PCT_USABLE_BASES        PCT_CORRECT_STRAND_READS
        MEDIAN_CV_COVERAGE      MEDIAN_5PRIME_BIAS      MEDIAN_3PRIME_BIAS      MEDIAN_5PRIME_TO_3PRIME_BIAS    SAMPLE  LIBRARY READ_GROUP
187271846       138942331               51910432        64431567        16125536        6474796 0       0       0       319937  311090  76023   0.50701 0.49299         0.373611        0.463729        0.116059        0.046601        0.83734 0.621247        0       0.543281        0.464354        0.812246        0.433234

## HISTOGRAM    java.lang.Integer
normalized_position     All_Reads.normalized_coverage
0       0.335047
1       0.382586
2       0.427105
3       0.46859
4       0.510559
5       0.540059
6       0.57151
7       0.595478
8       0.617553
9       0.632466
10      0.660029
11      0.688365
12      0.704931
13      0.725365
14      0.749698
15      0.767501
16      0.791717
17      0.813551
18      0.831106
19      0.85886
20      0.878134
21      0.88754
22      0.906961
23      0.930191
24      0.944179
25      0.956291
26      0.975839
######### coverage histograms