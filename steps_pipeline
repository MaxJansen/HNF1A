#This is all based on the 'README.txt' in /well/mccarthy/production/rna-seq/code/scRNA_pipeline
#See that script for details and explanations.
#Step 1. Load dependencies and run trimming
module load java/1.8.0_latest
module load samtools
module load python/3.5.2-gcc5.4.0/

perl  /well/mccarthy/production/rna-seq/code/scRNA_pipeline/1.trim_map_dups.pl -d StemBANCC/HNF1A_scRNA_Max/Batch1/180410_K00150_0313_BHTV35BBXX -q c
perl  /well/mccarthy/production/rna-seq/code/scRNA_pipeline/1.trim_map_dups.pl -d StemBANCC/HNF1A_scRNA_Max/Batch2/180524_K00198_0314_BHV55CBBXX -q c
perl  /well/mccarthy/production/rna-seq/code/scRNA_pipeline/1.trim_map_dups.pl -d StemBANCC/HNF1A_scRNA_Max/Batch2/181205_K00198_0374_BH2NJKBBXY -q c
perl  /well/mccarthy/production/rna-seq/code/scRNA_pipeline/1.trim_map_dups.pl -d StemBANCC/HNF1A_scRNA_Max/Batch3/190322_K00181_0125_AH25MHBBXY -q c
perl  /well/mccarthy/production/rna-seq/code/scRNA_pipeline/1.trim_map_dups.pl -d StemBANCC/HNF1A_scRNA_Max/Batch4/190328_K00181_0127_BH25MJBBXY -q c

#Step 2. Merge project. Add "--lib-type unstranded", don't do "keep"
perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/2.merge_project_and_featureCount_quantify.pl -q c -d StemBANCC/HNF1A_scRNA_Max/Batch1 --lib-type unstranded
perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/2.merge_project_and_featureCount_quantify.pl -q c -d StemBANCC/HNF1A_scRNA_Max/Batch2 --lib-type unstranded
perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/2.merge_project_and_featureCount_quantify.pl -q c -d StemBANCC/HNF1A_scRNA_Max/Batch3 --lib-type unstranded
perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/2.merge_project_and_featureCount_quantify.pl -q c -d StemBANCC/HNF1A_scRNA_Max/Batch4 --lib-type unstranded

#Step 3.
module load R/3.2.2
/apps/well/perl/5.16.3/bin/perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/5.combine_project_quants.pl --data StemBANCC/HNF1A_scRNA_Max/Batch1 --genome GRCh37
/apps/well/perl/5.16.3/bin/perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/5.combine_project_quants.pl --data StemBANCC/HNF1A_scRNA_Max/Batch2 --genome GRCh37
/apps/well/perl/5.16.3/bin/perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/5.combine_project_quants.pl --data StemBANCC/HNF1A_scRNA_Max/Batch3 --genome GRCh37
/apps/well/perl/5.16.3/bin/perl /well/mccarthy/production/rna-seq/code/scRNA_pipeline/5.combine_project_quants.pl --data StemBANCC/HNF1A_scRNA_Max/Batch4 --genome GRCh37

#Redo! Step 3.

> count_f_name[!(count_f_name %in% summary_names)]
  [1] "gene0"      "gene4"      "gene6"      "gene15"     "gene14"
  [6] "gene17"     "gene20"     "gene22"     "gene23"     "gene24"
 [11] "gene34"     "rna4"       "rna5"       "gene67"     "rna6"
 [16] "gene71"     "gene70"     "gene73"     "rna9"       "gene126"
 [21] "gene132"    "rna12"      "gene144"    "gene145"    "gene151"
 [26] "ERCC-00002" "ERCC-00003" "ERCC-00004" "ERCC-00009" "ERCC-00012"
 [31] "ERCC-00013" "ERCC-00014" "ERCC-00016" "ERCC-00017" "ERCC-00019"
 [36] "ERCC-00022" "ERCC-00024" "ERCC-00025" "ERCC-00028" "ERCC-00031"
 [41] "ERCC-00033" "ERCC-00034" "ERCC-00035" "ERCC-00039" "ERCC-00040"
 [46] "ERCC-00041" "ERCC-00042" "ERCC-00043" "ERCC-00044" "ERCC-00046"
 [51] "ERCC-00048" "ERCC-00051" "ERCC-00053" "ERCC-00054" "ERCC-00057"
 [56] "ERCC-00058" "ERCC-00059" "ERCC-00060" "ERCC-00061" "ERCC-00062"
 [61] "ERCC-00067" "ERCC-00069" "ERCC-00071" "ERCC-00073" "ERCC-00074"
 [66] "ERCC-00075" "ERCC-00076" "ERCC-00077" "ERCC-00078" "ERCC-00079"
 [71] "ERCC-00081" "ERCC-00083" "ERCC-00084" "ERCC-00085" "ERCC-00086"
 [76] "ERCC-00092" "ERCC-00095" "ERCC-00096" "ERCC-00097" "ERCC-00098"
 [81] "ERCC-00099" "ERCC-00104" "ERCC-00108" "ERCC-00109" "ERCC-00111"
 [86] "ERCC-00112" "ERCC-00113" "ERCC-00116" "ERCC-00117" "ERCC-00120"
 [91] "ERCC-00123" "ERCC-00126" "ERCC-00130" "ERCC-00131" "ERCC-00134"
 [96] "ERCC-00136" "ERCC-00137" "ERCC-00138" "ERCC-00142" "ERCC-00143"
[101] "ERCC-00144" "ERCC-00145" "ERCC-00147" "ERCC-00148" "ERCC-00150"
[106] "ERCC-00154" "ERCC-00156" "ERCC-00157" "ERCC-00158" "ERCC-00160"
[111] "ERCC-00162" "ERCC-00163" "ERCC-00164" "ERCC-00165" "ERCC-00168"
[116] "ERCC-00170" "ERCC-00171"
