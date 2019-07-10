#!/bin/bash
origins="/well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_SS2_scRNA_remap"
#jobs
./2.start_large_table.sh "${origins}/Batch1/180410_K00150_0313_BHTV35BBXX/fastq" \
/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch1

./2.start_large_table.sh "${origins}/Batch2/180524_K00198_0314_BHV55CBBXX/fastq" \
/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch2

./2.start_large_table.sh "${origins}/Batch2/181205_K00198_0374_BH2NJKBBXY/fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch2second

./2.start_large_table.sh "${origins}/Batch3/190322_K00181_0125_AH25MHBBXY/fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch3

./2.start_large_table.sh "${origins}/Batch4/190328_K00181_0127_BH25MJBBXY/fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch4
