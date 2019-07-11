#!/bin/bash
origins="/well/mccarthy/production/rna-seq/data/StemBANCC/HNF1A_SS2_scRNA_remap"
#jobs
./start_large_table.sh "${origins}/Batch1/180410_K00150_0313_BHTV35BBXX/trim_fastq" \
/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch1_trim

./start_large_table.sh "${origins}/Batch2/180524_K00198_0314_BHV55CBBXX/trim_fastq" \
/well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch2_trim

./start_large_table.sh "${origins}/Batch2/181205_K00198_0374_BH2NJKBBXY/trim_fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch2second_trim

./start_large_table.sh "${origins}/Batch3/190322_K00181_0125_AH25MHBBXY/trim_fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch3_trim

./start_large_table.sh "${origins}/Batch4/190328_K00181_0127_BH25MJBBXY/trim_fastq" \
 /well/mccarthy/users/maxlouis/oxford2/HNF1A_project/data/Batch4_trim
