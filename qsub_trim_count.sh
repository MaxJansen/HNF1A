#!/bin/bash
qsub -V -cwd -N count_fastq -q short.qc -e count_fastq.err -o count_fastq.out batches_count_fastq.sh
