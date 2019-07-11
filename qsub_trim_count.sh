#!/bin/bash
qsub -V -cwd -N count_trim -q short.qc -e count_trim.err -o count_trim.out batches_trim_fastq.sh
