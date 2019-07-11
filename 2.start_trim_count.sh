#!/bin/bash
#Go to directory, based on arguments
DIRECTORY1=$1
cd $DIRECTORY1
#Determine output
OUTPUT1=$2
#zcat multiple files and count the lines
for file in ./*_1.fastq.gz; do
  echo $file >> "${OUTPUT1}_1.txt"
  zcat $file | wc -l >> "${OUTPUT1}_1.txt"
  file_2=${file/_1/_2}
  echo $file_2 >> "${OUTPUT1}_2.txt"
  zcat $file_2 | wc -l >> "${OUTPUT1}_2.txt"
done;
