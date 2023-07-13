#!/bin/bash

module load gcc/4.8.2 samtools/1.8

date

for i in ./HISAT2/*.bam
do
  echo "processing $i"
  samtools view "$i" | cut -f3 | uniq | cut  -d"_" -f1 | sort | uniq >> histat2_matched_genes.list 
done

sort histat2_matched_genes.list | uniq > histat2_matched_genes.list.uniq
