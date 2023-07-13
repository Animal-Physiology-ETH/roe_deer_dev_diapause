#!/bin/bash

module load gcc/4.8.2 samtools/1.8

date

for i in ./HISAT2/*.bam
do
  echo "processing $i"
  samtools fastq -@ 4 -F 4 "$i" >> hisat2_mapped_reads.fastq
done

echo "done extracting fastq - compressing file"
gzip hisat2_mapped_reads.fastq
