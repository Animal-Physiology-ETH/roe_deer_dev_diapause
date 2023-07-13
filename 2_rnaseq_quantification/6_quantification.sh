#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2048

parallel -j5 -k --plus 'rsem-calculate-expression -p 8 {} /cluster/work/nme/data/josousa/projects/roe_deer_embryo/genomes/rsem_genome_roe_deer/joined_em_en_mi_dev_tpm_2_filtered {/..}' ::: $(ls /cluster/work/nme/data/josousa/projects/roe_deer_embryo/data/raw/p2583_1/fastq/plate*/*.fastq.gz /cluster/work/nme/data/josousa/projects/roe_deer_embryo/data/raw/p2583_1/fastq/plate*/*.fastq.gz)