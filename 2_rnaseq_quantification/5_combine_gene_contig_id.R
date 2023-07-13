library(tidyverse)

setwd("/Volumes/green_groups_nme_public/JoaoSousa/tmp/")

# Load contig annotation
devNiklas_gene_contig_id      <- read_tsv(file = "devNiklas_gene_contig_id.txt", col_names = F)
Jochen_annotation_contig_gene <- read_tsv(file = "Jochen_annotation_contig_gene.txt", col_names = F)

# Combine Jochen and devNiklas contig IDs
joined_em_en_mi_dev_tpm_2_filtered <- rbind(Jochen_annotation_contig_gene, devNiklas_gene_contig_id) %>% unique

# Save the combined contig IDs
write_tsv(joined_em_en_mi_dev_tpm_2_filtered, file = "joined_em_en_mi_dev_tpm_2_filtered_genes.txt", col_names = F)