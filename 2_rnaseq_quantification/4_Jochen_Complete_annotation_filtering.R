library(tidyverse)

setwd("/Volumes/green_groups_nme_public/JoaoSousa/tmp/")

# Loading Jochen contig annotation
Jochen_annotation <- 
  read_tsv("Jochen_Complete_annotation_hsa_bta_mmu_Unique_on_data_212.txt", col_names = F) %>% 
  unique %>% 
  dplyr::rename(contig_id = X1, bovine_gene_symbol = X9, human_gene_symbol = X16, mouse_gene_symbol = X23) %>%
  select(contig_id, bovine_gene_symbol, human_gene_symbol, mouse_gene_symbol) %>% 
  mutate(mouse_gene_symbol = str_to_upper(mouse_gene_symbol)) 

# Filter genes after BLASTN
blastn <- read_tsv(file = "devNiklas_JochenAssembly_blastn.out", col_names = F, skip = 2)

genes_to_exclude_blastn <- 
Jochen_annotation %>% 
  mutate(gene_symbol = case_when(bovine_gene_symbol != "-" ~ bovine_gene_symbol, 
                                 bovine_gene_symbol == "-" & human_gene_symbol != "-" ~ human_gene_symbol,
                                 bovine_gene_symbol == "-" & human_gene_symbol == "-" & mouse_gene_symbol != "-" ~ mouse_gene_symbol,
                                 TRUE ~ "")) %>% 
  mutate(species = case_when(bovine_gene_symbol != "-" ~ "bovine", 
                             bovine_gene_symbol == "-" & human_gene_symbol != "-" ~ "human",
                             bovine_gene_symbol == "-" & human_gene_symbol == "-" & mouse_gene_symbol != "-" ~ "mouse",
                             TRUE ~ "")) %>% 
  select(gene_symbol, contig_id, species) %>% 
  unique %>% 
  filter(contig_id %in% unique(pull(blastn, X2))) %>% 
  pull(gene_symbol) %>% 
  unique

# Filter out genes with the same name as in the developmental gene set
trinity_filtered_table <- read_tsv(file = "devNiklas_gene_contig_id.txt", col_names = F)
dev_gene_symbol <- trinity_filtered_table %>% pull(X1) %>% unique

# Combining developmental gene set and BLASTN genes to exclude
genes_to_exclude <- unique(c(dev_gene_symbol, genes_to_exclude_blastn))

# Save filtered table
Jochen_annotation %>% 
  mutate(gene_symbol = case_when(bovine_gene_symbol != "-" ~ bovine_gene_symbol, 
                                 bovine_gene_symbol == "-" & human_gene_symbol != "-" ~ human_gene_symbol,
                                 bovine_gene_symbol == "-" & human_gene_symbol == "-" & mouse_gene_symbol != "-" ~ mouse_gene_symbol,
                                 TRUE ~ "")) %>% 
  mutate(species = case_when(bovine_gene_symbol != "-" ~ "bovine", 
                             bovine_gene_symbol == "-" & human_gene_symbol != "-" ~ "human",
                             bovine_gene_symbol == "-" & human_gene_symbol == "-" & mouse_gene_symbol != "-" ~ "mouse",
                             TRUE ~ "")) %>% 
  select(gene_symbol, contig_id, species) %>% 
  filter(!gene_symbol %in% genes_to_exclude) %>% 
  unique %>% 
  write_tsv("Jochen_annotation_contig_gene.txt", col_names = F)
