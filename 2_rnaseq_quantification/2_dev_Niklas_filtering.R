library("tidyverse")

setwd("/Volumes/green_groups_nme_public/JoaoSousa/tmp/")

# Loading the annotation table of developmental genes
Galaxy898_bta_hsa_mmu <- 
read_tsv(file = "/Volumes/green_groups_nme_public/JoaoSousa/tmp/devNiklas_bta_hsa_mmu.tsv", col_names = F)

# Renaming columns
Galaxy898_bta_hsa_mmu <-
  Galaxy898_bta_hsa_mmu %>% 
  dplyr::rename(contig_id = X1, gene_species_tx = X2, pident = X3, evalue = X11, bitscore = X12, qseq = X21, sseq = X22) %>%
  mutate(gene_id = str_remove(contig_id, "_i\\d+$")) %>% 
  mutate(gene_symbol = str_extract(gene_species_tx, "^[^_]*")) %>% 
  mutate(species = str_extract(gene_species_tx, "_.*_") %>% str_remove_all("_"))

# Selecting unique trinity_gene_id and gene symbol
trinity_gene_id_selected <- 
  Galaxy898_bta_hsa_mmu %>% 
  select(gene_id, gene_symbol, species, evalue, pident, bitscore) %>% 
  unique %>% 
  arrange(gene_id, evalue, dplyr::desc(bitscore), dplyr::desc(pident)) %>% 
  group_by(gene_id) %>% 
  filter(evalue < 0.05) %>% 
  filter(pident == max(pident)) %>%
  filter(bitscore == max(bitscore)) %>%
  unique %>% 
  ungroup %>% 
  group_by(gene_symbol) %>% 
  filter(evalue < 0.05) %>% 
  filter(pident == max(pident)) %>%
  filter(bitscore == max(bitscore)) %>%
  unique %>% 
  ungroup %>% 
  arrange(gene_id) %>%
  select(gene_id, gene_symbol, species) %>% 
  dplyr::rename(gene_symbol_selected = gene_symbol) %>% 
  unique %>% 
  group_by(gene_id) %>%
  mutate(priority = case_when(species == "bovine" ~ 1, species == "human" ~ 1, species == "mouse" ~ 1, TRUE ~ 4)) %>% 
  filter(species == min(species)) %>% 
  ungroup %>% 
  select(-priority)

# Checking if the trinity id and gene symbol are unique
stopifnot((trinity_gene_id_selected$gene_id %>% duplicated() %>% sum) == 0)
stopifnot((trinity_gene_id_selected$gene_symbol_selected %>% duplicated() %>% sum) == 0)

# Selecting unique contig ids and creating the final table
trinity_filtered_table <-
  Galaxy898_bta_hsa_mmu %>% 
  select(gene_id, contig_id, gene_symbol, species, evalue, pident, bitscore, qseq) %>% 
  left_join(trinity_gene_id_selected) %>% 
  unique %>% 
  arrange(contig_id, evalue, dplyr::desc(bitscore), dplyr::desc(pident)) %>% 
  group_by(contig_id) %>% 
  filter(gene_symbol == gene_symbol_selected) %>% 
  filter(evalue < 0.05) %>%
  filter(pident == max(pident)) %>%
  filter(bitscore == max(bitscore)) %>% 
  unique %>% 
  ungroup

# Checking if the contig ids are unique
stopifnot((trinity_filtered_table$contig_id %>% duplicated() %>% sum) == 0)
stopifnot((trinity_filtered_table %>% select(contig_id, gene_symbol) %>% unique %>% duplicated() %>% sum) == 0)

# Saving query as fasta file
trinity_filtered_table %>% 
  select(contig_id, qseq) %>% 
  unique %>% 
  mutate(qseq = str_replace_all(qseq, "-", "")) %>% 
  mutate(qseq = str_replace_all(qseq, " ", "")) %>% 
  dplyr::rename(seq.name = contig_id, seq.text = qseq) %>% 
  phylotools::dat2fasta(outfile = "devNiklas_bta_hsa_mmu_filtered.fasta")

# Saving gene symbol and contig id
write_tsv(trinity_filtered_table %>% select(gene_symbol, contig_id, species), file = "devNiklas_gene_contig_id.txt", col_names = F)