# SETUP ------------------------------------------------------------------------
## Loading packages
library(tidyverse)
library(here)
library(ape)

## Input files
annot_info_file <- "rnaseq2/ref_genomes/A4/annotation/Gmax_508_Wm82.a4.v1.annotation_info.txt"
gff_file <- "rnaseq2/ref_genomes/A4/annotation/Gmax_508_Wm82.a4.v1.gene.gff3"

## Output files
GO_map_file <- "rnaseq2/results/GO/GO_map.txt"
gene_len_file <- "rnaseq2/results/GO/gene_lengths.txt"

## Read input files
annot_info <- read_tsv(annot_info_file)
gff <- read.gff(gff_file)

# Create GO map ----------------------------------------------------------------
## Count the max nr of GO terms for a gene
max_terms <- max(str_count(annot_info$GO, pattern = ",") , na.rm = TRUE)+ 1
## Create a vector with
go_colnames <- paste0("term", seq_len(max_terms))

go_map <- annot_info %>%
  select(gene_id = locusName, GO) %>%  # Get rid of other columns
  # Spread all GO terms out into separate columns:
  separate(col = "GO", sep = ",", into = go_colnames, fill = "right") %>%
  # Pivot from wide to long format with one row per gene-GO-term combination:
  pivot_longer(-gene_id, names_to = "term_idx", values_to = "GO_term") %>%
  drop_na() %>%                     # Remove rows with NAs
  select(-term_idx) %>%             # Remove the term-index column
  distinct()                        # Remove duplicate rows (due to transcripts)

# Create gene length dataframe -------------------------------------------------
gene_len <- gff %>%
  filter(type == "gene") %>%        # Only select genes (not exons etc)
  mutate(length = end - start + 1,  # Compute gene length
         # Extract gene ID from attributes column:
         gene_id = sub(".*Name=(Glyma[^;]*).*", "\\1", attributes)) %>%
  select(gene_id, length) %>%       # Only select these 2 columns
  arrange(gene_id)                  # Sort by gene ID


# Save output ------------------------------------------------------------------
write_tsv(go_map, GO_map_file)
write_tsv(gene_len, gene_len_file)
