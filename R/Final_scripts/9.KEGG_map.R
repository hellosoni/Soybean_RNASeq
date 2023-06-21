# Load packages
library(tidyverse)
library(KEGGREST) #BiocManager::install("KEGGREST")
library(rentrez)

# Source script with functions
source("/fs/ess/PAS0471/nghi/rnaseq2/R/enrich_funs.R")
#source("scripts/KEGG_fun.R")

# Define input files
geneID_map_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/geneID_map.txt"  # Created by 'scripts/geneID_map.R'

# Define output files
outdir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ko_map_file <- file.path(outdir, "ko_map.txt")
pw_map_file <- file.path(outdir, "pw_map.txt")
pw_df_file <- file.path(outdir, "pw_info.txt")


# GET NCBI GENE_IDS ------------------------------------------------------------
# It is necessary to get NCBI gene IDs because those will be reported by KEGG
# ...when getting the genes associated with each pathway
geneID_map <- read_tsv(geneID_map_file, show_col_types = FALSE) |>
  drop_na() |>
  group_by(geneID) |>
  filter(n() == 1) |> # Remove IDs with multiple NCBI matches
  ungroup()


# GET KEGG PATHWAYS -------------------------------------------------------------
# Get KEGG pathways ("gmx" = Glycine max)
pw_list <- keggList("pathway", "gmx")
pw_ids <- sub("path:", "", names(pw_list))
pw_df <- tibble(category = pw_ids, description = unname(pw_list))
write_tsv(pw_df, pw_df_file)

# Get all genes (NCBI IDs) associated with each pathway
pw_genes_raw <- sapply(pw_ids, get_pw_genes)
pw_genes <- tibble(KEGG_pathway = names(pw_genes_raw),
                   geneID_NCBI = pw_genes_raw) %>%
  unnest(cols = geneID_NCBI)


# CREATE FINAL KEGG MAP --------------------------------------------------------
pw_map <- merge(pw_genes, geneID_map, by = "geneID_NCBI") %>%
  select(category = KEGG_pathway,
         gene_id = geneID) %>%
  arrange(category)

# Report
message("Nr of unique KEGG pathways:  ", length(unique(pw_map$category)))
message("Nr of gene2pathway links:    ", dim(pw_map))

# Save to file
write_tsv(pw_map, pw_map_file)
