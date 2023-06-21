# SETUP ------------------------------------------------------------------------
# Load Packages
library(tidyverse)
library(clusterProfiler) # BiocManager::install("clusterProfiler")

# Source script with functions
source("/fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/rnaseq/rfuns/enrich_funs.R")
# Define input files
kegg_map_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/pw_map.txt"
kegg_info_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/pw_info.txt"

FILE_PATTERN <- "all-res.*txt"
DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")
# Output files
outdir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
kegg_res_file <- file.path(outdir, "kegg_res.txt")
kegg_sig_file <- file.path(outdir, "kegg_sig.txt")

# Read the gene length df - should have columns "gene_id" and "length"
gene_lens <- read_tsv(gene_len_file, show_col_types = FALSE)

# Read the DE results
DE_res <- read_tsv(DE_files, id = "time", show_col_types = FALSE) %>%
  mutate(time = as.integer(gsub(".*DE_T(\\d\\d).tsv", "\\1", time)))

DE_res <- DE_res |>
  mutate(group1 = sub("M92_220", "M92", group1),
         group2 = sub("M92_220", "M92", group2)) |>
  separate_wider_delim(cols = group1, names = c("line1", "treat1"), delim = "_") |>
  separate_wider_delim(cols = group2, names = c("line2", "treat2"), delim = "_") |>
  mutate(gene_id = sub(".Wm82.a4.v1", "", gene),
         treat1 = sub("Inoculated", "Inoc", treat1),
         treat2 = sub("Inoculated", "Inoc", treat2),
         contrast = paste0(line1, treat1, "_", line2, treat2, "_", time)) %>%
  select(-gene)
# Read input files
kegg_info <- read_tsv(kegg_info_file, show_col_types = FALSE)
kegg_map <- read_tsv(kegg_map_file, show_col_types = FALSE) |>
  left_join(kegg_info, by = "category")

# RUN THE KEGG ENRICHMENT ANALYSIS ---------------------------------------------
# Run for all contrasts
# Run the kegg analysis for all contrasts
kegg_all <- map_dfr(.x = unique(DE_res$contrast), .f = run_enrich,
                    DE_res = DE_res, cat_map = kegg_map, return_df = TRUE)

# Only keep significant terms
kegg_sig <- kegg_all %>% filter(sig == 1)
kegg_sig %>% count(contrast)

# Write to file
write_tsv(kegg_all, kegg_res_file)
write_tsv(kegg_sig, kegg_sig_file)


