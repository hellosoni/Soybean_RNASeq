# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(goseq)
library(enrichplot)
library(ggupset)
library(ggridges)
library(clusterProfiler)


# Source script with functions
source("/fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/rnaseq/rfuns/enrich_funs.R")

# Input files
GO_map_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO/GO_map.txt"
gene_len_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO/gene_lengths.txt"

DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")

# Output files
outdir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
GO_res_file <- file.path(outdir, "GO_res.txt")
GO_sig_file <- file.path(outdir, "GO_sig.txt")

# Read the GO map - should have columns "gene_id" and "go_term", and needs to be a regular df (not a tibble)
GO_map <- read.delim(GO_map_file, sep = "\t")

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

# USE THE WRAPPER FUNCTION -----------------------------------------------------
# Run the GO analysis for all contrasts
GO_res <- map_dfr(.x = unique(DE_res$contrast), .f = run_GO,
                  DE_res = DE_res, GO_map = GO_map, gene_lens = gene_lens)

# Only keep significant terms
GO_sig <- GO_res %>% filter(sig == 1)
GO_sig %>% count(contrast)

# Write to file
write_tsv(GO_res, GO_res_file)
write_tsv(GO_sig, GO_sig_file)



