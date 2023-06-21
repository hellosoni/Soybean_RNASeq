# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(ggridges)

# Source script with functions
source("/fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/rnaseq/rfuns/enrich_funs.R")

# Define input files
GO_map_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO/GO_map.txt"
gene_len_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO/gene_lengths.txt"
GO_res_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/GO/GO_res.txt"

FILE_PATTERN <- "all-res.*txt"

DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")

# Define output files
outdir <- "rnaseq2/results/gsea"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "XX.txt")


# PREP INPUT DATA --------------------------------------------------------------
# Read the GO results - for comparison, and to get the GO term descriptions
GO_res <- read_tsv(GO_res_file, show_col_types = FALSE)

# Make the GO map
GO_names <- GO_res |> select(category, description) |> distinct(.keep_all = TRUE)
GO_map <- read.delim(GO_map_file, sep = "\t") |>
  select(GO_term, gene_id) |>
  dplyr::rename(category = GO_term) |>
  left_join(GO_names, by = "category")

# Read the DE results
DE_res <- read_tsv(DE_files, id = "time", show_col_types = FALSE) %>%
  mutate(time = as.integer(gsub(".*DE_T(\\d\\d).tsv", "\\1", time))) %>%
  rename(log2FoldChange = lfc)

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


# RUN GSEA ---------------------------------------------------------------------

# Run all and return a regular dataframe
gsea_res <- map_dfr(.x = unique(DE_res$contrast), .f = run_gsea,
                    DE_res = DE_res, cat_map = GO_map,p_enrich = 0.05, return_df = TRUE)


# PLOT GSEA RESULTS ------------------------------------------------------------
# Run GSEA for one contrast, and keep gsea format for plots
#Conrad
gsea_res <- run_gsea("ConradInoc_ConradMock_24", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Conrad24_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 24hpi")


gsea_res <- run_gsea("ConradInoc_ConradMock_48", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Conrad48_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 48hpi")

gsea_res <- run_gsea("ConradInoc_ConradMock_72", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Conrad72_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots -
nsig <- sum(gsea_res$p.adjust < 0.05)

gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 72hpi")

#M92-220
gsea_res <- run_gsea("M92Inoc_M92Mock_48", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "M92_48_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("M92_220 at 48hpi")

gsea_res <- run_gsea("M92Inoc_M92Mock_72", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "M9272_GO_sig_72.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots -
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("M92_220 at 72hpi")


#Sloan
gsea_res <- run_gsea("SloanInoc_SloanMock_24", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Sloan_24_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)

treeplot(gsea_res2) + ggtitle("Sloan at 24hpi")



gsea_res <- run_gsea("SloanInoc_SloanMock_48", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Sloan_48_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)

treeplot(gsea_res2) + ggtitle("Sloan at 48hpi")

gsea_res <- run_gsea("SloanInoc_SloanMock_72", DE_res,p_enrich = 0.05, GO_map, return_df = FALSE)
gseaGO_results <- gsea_res@result
#Write result to file
GO_sig_file <-file.path(outdir, "Sloan72_GO_sig.txt")
write_tsv(gseaGO_results,GO_sig_file)

# Plots -
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("Sloan at 72hpi")


##PLot 2 lines comparison
#Note: Need to run all the contrast first
# Comparing M92 & Conrad at 48hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "ConradInoc_ConradMock_48")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 60)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)


p + facet_grid(.~type) + ggtitle("GO pathway enrichment at 48hpi")

# Comparing M92 & Conrad at 72hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_72", "ConradInoc_ConradMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 60)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)


p + facet_grid(.~type) + ggtitle("GO pathway enrichment at 72hpi")


### Comparing M92 & Sloan at 48hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "SloanInoc_SloanMock_48")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)


p + facet_grid(.~type) + ggtitle("GO pathway enrichment at 48hpi")

#Compare M92 & Sloan at 72hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_72", "SloanInoc_SloanMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)


p + facet_grid(.~type) + ggtitle("GO pathway enrichment at 72hpi")
#Compare Sloan at 3 timepoints
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("SloanInoc_SloanMock_24", "SloanInoc_SloanMock_48", "SloanInoc_SloanMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 40)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

#Compare Conrad at 3 timepoints
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("ConradInoc_ConradMock_24", "ConradInoc_ConradMock_48", "ConradInoc_ConradMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 40)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

#Compare M92_220 at 2 timepoints
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c( "M92Inoc_M92Mock_48", "M92Inoc_M92Mock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 40)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

