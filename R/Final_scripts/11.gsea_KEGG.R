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
kegg_map_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/pw_map.txt"
kegg_info_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/pw_info.txt"
kegg_res_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/kegg/kegg_res.txt"
FILE_PATTERN <- "*all-res.*txt"

DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")

# Define output files
outdir <- "rnaseq2/results/gsea"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "XX.txt")


# PREP INPUT DATA --------------------------------------------------------------
# Read the kegg results - for comparison, and to get the kegg term descriptions
kegg_res <- read_tsv(kegg_res_file, show_col_types = FALSE)


# Make the kegg map
kegg_names <- kegg_res |> select(category, description) |> distinct(.keep_all = TRUE)

# Make the kegg map
kegg_names <- kegg_res |> select(category, description) |> distinct(.keep_all = TRUE)
kegg_map <- read.delim(kegg_map_file, sep = "\t") |>
  left_join(kegg_names, by = "category")


# Read the DE results
col_names <- c("gene_id", read_tsv(DE_files[1]) %>% colnames())

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
#TODO - Could also try with shrunken LFC?

# Run all and return a regular dataframe
gsea_res <- map_dfr(.x = unique(DE_res$contrast), .f = run_gsea,
                    DE_res = DE_res, cat_map = kegg_map,p_enrich = 0.05, return_df = TRUE)


# PLOT GSEA RESULTS ------------------------------------------------------------
# Run GSEA for one contrast, and keep gsea format for plots
# Filter results to be less than pvalue cutoff of 0.05
#Conrad
gsea_res <- run_gsea("ConradInoc_ConradMock_24", DE_res, kegg_map,p_enrich = 0.05, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Conrad24_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
#Plot sig terms
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 24hpi")
#


gsea_res <- run_gsea("ConradInoc_ConradMock_48", DE_res, kegg_map,p_enrich = 0.05, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Conrad48_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
#Plot sig terms
nsig <- sum(gsea_res$p.adjust < 0.05)
gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 48hpi")


#

gsea_res <- run_gsea("ConradInoc_ConradMock_72", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Conrad72_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)

# Plots - https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
nsig <- sum(gsea_res$p.adjust < 0.05)

gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2)+ ggtitle("Conrad at 72hpi")

ridgeplot(gsea_res) +
  theme(axis.text.y = element_text(size = 10))


#M92-220
gsea_res <- run_gsea("M92Inoc_M92Mock_48", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "M92_48_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)

gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("M92_220 at 48hpi")

#M92-220
gsea_res <- run_gsea("M92Inoc_M92Mock_72", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "M92_72_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)

# Plots - https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
nsig <- sum(gsea_res$p.adjust < 0.05)

gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("M92_220 at 72hpi")

ridgeplot(gsea_res) +
  theme(axis.text.y = element_text(size = 10))


#Sloan
gsea_res <- run_gsea("SloanInoc_SloanMock_24", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Sloan24_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
gsea_res2 <- pairwise_termsim(gsea_res)
nsig <- sum(gsea_res$p.adjust < 0.05)
treeplot(gsea_res2) + ggtitle("Sloan at 24hpi")



gsea_res <- run_gsea("SloanInoc_SloanMock_48", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Sloan48_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
gsea_res2 <- pairwise_termsim(gsea_res)
nsig <- sum(gsea_res$p.adjust < 0.05)
treeplot(gsea_res2) + ggtitle("Sloan at 48hpi")

gsea_res <- run_gsea("SloanInoc_SloanMock_72", DE_res, p_enrich = 0.05, kegg_map, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "Sloan72_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)

# Plots - https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
nsig <- sum(gsea_res$p.adjust < 0.05)
dotplot(gsea_res, showCategory = nsig)
heatplot(gsea_res, showCategory = nsig) #foldChange = geneList,

gsea_res2 <- pairwise_termsim(gsea_res)
treeplot(gsea_res2) + ggtitle("Sloan at 72hpi")

upsetplot(gsea_res)

ridgeplot(gsea_res) +
  theme(axis.text.y = element_text(size = 10))


#Get a certain pathways with names
#Note: Run all gsea contrasts first

#focal_pathways <- gsea_res$Description[grep("pathogen", gsea_res$Description)]
#Compare M92 & Conrad at 48hpi-----

gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "ConradInoc_ConradMock_48")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 60)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 48hpi")


# Compare M92_220 & Conrad at 72hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_72", "ConradInoc_ConradMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 60)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 72hpi")

#Compare M92_220 and Sloan at 48hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "SloanInoc_SloanMock_48")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 48hpi")
#Compare M92_220 and Sloan at 72hpi
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_72", "SloanInoc_SloanMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 72hpi")

#enrich_plot(enrich_res = gsea_res_plot, x_var = contrast)

#Compare Conrad and Sloan
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c( "ConradInoc_ConradMock_72","SloanInoc_SloanMock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 72hpi")



#Plot 1 line at 3 time points
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("SloanInoc_SloanMock_24", "SloanInoc_SloanMock_48", "SloanInoc_SloanMock_72")) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)
p + facet_grid(.~type)

#

gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("ConradInoc_ConradMock_48", "ConradInoc_ConradMock_72")) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type)

#
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "M92Inoc_M92Mock_72")) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type)

#FN0170228
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "FN0170228Inoc_FN0170228Mock_48")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 48hpi")

##
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_72", "FN0170228Inoc_FN0170228Mock_72")) %>%
  #filter(Description %in% focal_pathways) %>%
  mutate(Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 30)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3)

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"


p <-ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(direction = -1)

p + facet_grid(.~type) +   ggtitle("KEGG pathway enrichment at 72hpi")


# Comparing M92 & FN0140856----
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "M92Inoc_M92Mock_72","FN0140856Inoc_FN0140856Mock_48", "FN0140856Inoc_FN0140856Mock_72")) %>%
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


p + facet_grid(.~type) + ggtitle("KEGG pathway enrichment for FN0140856 vs M92-220")

# Comparing M92 & FN0170228---
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "M92Inoc_M92Mock_72","FN0170228Inoc_FN0170228Mock_48", "FN0170228Inoc_FN0170228Mock_72")) %>%
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


p + facet_grid(.~type) + ggtitle("KEGG pathway enrichment for FN0170228 vs M92-220")

# Comparing M92 & FNMN0329----
gsea_res_plot <- gsea_res %>%
  filter(contrast %in% c("M92Inoc_M92Mock_48", "M92Inoc_M92Mock_72","FNMN0329Inoc_FNMN0329Mock_48", "FNMN0329Inoc_FNMN0329Mock_72")) %>%
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


p + facet_grid(.~type) + ggtitle("KEGG pathway enrichment for FNMN0329 vs M92-220")




