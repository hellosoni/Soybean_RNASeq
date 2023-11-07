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


FILE_PATTERN <- "res_mocks.*txt"
DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE/Mocks"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = FILE_PATTERN)
DE_files <- DE_files[-c(31, 32)] # Problems with reading these files
# Source script with functions

# Define output files
outdir <- "rnaseq2/results/gsea"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "XX.txt")

# PREP INPUT DATA --------------------------------------------------------------
# PREP INPUT DATA --------------------------------------------------------------
# Read the kegg results - for comparison, and to get the kegg term descriptions
kegg_res <- read_tsv(kegg_res_file, show_col_types = FALSE)


# Make the kegg map
kegg_names <- kegg_res |> select(category, description) |> distinct(.keep_all = TRUE)

# Make the kegg map
kegg_names <- kegg_res |> select(category, description) |> distinct(.keep_all = TRUE)
kegg_map <- read.delim(kegg_map_file, sep = "\t") |>
  left_join(kegg_names, by = "category")

col_names <- c("gene_id", read_tsv(DE_files[1]) %>% colnames())
DE_res <- read_tsv(DE_files,
                   skip = 1,
                   col_names = col_names,
                   id = "file",
                   show_col_types = FALSE)

DE_res <- DE_res %>% select (-file) |>
  mutate(level1 = sub("M92_220", "M92", level1),
         level2 = sub("M92_220", "M92", level2)) |>
  separate_wider_delim(cols = level1, names = c("line1", "treat1"), delim = "_") |>
  separate_wider_delim(cols = level2, names = c("line2", "treat2"), delim = "_") |>
  mutate(gene_id = sub(".Wm82.a4.v1", "", gene_id),
         treat1 = sub("Mocks", "Mock", treat1),
         treat2 = sub("Mocks", "Mock", treat2),
         contrast = paste0(line1, treat1, "_", line2, treat2))


#M92Mock_SloanMock
gsea_res <- run_gsea("M92Mock_SloanMock", DE_res, kegg_map,p_enrich = 0.05, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "M92SloanMock_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
#Plot sig terms

nsig <- sum(gsea_res$p.adjust < 0.05)
dotplot(gsea_res, showCategory = nsig)

ridgeplot(gsea_res) +
  theme(axis.text.y = element_text(size = 10))
#Conrad Mock vs Sloan Mock
gsea_res <- run_gsea("ConradMock_SloanMock", DE_res, kegg_map,p_enrich = 0.05, return_df = FALSE)
gseaKEGG_results <- gsea_res@result
#Write result to file
kegg_sig_file <-file.path(outdir, "ConradSloanMock_kegg_sig.txt")
write_tsv(gseaKEGG_results,kegg_sig_file)
#Plot sig terms

nsig <- sum(gsea_res$p.adjust < 0.05)
dotplot(gsea_res, showCategory = nsig)

ridgeplot(gsea_res) +
  theme(axis.text.y = element_text(size = 10))


#Comparing two contrasts
#Run gsea for all contrasts
gsea_res <- map_dfr(.x = unique(DE_res$contrast), .f = run_gsea,
                    DE_res = DE_res, cat_map = kegg_map,p_enrich = 0.05, return_df = TRUE)
# Defense pathways
def_pw <- c("Terpenoid backbone biosynthesis", "Phenylpropanoid biosynthesis", "Plant-pathogen interaction", "MAPK signaling pathway", "Flavonoid biosynthesis", "Isoflavonoid biosynthesis","Isoquinoline alkaloid biosynthesis","Tropane, piperidine and pyridine alkaloid biosynthesis","Zeatin biosynthesis", "Thiamine metabolism","Isoquinoline alkaloid biosynthesis", "Glucosinolate biosynthesis","Sesquiterpenoid and triterpenoid biosynthesis","Glutathione metabolism" )

gsea_res_plot <- gsea_res %>%
  filter(!is.na(Description)) %>%
  filter(contrast %in% c("ConradMock_SloanMock", "M92Mock_SloanMock")) %>%
  mutate(contrast = ifelse(grepl("Conrad", contrast), "Conrad", "M92-220"),
         Description = sub(" - .*", "", Description),
         Description = str_trunc(Description, width = 50)) %>%
  filter(p.adjust < 0.05) %>%
  add_count(category) %>%
  filter(n < 3) %>%
  mutate(shared = ifelse(n == 2, "shared", contrast)) %>%
  mutate(Description = ifelse(Description %in% def_pw,
                              paste0("*** ", Description),
                              Description))

gsea_res_plot$type = "Upregulated"
gsea_res_plot$type[gsea_res_plot$NES < 0] = "Downregulated"

ggplot(data = gsea_res_plot) +
  aes(x = contrast, y = Description, color = p.adjust, size = NES) +
  geom_point() +
  facet_grid(shared~type, scales = "free_y", space = "free_y") +
  scale_color_viridis_c(direction = -1) +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey95"))




