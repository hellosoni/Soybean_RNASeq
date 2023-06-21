##Create Venn diagrams for comparison##

## Create Venn diagrams for differently expressed genes results
#Load packages:
library(VennDiagram)
library(here)
library(tidyverse)

# Make the dir for the figures
dir.create("rnaseq2/results/venn/", recursive = TRUE, showWarnings = FALSE)

DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")

# Read the DE results
DE_res <- read_tsv(DE_files, id = "time", show_col_types = FALSE) %>%
  mutate(time = as.integer(gsub(".*DE_T(\\d\\d).tsv", "\\1", time))) %>%
  filter(isDE == TRUE) %>%
  mutate(DE_direction = ifelse(lfc > 0, "up", "down"),
         line = sub("_.*", "", group1))

# Venn diagram comparing Conrad & M92 treatment vs mock at 72 h
genes_A <- DE_res %>% filter(line == "Conrad", time == 72) %>% pull(gene)
genes_B <- DE_res %>% filter(line == "M92", time == 72) %>% pull(gene)

# If we only want upregulated
#genes_A <- DE_res %>% filter(line == "Conrad", time == 72, DE_direction == "up") %>% pull(gene)

#Make function
make_venn <- function(line1, line2, timepoint = 72, direction = "both") {

  if (direction == "both") {
    DE_directions <- c("up", "down")
  } else {
    DE_directions <- direction
  }

  # Get the gene sets
  genes_A <- DE_res %>%
    filter(line == line1, time == timepoint, DE_direction %in% DE_directions) %>%
    pull(gene)
  genes_B <- DE_res %>%
    filter(line == line2, time == timepoint, DE_direction %in% DE_directions) %>%
    pull(gene)
  venn_list <- list(genes_A, genes_B)
  names(venn_list) <- c(line1, line2)

  # Define the plot title
  plot_title <- paste0("Time: ", timepoint, "& DE direction:", direction)

  # Define the output file
  outfile <- paste0("rnaseq2/results/venn/", line1, "_", line2, "_", timepoint,
                    "_", direction, ".jpg")

  v1 <- venn.diagram(venn_list,
                     fill = c("lightblue", "purple"),
                     #fill = c("purple", "white"), #For the M92 & mutants comparison
                     main = plot_title,
                     alpha = c(0.5, 0.5),
                     filename=NULL,
                     cex = 2.20,
                     cat.cex = 3,
                     cat.dist = 0.07,
                     margin = 0.25)
  jpeg(outfile)
  grid.newpage()
  grid.draw(v1)
  dev.off()
}


make_venn(line1 = "Conrad", line2 = "M92", timepoint = 72, direction = "up")
make_venn(line1 = "Conrad", line2 = "M92", timepoint = 72, direction = "down")
make_venn(line1 = "Conrad", line2 = "M92", timepoint = 48, direction = "up")
make_venn(line1 = "Conrad", line2 = "M92", timepoint = 48, direction = "down")

make_venn(line1 = "M92", line2 = "FN0170228", timepoint = 48, direction = "both")
make_venn(line1 = "M92", line2 = "FN0170228", timepoint = 72, direction = "both")

make_venn(line1 = "M92", line2 = "FN0140856", timepoint = 48, direction = "both")
make_venn(line1 = "M92", line2 = "FN0140856", timepoint = 72, direction = "both")

make_venn(line1 = "M92", line2 = "FNMN0329", timepoint = 48, direction = "both")
make_venn(line1 = "M92", line2 = "FNMN0329", timepoint = 72, direction = "both")

