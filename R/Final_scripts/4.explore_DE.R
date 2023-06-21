#Explore DE results#

# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)

#Loading DE results
DE_dir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE"
DE_files <- list.files(DE_dir, full.names = TRUE, pattern = "DE_")

# Read the DE results
DE_res <- read_tsv(DE_files, id = "time", show_col_types = FALSE) %>%
  mutate(time = as.integer(gsub(".*DE_T(\\d\\d).tsv", "\\1", time)))

# Check the number of DE genes by contrast, timepoint and DE direction
DE_res %>%
  filter(isDE == TRUE) %>%
  count(time, contrast, lfc > 0) %>%
  mutate(DE_direction = ifelse(`lfc > 0` == TRUE, "up", "down")) %>%
  select(-`lfc > 0`)

#Put results into a dataframe
results <- DE_res %>%
  filter(isDE == TRUE) %>%
  count(time, contrast, lfc > 0) %>%
  mutate(DE_direction = ifelse(`lfc > 0` == TRUE, "up", "down")) %>%
  select(-`lfc > 0`)


#Visually explore the DE results with volcano plots and plotting specific genes-----

#Volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



#Plot specific genes
#select the 5 genes with the lowest adjusted p-value
top5 <- row.names(res[order(res$padj)[1:5], ])
#create a function to make a plot for a single gene:

plotgene <- function(geneID, dds24) {

  d <- plotCounts(dds24,
                  gene = geneID,
                  intgroup = "group",
                  returnData = TRUE)

  p <- ggplot(d, aes(x = group, y = count)) +
    geom_point(position = position_jitter(w = 0.1, h = 0)) +
    labs(title = geneID) +
    theme_bw()

  print(p)
}

none <- sapply(top5, plotgene, dds24)
