# SETUP ------------------------------------------------------------------------
# Packages
library(tidyverse)
library(DESeq2)

# Hardcoded settings
MAX_COUNTMEAN <- 1         # Max. gene count mean in FN line for a potentially deleted gene
MIN_LFC <- 3               # Min. LFC between M92 and FN line

#Note: Edit this script with correct genes

# Input files
dds_file <- "/fs/ess/PAS0471/nghi/rnaseq2/results/DE/dds.rds"

functions_script <- "/fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/rnaseq/rfuns/DE_funs.R"

# reported deleted genes list
FN0170228_deleted_file <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FN0170228_deleted.txt"
FN0140856_deleted_file <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FN0140856_deleted.txt"
FNMN0329_deleted_file <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FNMN0329_deleted.txt"

# Output files
outdir <- "/fs/ess/PAS0471/nghi/rnaseq2/results/deleted_genes"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Source script with functions
source(functions_script)


# PREP FOR ANALYSIS ------------------------------------------------------------
# Read the DESeq object
dds <- readRDS(dds_file)

# Remove specific samples
dds <- dds[, rownames(colData(dds)) != "ICR348_S81"]
dds <- dds[, rownames(colData(dds)) != "M6R172_S34"]
dds <- dds[, rownames(colData(dds)) != "I9R372_S104"]

# Make gene names compatible with Soybase names
rownames(dds) <- tolower(sub("\\.Wm82.*", "", rownames(dds)))

# Set the appropriate reference levels
colData(dds)$Line <- relevel(colData(dds)$Line, ref = "M92_220")
colData(dds)$Treatment <- relevel(colData(dds)$Treatment, ref = "Mock")

# Get a normalized, long-format count df with metadata
count_df <- norm_counts(dds, transform = 'lib_size')

# And separate matrix and df for heatmap
count_mat <- norm_counts(dds, transform = 'lib_size', return_matrix = TRUE)
meta <- as.data.frame(colData(dds))

# Run the DE analysis
design(dds) <- formula(~ Line)
dds <- DESeq(dds)


# FN0170228 --------------------------------------------------------------------
# Define line and deleted genes
FOCAL_LINE <- "FN0170228"
deleted <- readLines(FN0170228_deleted_file)

# Check whether all reported-to-be-deleted genes are found in our df
stopifnot(deleted %in% unique(count_df$gene))

# Heatmap of genes reported to be deleted
pheat(genes = deleted, count_mat = count_mat, meta_df = meta, groups = "Line")


# Get a df with DE results including group means
DE_res <- extract_DE(comp = c("M92_220", FOCAL_LINE),
                     fac = "Line",
                     dds = dds,
                     count_df = count_df)

# To find potentially deleted genes based on expression levels,
# select delection-candidate genes: almost no expression in FN line, and large-ish LFC
DE_deleted <- DE_res |>
  filter(mean_B <= MAX_COUNTMEAN, lfc > MIN_LFC) |>
  pull(gene)
message("Found ", length(DE_deleted), " potentially deleted genes")

# Heatmap of genes from our DE analysis
pheat(genes = DE_deleted, count_mat = count_mat, meta_df = meta, groups = "Line")
cheatmap(genes = DE_deleted, count_mat = count_mat, meta_df = meta, groups = "Line", scale = TRUE)

# Boxplots of reported deleted genes from our DE analysis
map(.x = deleted, .f = pbox,
    count_df = count_df, x_by = "Treatment", col_by = "Line")

# Check individual genes using boxplots
FN0170228_genes <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FN0170228_genes.txt"

genes <- readLines(FN0170228_genes)

map(.x = genes, .f = pbox,
    count_df = count_df, x_by = "Treatment", col_by = "Line")

# FN0140856 --------------------------------------------------------------------
# Define line and deleted genes
FOCAL_LINE <- "FN0140856"
deleted <- readLines(FN0140856_deleted_file)

# Check whether all reported-to-be-deleted genes are found in our df
stopifnot(deleted %in% unique(count_df$gene))

# Heatmap of genes reported to be deleted
pheat(genes = deleted, count_mat = count_mat, meta_df = meta, groups = "Line")

# Get a df with DE results including group means
DE_res <- extract_DE(comp = c("M92_220", FOCAL_LINE),
                     fac = "Line",
                     dds = dds,
                     count_df = count_df)

# To find potentially deleted genes based on expression levels,
# select delection-candidate genes: almost no expression in FN line, and large-ish LFC
DE_deleted <- DE_res |>
  filter(mean_B <= MAX_COUNTMEAN, lfc > MIN_LFC) |>
  pull(gene)
message("Found ", length(DE_deleted), " potentially deleted genes")

# Heatmap of genes from our DE analysis
pheat(genes = DE_deleted, count_mat = count_mat, meta_df = meta, groups = "Line")


# Boxplots of reported deleted genes from our DE analysis
map(.x = deleted, .f = pbox,
    count_df = count_df, x_by = "Treatment", col_by = "Line")

# Check individual genes using boxplots
FN0140856_genes <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FN0140856_genes.txt"

genes <- readLines(FN0140856_genes)

map(.x = genes, .f = pbox,
    count_df = count_df, x_by = "Treatment", col_by = "Line")


# FNMN0329 --------------------------------------------------------------------
# Define line and deleted genes
FOCAL_LINE <- "FNMN0329"  # Focal FN line
deleted <- readLines(FNMN0329_deleted_file)

# Check whether all reported-to-be-deleted genes are found in our df
stopifnot(deleted %in% unique(count_df$gene))

# Heatmap of genes reported to be deleted
pheat(genes = deleted, count_mat = count_mat, meta_df = meta, groups = "Line")

# Get a df with DE results including group means
DE_res <- extract_DE(comp = c("M92_220", FOCAL_LINE),
                     fac = "Line",
                     dds = dds,
                     count_df = count_df)

# To find potentially deleted genes based on expression levels,
# select delection-candidate genes: almost no expression in FN line, and large-ish LFC
DE_deleted <- DE_res |>
  filter(mean_B <= MAX_COUNTMEAN, lfc > MIN_LFC) |>
  pull(gene)
message("Found ", length(DE_deleted), " potentially deleted genes")

# Heatmap of genes from our DE analysis
pheat(genes = DE_deleted, count_mat = count_mat, meta_df = meta, groups = "Line")

# Check individual genes using boxplots
FNMN0329_genes <- "/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/FNMN0329_genes.txt"

genes <- readLines(FNMN0329_genes)

map(.x = genes, .f = pbox,
    count_df = count_df, x_by = "Treatment", col_by = "Line")

