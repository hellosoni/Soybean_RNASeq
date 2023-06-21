# SET UP -----------------------------------------------------------------------
# Load necessary packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here")            # Managing file paths
pacman::p_load(char = packages)

# Set a ggplot theme
theme_set(theme_bw())

# Define input files
count_table_file <- here("rnaseq2/results/nfc_rnaseq/star_salmon/salmon.merged.gene_counts_length_scaled.rds")
metadata_file <- here("rnaseq2/raw_data/meta/metadata.csv")

# Define output files
dds_file <- here("rnaseq2/results/dds.rds")
norm_counts_file <- "rnaseq2/results/normalized_counts.tsv"

# PREPARE COUNTS AND METADATA --------------------------------------------------
# Read the counts and metadata
counts <- readRDS(count_table_file)
counts <- round(assays(counts)$counts)
metadata <- read.csv(metadata_file, header = TRUE)

# Get rid of the "_S" suffices to sample names in the count matrix
colnames(counts) <- sub("_S\\d+$", "", colnames(counts))

# Check if the sample IDs of metadata and the count matrix match
stopifnot(all(colnames(counts) %in% metadata$Sample.ID))

# Turn the sample ID column into the row names
rownames(metadata) <- metadata$Sample.ID
metadata$Sample.ID <- NULL

# Order by sample ID, or DESeq will complain
metadata <- metadata[order(rownames(metadata)), ]
counts <- counts[, order(colnames(counts))]

# Change "Line", "Time", and "Treatment" into factors
metadata <- metadata %>%
  mutate(Line = as.factor(Line),
         Time = as.factor(Time),
         Treatment = as.factor(Treatment))


# CREATE THE DESEQ OBJECT ------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ 1)

# Remove samples with low counts
dds <- dds[, rownames(colData(dds)) != "ICR348"]
dds <- dds[, rownames(colData(dds)) != "M6R172"]
dds <- dds[, rownames(colData(dds)) != "I9R372"]

# Merge factors Line & treatment into a single factor called "group"
dds$group <- factor(paste(dds$Line, dds$Treatment, sep = "_"))

# Write to file
saveRDS(dds, dds_file)


# GET NORMALIZED COUNTS --------------------------------------------------------
# Normalization using median of ratios method

# generate size factors
dds <- estimateSizeFactors(dds)

#take a look at the normalization factor applied to each sample
sizeFactors(dds)

#Retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE
normalized_counts <- counts(dds, normalized=TRUE)

#save this normalized data matrix to file for later use in visualization
write.table(normalized_counts, norm_counts_file,
            sep = "\t", quote = FALSE, col.names = NA)

#NOTE: DESeq2 doesn't actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM).
#These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential expression analysis which use the negative binomial model.

