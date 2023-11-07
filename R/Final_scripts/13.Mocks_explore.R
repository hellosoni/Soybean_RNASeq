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
#Create output directories
outdir <- here("/fs/ess/PAS0471/nghi/rnaseq2/results/DE/Mocks")

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
#Rerun the Deseq
colData(dds)$group
dds$group <- factor(paste(dds$Line, dds$Treatment, sep = "_"))
table(dds$group)
levels(dds$group)

design(dds) <- formula(~ group)

# Differential expression analysis for mocks
dds <- DESeq(dds)
#The results table
res <- results(dds)
resultsNames(dds)

my_contrast <- c("M92_220_Mock","Sloan_Mock")
#my_contrast <- c("Conrad_Mock","Sloan_Mock")
res <- results(dds,
               contrast = c("group", my_contrast))

#How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm = TRUE)
#Convert res to regular dataframe
res <- res %>% as.data.frame()
#Upregulated genes
#put into seperate dataframe
upreg <- filter(res, res$padj < 0.05,res$log2FoldChange < 0)
downreg <- filter(res, res$padj < 0.05,res$log2FoldChange > 0)

upreg_genes <- row.names(upreg)
downreg <- row.names(downreg)


#Export the results
my_contrast
my_contrast_pasted<- paste0(my_contrast, collapse = "_vs_")
my_contrast_pasted

res$level1 <- my_contrast[1]
res$level2 <- my_contrast[2]

#save a result file
res_file <- file.path(outdir, paste0(my_contrast_pasted, '_res_mocks.txt'))
write.table(res,res_file ,
            sep = '\t', row.names = TRUE, quote = FALSE)

#save a separate table with only significant results
res_sig_file <- file.path(outdir, paste0(my_contrast_pasted, '_sig-res mocks.txt'))

res_sig <- res %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05)

write.table(res_sig, res_sig_file,
            sep = '\t', row.names = TRUE, quote = FALSE)
#save file for up and down regulated lists
up <-file.path(outdir, paste0(my_contrast_pasted, '_sig-res mocks up.txt'))
write.table(upreg,up,
            sep = '\t', row.names = TRUE, quote = FALSE)
down <-file.path(outdir, paste0(my_contrast_pasted, '_sig-res mocks down.txt'))
write.table(downreg,down,
            sep = '\t', row.names = TRUE, quote = FALSE)

