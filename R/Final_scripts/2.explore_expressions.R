# SET UP -----------------------------------------------------------------------
# Load necessary packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "ggrepel")
pacman::p_load(char = packages)

# Define input files
dds_file <- "rnaseq2/results/DE/dds.rds"

# Read the DESeq object
dds <- readRDS(dds_file)


# EXPLORE GENE COUNTS  ---------------------------------------------------------
##Explore the count data
counts_dds <- assay(dds)
dim(counts_dds) #What are number of rows and columns of the count matrix?
dim(counts_dds[rowSums(counts_dds) > 0, ]) #How many genes have non-zero counts?
dim(counts_dds[rowSums(counts_dds) >= 10, ]) #How many genes have total counts of at least 10?

#Histogram of gene count
theme_set(theme_bw())

summed_gene_counts <- data.frame(gene_count = rowSums(counts_dds)) %>%  rownames_to_column("gene_id")

ggplot(data = summed_gene_counts) +
  geom_histogram(aes(x = gene_count), binwidth = 10000) +
  scale_y_log10(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0))

ggplot(data = summed_gene_counts) +
  geom_histogram(aes(x = gene_count), binwidth = 1000) +
  scale_y_log10(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 200000), expand = c(0,0)) +
  theme(plot.margin = margin(0.5, 0.7, 0.5, 0.5, "cm"))

#sum of counts for each column
apply(X = counts_dds, MARGIN = 2, FUN = sum)


# PCA --------------------------------------------------------------------------
#Run the PCA and prepare for plotting
vsd <- varianceStabilizingTransformation(dds, blind = TRUE) # normalize the count data to have even sampling across samples (with respect to library size) and approximately even variance
#run the PCA and retrieve the data to plot with ggplot2
pcaData <- plotPCA(vsd,
                   ntop = 500,
                   intgroup = c("Line","Time", "Treatment" ),
                   returnData = TRUE)
#extract the percentage of variance explained by different principal components, so we can later add this information to the plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar
#create a plot title with the species name in italic using the  expression() function
plot_title <- expression("PCA of " * italic(Glycine ~ max) * " transcript count data")

##Plot the PCA results------------------------
ggplot(pcaData,
       aes(x = PC1, y = PC2, color = Line, shape = Treatment)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title)

pca_plot2 <- ggplot(pcaData,
                    aes(PC1, PC2, color = Line, shape = Treatment)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle(plot_title)

print(pca_plot2)

