## DE ANALYSIS ##

# SET UP -----------------------------------------------------------------------
# Load necessary packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("DESeq2",          # Differential expression analysis
              "tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "apeglm")          # For LFC shrinkage
pacman::p_load(char = packages)

# Define input files
dds_file <- "rnaseq2/results/dds.rds"
functions_script <- "rnaseq2/mcic-scripts/rnaseq/rfuns/DE_funs.R"

# Define contrasts to analyze
my_timepoints <- c( 24,48, 72)
my_contrasts <- list(
  c("Conrad_Inoculated","Conrad_Mock"),
  c("M92_220_Inoculated","M92_220_Mock"),
  c("Sloan_Inoculated", "Sloan_Mock"),
  c("FN0140856_Inoculated", "FN0140856_Mock"),
  c("FN0170228_Inoculated", "FN0170228_Mock" ),
  c("FNMN0329_Inoculated", "FNMN0329_Mock")
)

# Define and create the output directories
outdir <- here("rnaseq2/results/DE")
plotdir <- here("rnaseq2/results/DE/fig")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

# Source script with functions
source(functions_script)

# Prep the DESeq object
dds <- readRDS(dds_file)
design(dds) <- formula(~ group)

# Function to run the DE analysis for 1 timepoint
run_DE <- function(dds,
                   timepoint,
                   contrasts,
                   norm_transform = "rlog",
                   lfc_tres = 0) {
  message("\n-------------\nTime point: ", timepoint)

  # Subset dds object to 1 timepoint
  dds <- subset(dds, select=colData(dds)$Time == timepoint)

  # Just for timepoint 24hpi, remove M92_220 because there are no mocks
  if (timepoint == 24) {
    # Remove M92 from the DESeq object
    dds <- subset(dds, select=colData(dds)$Line != "M92_220")
    colData(dds)$group <- droplevels(colData(dds)$group)
    # Remove M92 from the list of contrasts
    contrasts <- contrasts[-grep("M92_220", contrasts)]
  }

  # Get a dataframe with counts
  count_df <- norm_counts(dds, transform = norm_transform)

  # Run the DE analysis
  dds <- suppressMessages(DESeq(dds))

  # Extract the DE results -- for each contrast
  DE_res <- map_dfr(.x = contrasts, .f = function(contrast) {
    extract_DE(comp = contrast, fac = "group", dds = dds, count_df = count_df,
               lfc_tres = lfc_tres)
  })

  # Check nr's of DE genes, incl. up- and downregulated
  DE_res %>%
    filter(isDE == TRUE) %>%
    count(contrast, lfc < 0)

  # Write to file -- one for all contrasts
  DE_res_file <- here(outdir, paste0("DE_T", timepoint, ".tsv"))
  write_tsv(DE_res, DE_res_file)
}


# RUN THE DE ANALYSIS ----------------------------------------------------------
walk(.x = my_timepoints, .f = function(timepoint) {
  run_DE(dds = dds, timepoint = timepoint, contrasts = my_contrasts)
})




