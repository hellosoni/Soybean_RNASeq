## Workkflow rnaseq 2 ###

##Download the data
sftp -P 22 fiorella_cisneros_carter@3.17.199.203 # password: XDwc4vgiae*yf0k
get -R 221123_Meulia_GSL-FCC-3127

## Get the reference genome files
cp ../01_archive/2021-10_nghi/data/ref/jgi/phytozome/Gmax/Wm82.a4.v1/assembly/Gmax_508_v4.0.fa data/ref
cp ../01_archive/2021-10_nghi/data/ref/jgi/phytozome/Gmax/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 data/ref

## Processing raw reads
# Concatenate FASTQ files from different lanes 
indir=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/221123_Meulia_GSL-FCC-3127/
outdir=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/fastq/concat
sbatch /fs/project/PAS0471/nghi/rnaseq/mcic-scripts/misc/fqconcat.sh -i "$indir" -o "$outdir"

## Running fastqc for the concatenated files: (No need for the pipeline)
indir=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/fastq/concat
outdir=/fs/ess/PAS0471/nghi/rnaseq2/results/qc_concat
for infile in $in_dir/*.fastq.gz; do 
    sbatch /fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/qc/fastqc.sh -i "$infile" -o "$outdir"
done

#Subset the fastq files (To test the pipeline) using fqsub.sh
R1_in=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/fastq/concat
outdir=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/subset
for infile in $R1_in/I8R1*_R1_001.fastq.gz; do
    sbatch /fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/utils/fqsub.sh  -i  "$infile" -o "$outdir"
done

#Running nextflow pipeline for the subset files
samplesheet=/fs/ess/PAS0471/nghi/rnaseq2/samplesheet.csv
outdir=/fs/ess/PAS0471/nghi/rnaseq2/results/subset
ref_fasta=/fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/assembly/Gmax_508_v4.0.fa
ref_annot=/fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff
#Note: can directly put in the sbatch
sbatch /fs/ess/PAS0471/nghi/rnaseq2/mcic-scripts/rnaseq/nfcore_rnaseq.sh --samplesheet /fs/ess/PAS0471/nghi/rnaseq2/samplesheet.csv --ref_fasta /fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/assembly/Gmax_508_v4.0.fa --ref_annot /fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff -o /fs/ess/PAS0471/nghi/rnaseq2/results/subset

# RUN THE NFCORE RNASEQ PIPELINE (for all the sample files)
# ==============================================================================
# Make output directories
mkdir -p scripts /fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta

# Define input files
ref_fasta=/fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/assembly/Gmax_508_v4.0.fa
ref_gtf=/fs/ess/PAS0471/nghi/rnaseq2/ref_genomes/A4/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff
fqdir=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/fastq/concat
samplesheet=/fs/ess/PAS0471/nghi/rnaseq2/raw_data/meta/samplesheet.csv   # This file will be produced by the script below



# Download the samplesheet prep script
wget -P scripts -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py

# Prep the samplesheet
python3 /fs/project/PAS0471/nghi/rnaseq2/raw_data/meta/scripts/fastq_dir_to_samplesheet.py \
    "$fqdir" \
    "$samplesheet" \
    --strandedness reverse \
    --read1_extension "_R1_001.fastq.gz" \
    --read2_extension "_R2_001.fastq.gz"
cat $samplesheet # Check if it looks good

# Run the pipeline
sbatch mcic-scripts/rnaseq/nfcore_rnaseq.sh \
    --samplesheet "$samplesheet" \
    --ref_fasta "$ref_fasta" \
    --ref_annot "$ref_gtf" \
    -o results/nfc_rnaseq \
    --more_args "--skip_markduplicates"
    

# You might want to skip this step, since it takes long time, and may time out the Slurm job:
#--more_args "--skip_markduplicates"

# Process the gene counts from the workflow for DE analysis etc
#> See the script mcic-scripts/rnaseq/R_templates/nfcore-rnaseq_load-counts.R

