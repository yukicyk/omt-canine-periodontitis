#
# SCRIPT: 02_dada2_and_phyloseq.R
#
# PURPOSE: Core script for processing trimmed reads. It performs quality filtering,
#          error rate learning, ASV inference, merging, and chimera removal.
#          It then assigns taxonomy and combines all data into a phyloseq object.
#
# INPUTS:
#   - Trimmed FASTQ files from `01_preprocessing.sh` (e.g., in `data/trimmed_fastq/`)
#   - `data/metadata.txt`: Sample metadata file.
#   - `data/taxonomy/`: Folder with custom taxonomy reference files.
#   - `data/ASV_phylogeny.nwk`: Phylogenetic tree built externally.
#
# OUTPUTS:
#   - `results/phyloseq_object.rds`: The final, comprehensive phyloseq object.
#   - `results/dada2_seqtab.rds`: The final ASV table.
#   - `results/dada2_taxa.rds`: The final taxonomy table.
#

# --- 1. Load Libraries and Set Paths ---

# It's good practice to install/load all packages at the beginning.
# renv will manage versions automatically.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dada2, phyloseq, ggplot2, seqinr, ape, BiocManager)

# Define relative paths for clarity and reproducibility
fastq_path <- "data/trimmed_fastq" # Path to your trimmed .fastq.gz files
filt_path <- file.path(fastq_path, "filtered") # Subdirectory for filtered files
tax_path <- "data/taxonomy" # Path to taxonomy reference files
results_path <- "results"

# Create output directories if they don't exist
if (!dir.exists(filt_path)) dir.create(filt_path)
if (!dir.exists(results_path)) dir.create(results_path)

# --- 2. Quality Filtering and Trimming ---

# Get matched lists of forward and reverse fastq files
fnFs <- sort(list.files(fastq_path, pattern = "_1.fq.gz.out", full.names = TRUE))
fnRs <- sort(list.files(fastq_path, pattern = "_2.fq.gz.out", full.names = TRUE))

# Extract sample names, assuming format: SAMPLENAME_1.fq.gz.out
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Define filenames for filtered fastq files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Plot quality profiles to inform trimming parameters
# This saves all plots into one PDF for easier review
pdf(file.path(results_path, "figures/quality_profiles.pdf"))
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
dev.off()

# Perform filtering and trimming
# NOTE: truncLen values were chosen based on the quality plots.
# They should be adjusted if running on new data.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(250, 230),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# --- 3. DADA2 Core Algorithm ---

# Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot errors to verify
pdf(file.path(results_path, "figures/dada2_error_plots.pdf"))
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer sample sequences (the core DADA2 algorithm)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Construct sequence table (ASV table)
seqtab <- makeSequenceTable(mergers)

# Filter sequences to expected length range for V3-V4 region
seqtab_filt <- seqtab[, nchar(colnames(seqtab)) %in% 400:452]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab_filt, method = "consensus", multithread = TRUE, verbose = TRUE)
cat("Fraction of non-chimeric reads:", sum(seqtab.nochim) / sum(seqtab_filt), "\n")

# --- 4. Taxonomy Assignment ---

# Assign taxonomy using the custom reference database
# NOTE: Ensure these reference files are in the `data/taxonomy/` directory
taxa <- assignTaxonomy(seqtab.nochim, file.path(tax_path, "COT_HOMD_RDP_train_set_v2.txt"), multithread = TRUE)
taxa <- addSpecies(taxa, file.path(tax_path, "COT_HOMD_RDP_spec_assign_v2.txt"), allowMultiple = FALSE) # Using single assignment as in original analysis

# --- 5. Create Phyloseq Object ---

# Import metadata
samdf <- read.delim(file.path("data", "metadata.txt"), row.names = 1)

# Import phylogenetic tree
phy_tree <- read.tree(file.path("data", "ASV_phylogeny.nwk"))

# Create the phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa),
  phy_tree(phy_tree)
)

# Remove the obsolete sample found in the original analysis
ps <- subset_samples(ps, sample_names(ps) != "MT-E-V-2")

# --- 6. Save Outputs ---

# Save the final phyloseq object and intermediate files
saveRDS(ps, file = file.path(results_path, "phyloseq_object.rds"))
saveRDS(seqtab.nochim, file = file.path(results_path, "dada2_seqtab.rds"))
saveRDS(taxa, file = file.path(results_path, "dada2_taxa.rds"))

print("DADA2 and Phyloseq object creation complete.")
print(paste("Final phyloseq object saved to:", file.path(results_path, "phyloseq_object.rds")))