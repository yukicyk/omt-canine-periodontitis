# Data Directory README

This directory contains all the necessary input files to run the analysis pipeline, from raw data processing to final statistical analysis.

---

### 1. Raw Sequencing Data

The raw, demultiplexed, paired-end FASTQ files are the primary input for this project.

-   **Availability**: Raw sequencing data has been deposited in the NCBI Sequence Read Archive (SRA) and is publicly available.
    -   **Accession Number**: `PRJNA598540`
    -   **Direct Link**: [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA598540]

-   **Quality Control**: The quality of the raw reads was assessed prior to analysis. The quality score plots for both forward (R1) and reverse (R2) reads are available at:
    -   `data/raw_fastq/Canine_transplant_MiSeq_R1_QS.pdf`
    -   `data/raw_fastq/Canine_transplant_MiSeq_R2_QS.pdf`

-   **Preprocessing**: The exact commands used for quality trimming and filtering of the raw FASTQ files are documented in the `R/01_preprocessing.sh` script. To fully reproduce the analysis from scratch, download the raw data from the SRA and execute this script. The output will be the trimmed FASTQ files used in the DADA2 pipeline.

---

### 2. Core Input Files for Analysis

These files are required to run the main DADA2 pipeline (`R/02_dada2_and_phyloseq.R`) and subsequent downstream analysis scripts.

-   **`data/metadata.txt`**: The primary sample metadata file in **tab-separated format**. It maps each `SampleID` to its corresponding experimental variables, such as `Dog`, `Treatment`, and `Timepoint`.

-   **`data/clinical_metadata.csv`**: A supplementary **comma-separated** file containing the clinical parameters measured for the samples, including `PPD`, `BOP`, and `PlaqueIndex`. This file is used in the statistical analysis script (`R/04_statistical_and_clinical_analysis.R`).

-   **`data/ASV_phylogeny.nwk`**: A phylogenetic tree of the final Amplicon Sequence Variants (ASVs) in **Newick format**. This tree was built externally after aligning the ASV sequences and is required for phylogenetic diversity calculations (e.g., UniFrac).

---

### 3. Taxonomy Reference Files

This section describes the reference files used for taxonomic classification.

#### Original Project Methodology

The original analysis, as detailed in `docs/taxonomy.md`, used a custom database created by combining three sources:
1.  **RDP trainset 16** (Cole et al., 2014)
2.  **eHOMD 16S rRNA RefSeq v15.1** (Chen et al., 2010)
3.  **Canine Oral Microbiome (COT)** database (Dewhirst et al., 2012)

The `R/00_setup_taxonomy_reference.R` script processes these source files from `data/raw_taxonomy/` into DADA2-compatible formats in `data/taxonomy/`.

#### Recommendation for Modern Reproducibility

The RDP and Greengenes databases are now largely considered outdated. For any new or reproduced analysis, we strongly recommend using the **SILVA database**, which is the current standard in microbiome research due to its regular updates and high-quality curation.

-   **Recommended Database**: SILVA (e.g., v138.1 or newer)
-   **Download Location**: Pre-formatted versions for DADA2 are available here: [https://benjjneb.github.io/dada2/training.html](https://benjjneb.github.io/dada2/training.html)
-   **To Reproduce**: A user should download the SILVA training set (`silva_*_train_set.fa.gz`) and the species assignment file (`silva_species_assignment_*.fa.gz`) and place them in the `data/taxonomy/` or `data/raw_taxonomy/` directory. The `R/02_dada2_and_phyloseq.R` script can then be pointed to these new files.

**Note on combining databases**: For the most robust analysis of this specific dataset, the ideal approach would be to combine the modern SILVA database with the specialized eHOMD and COT databases that contains most the relevant reference oral taxon from human and canine sources. This advanced task requires careful formatting of the taxonomy strings to ensure they are compatible across all files.