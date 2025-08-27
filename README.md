
# Oral Microbiota Transplant in Dogs with Naturally Occurring Periodontitis

This repository contains the complete analysis workflow for the study "Oral Microbiota Transplant in Dogs with Naturally Occurring Periodontitis," published in the *Journal of Dental Research*.

**Publication:** Beikler, T., Bunte, K., Chan, Y., et al. (2021). Oral Microbiota Transplant in Dogs with Naturally Occurring Periodontitis. *Journal of Dental Research*, 100(7), 764–770. [https://doi.org/10.1177/0022034521995423](https://doi.org/10.1177/0022034521995423) 

**Data Availability:** Raw sequencing data is available at the NCBI Sequence Read Archive (SRA) under accession number **PRJNA598540**.

---

### Project Overview

This study evaluated the safety and efficacy of an oral microbiota transplant (OMT) for treating naturally occurring periodontitis in beagle dogs. The analysis involved 16S rRNA gene amplicon sequencing to track changes in the oral microbiota of test (OMT) and control groups over 12 weeks. The workflow encompasses sequence data processing, taxonomic classification, diversity analysis, statistical testing, and correlation with clinical parameters like Probing Pocket Depth (PPD) and Bleeding on Probing (BOP).

---

### Compliance and Data Integrity (GLP, ALCOA+, ISO, GDPR)

This repository is designed to adhere to best practices for scientific reproducibility and data integrity, in line with several key standards:

*   **ALCOA+ (Attributable, Legible, Contemporaneous, Original, Accurate, Complete, Consistent, Enduring, Available):**
    *   **Attributable:** All changes to code are tracked via Git commits, linking work to the author.
    *   **Legible & Enduring:** All code and documentation are in human-readable text formats and stored permanently in this version-controlled repository.
    *   **Original:** Git preserves the original data and all subsequent changes, providing a complete audit trail.
    *   **Accurate & Consistent:** The documented, sequential workflow ensures that the analysis can be reproduced accurately and consistently.
    *   **Complete & Available:** The repository contains all scripts, documentation, and environment information necessary to replicate the analysis.

*   **GLP (Good Laboratory Practice) & ISO 9001/13485:** The structured workflow, version-controlled scripts, and detailed documentation serve as a standard operating procedure (SOP) for the data analysis. The use of `renv` for package management ensures a controlled and validated software environment, a key aspect of quality management and process control.

*   **GDPR / UK DPA:** All data used in this public repository is fully anonymized. The sample IDs are codes that do not contain personal identifying information. The raw sequence data, which contains no personal information, is hosted on a public archive (NCBI SRA).

---

### System and Software Requirements

*   **R:** Version 3.4.2 or later.
*   **R Packages:** A comprehensive list of required packages is provided in the `renv.lock` file. Key packages include:
    *   `dada2`
    *   `phyloseq`
    *   `ggplot2`
    *   `dplyr`
    *   `tidyr`
    *   `vegan` (for PERMANOVA and other ecological stats)
    *   `DESeq2`
    *   `randomForest`
    *   `Hmisc`
*   **External Tools:**
    *   `cutadapt` (v1.14 or later) for primer trimming.
    *   `vsearch` (for PICRUSt pre-processing).
    *   `mafft` (for multiple sequence alignment to build the phylogenetic tree).

#### Environment Setup with `renv`

To ensure perfect reproducibility, this project uses the `renv` package to lock package versions. To install all required packages at their correct versions, run the following command in R within the project's root directory:

```R
# Run this once to install all required packages
renv::restore()
```

---

### Proposed Repository Structure

```
.
├── R/
│   ├── 00_setup_taxonomy_reference.R
│   ├── 01_preprocessing.sh
│   ├── 02_dada2_and_phyloseq.R
│   ├── 03_diversity_and_ordination_analysis.R
│   ├── 04_statistical_and_clinical_analysis.R
│   ├── 05_functional_prediction_analysis.R
│   └── 06_animated_visualizations.R   <-- NEW SCRIPT
├── data/
│   ├── raw_fastq/              <- fastq QC graphs
│   ├── data_README.md
│   ├── metadata.txt
│   ├── ASV_phylogeny.nwk
│   └── clinical_metadata.csv
├── results/
│   ├── figures/
│   └── tables/
├── docs/
│   ├── analysis_methods_explained.md
│   ├── analysis_narrative.md
│   ├── taxonomy.md
│   ├── OMT_Dog_JDR2021.pdf
│   └── SI_OMT_Dog_JDR2021.pdf
├── archive/
│   ├── (All original, unedited logs, scripts and reports)
├── .gitignore
├── Canine_Periodontitis_OMT.Rproj
├── README.md
└── renv.lock
```

---

### Analysis Workflow: Step-by-Step Guide

Follow these steps in order to reproduce the analysis from the raw data.

#### **Step 0: Setup Taxonomy Reference Database**
*   **Script:** `R/00_setup_taxonomy_reference.R`
*   **Source Files:** `tidying_up_taxonomy_ref.R`, `tidying_up_taxonomy_species_ref.R`
*   **Description:** This script prepares the custom taxonomy reference files from the RDP, the Canine Oral Microbiome (COT) and Human Oral Microbiome Database (HOMD). This is a one-time setup step required before running the main DADA2 pipeline. See `docs/taxonomy.md` for details. 
*   **Inputs:** Raw FASTA and taxonomy files from HOMD/COT.
*   **Outputs:** Formatted `COT_HOMD_RDP_train_set.fa` and `COT_HOMD_RDP_spec_assign_v2.txt` files for DADA2.

#### **Step 1: Raw Sequence Pre-processing**
*   **Script:** `R/01_preprocessing.sh`
*   **Source Files:** `Canine_MiSeq_LOG copy.txt`, `Cutadapt.log.txt`
*   **Description:** This shell script uses `cutadapt` to trim the 341F and 806R 16S rRNA gene primers from the raw, paired-end `.fq.gz` files. See `data/data_README.md` for details.
*   **Inputs:** Raw `.fq.gz` files for forward and reverse reads.
*   **Outputs:** Trimmed `.fq.gz.out` files for each sample, ready for DADA2.

```bash
# Example commands from the script:
# Trim forward reads (R1)
find . -name "*_1.fq.gz" -exec cutadapt -g ACTCCTACGGGAGGCAGCAG -o {}.out {} \;

# Trim reverse reads (R2)
find . -name "*_2.fq.gz" -exec cutadapt -g GGACTACHVGGGTWTCTAAT -o {}.out {} \;
```

#### **Step 2: ASV Calling, Taxonomy Assignment, and Phyloseq Object Creation**
*   **Script:** `R/02_dada2_and_phyloseq.R`
*   **Source Files:** `Canine_transplant_MiSeq_v2.R`
*   **Description:** This is the core script for processing the trimmed reads. It performs quality filtering, error rate learning, dereplication, ASV inference, merging of paired-end reads, and chimera removal. It then assigns taxonomy using the custom reference files and combines the ASV table, taxonomy table, metadata, and phylogenetic tree into a single `phyloseq` object.
*   **Inputs:**
    *   Trimmed FASTQ files from Step 1.
    *   `mappingV2.txt` metadata file.
    *   Taxonomy reference files from Step 0.
    *   A phylogenetic tree (`ASV.nwk`), which was built externally using `mafft` and `MEGA7`.
*   **Outputs:**
    *   `seqtab_final.rds`: The final ASV table.
    *   `taxa_final.rds`: The final taxonomy table.
    *   `phyloseq_object.rds`: The final, comprehensive phyloseq object for downstream analysis.

#### **Step 3: Diversity and Ordination Analysis**
*   **Script:** `R/03_diversity_and_ordination_analysis.R`
*   **Source Files:** `canine_analysisv2.R`
*   **Description:** This script handles the core ecological analyses presented in the paper. It calculates alpha diversity metrics (Shannon, Observed ASVs) and performs beta diversity analysis using Bray-Curtis and weighted/unweighted UniFrac distances. It generates ordination plots (PCoA, NMDS) to visualize community differences (Figures 2 & 3 in the paper) and performs PERMANOVA tests to assess statistical significance. For a detailed explanation of PERMANOVA and the associated power analysis, see `docs/analysis_methods_explained.md`.
*   **Inputs:** `phyloseq_object.rds`.
*   **Outputs:** Alpha and beta diversity plots, ordination plots, and statistical results tables.

#### **Step 4: Advanced Statistical and Clinical Analysis**
*   **Script:** `R/04_statistical_and_clinical_analysis.R`
*   **Source Files:** `PRC.pdf`, `power.pdf`, `correlationAnalysis.R`, `cforestAnalysisWithMeta.R`
*   **Description:** This script brings together several advanced analyses:
    *   **Principal Response Curve (PRC):** To visualize the effect of the OMT on the community over time relative to the control (Figure 4 in the paper). For a detailed theoretical explanation of the PRC method, see `docs/analysis_methods_explained.md`.
    *   **PERMANOVA Power Analysis:** A retrospective analysis to confirm the statistical power of the study design.
    *   **Correlation Analysis:** Uses Pearson correlation to find significant associations between microbial taxa (ASVs/Species) and clinical parameters (BOP, PPD).
    *   **Predictive Modeling:** Employs Random Forest to identify which taxa are most predictive of clinical outcomes.
*   **Inputs:** `phyloseq_object.rds`, clinical data files (`Ergebnisse_BOP_PPD_summarized_BOP_PPD_3.txt`).
*   **Outputs:** PRC plots, correlation heatmaps, variable importance plots from Random Forest, and associated statistical tables.

#### **Step 5: Functional Prediction Analysis**
*   **Script:** `R/05_functional_prediction_analysis.R`
*   **Source Files:** `dada2_picrust.R`, `koAnalysis.R`
*   **Description:** This script details the workflow for predicting the functional potential of the microbial communities. It first prepares the DADA2 output for use with PICRUSt (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States). It then analyzes the resulting KEGG Orthology (KO) abundance tables to find functional pathways that differ between groups and timepoints.
*   **Inputs:** `seqtab_final.rds`, Greengenes reference database.
*   **Outputs:** Predicted KEGG pathway abundance tables and plots comparing functional profiles.

#### **Step 6: Create Animated Visualizations (Optional)**
*   **Script:** `R/06_animated_visualizations.R`
*   **Description:** This script uses the `gganimate` package to create animated ordination plots (MDS on Bray-Curtis dissimilarity). These animations visualize the trajectory of the microbial community for each dog over the 12-week study period, providing an intuitive view of community shifts in the control and recipient groups.
*   **Inputs:** `phyloseq_object.rds`.
*   **Outputs:** Animated `.gif` files (e.g., `control_group_shift.gif`, `recipient_group_shift.gif`) saved to the `results/figures/` directory.

---

### Data Availability and Use of Synthetic Data for Reproducibility

**Primary Data:** The raw 16S rRNA gene sequencing data is publicly available from the NCBI Sequence Read Archive (SRA) under accession number **PRJNA598540**.

**Metadata and Clinical Data:** To ensure full computational reproducibility, this repository includes all metadata required to execute the analysis workflow from start to finish. The data is provided in two key files:

  - `data/metadata.txt`: Links sample IDs to their respective treatment groups, dog IDs, and timepoints.
  - `data/clinical_metadata.csv`: Provides the deanonymized clinical parameters (e.g., mean PPD, BOP% and Plaque%) for each sample and timepoint.

All data in this repository has been fully deanonymized and contains the same information on experimental groups and clinical outcomes that was used for the analyses presented in the publication.

For compliance with the approved animal research ethics protocol **(Landesamt für Natur und Verbraucherschutz, #84-02.04.2014.A449)** and institutional data management policies, the original, non-anonymized source records and primary clinical charts are securely archived at the research institution.

This approach ensures that the computational methods are transparent and fully executable, respecting data governance policies while upholding the principles of open and reproducible science (e.g., FAIR, ALCOA+).
=======
