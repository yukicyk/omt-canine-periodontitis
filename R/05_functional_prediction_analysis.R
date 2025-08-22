#
# SCRIPT: 05_functional_prediction_analysis.R
#
# PURPOSE: Predict functional potential using PICRUSt and analyze the resulting
#          KEGG Orthology (KO) abundance tables.
#
# NOTE: This workflow uses the original PICRUSt (v1). For new analyses, the
#       updated PICRUSt2 is recommended, which has a different workflow.
#
# Part 1: Prepares DADA2 output for PICRUSt (run in a specific conda env).
# Part 2: Analyzes the predicted functional data in R.
#
# INPUTS:
#   - `results/dada2_seqtab.rds`: ASV table from script 02.
#   - `data/metadata.txt`: Sample metadata file.
#   - Greengenes reference files (required by PICRUSt).
#   - PICRUSt output files (e.g., `predicted_metagenomes.L2.txt`).
#
# OUTPUTS:
#   - Tables of significantly different KEGG pathways.
#   - Plots comparing functional profiles.
#

# --- 1. Load Libraries ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, biomformat, tibble)

# Define output paths
fig_path <- "results/figures"
table_path <- "results/tables"

# --- 2. Part 1: Preparing Data for PICRUSt (Conceptual Workflow) ---

# This part of the original script requires external command-line tools (vsearch, PICRUSt scripts)
# and is best run within a dedicated conda environment. The R code below is for context.

# conceptual_prepare_for_picrust <- function() {
#   library(dada2)
#   library(seqinr)
#   library(biom)
#
#   # 1. Load sequence table
#   seqtab.nochim <- readRDS("results/dada2_seqtab.rds")
#
#   # 2. Format for PICRUSt (assign ASV IDs, create FASTA)
#   seqs <- getSequences(seqtab.nochim)
#   ids <- paste0("ASV_", 1:length(seqs))
#   fasta_out <- "picrust_input/asv.fasta"
#   write.fasta(as.list(seqs), ids, fasta_out)
#
#   # 3. Create BIOM table
#   # The original PICRUSt requires OTU IDs that match Greengenes. The DADA2-to-PICRUSt
#   # workflow involves mapping ASVs to a reference (like Greengenes 97% OTUs).
#   # This step is complex and depends on the PICRUSt version.
#
#   # 4. Run PICRUSt command-line tools outside R:
#   #    - `normalize_by_copy_number.py`
#   #    - `predict_metagenomes.py`
#   #    - `categorize_by_function.py`
# }

# --- 3. Part 2: Analyzing PICRUSt Output ---

# Load PICRUSt output files (assuming they are in a `results/picrust/` directory)
picrust_path <- "results/picrust"
l2_data <- read_tsv(file.path(picrust_path, "predicted_metagenomes.L2.txt"), skip = 1)
l3_data <- read_tsv(file.path(picrust_path, "predicted_metagenomes.L3.txt"), skip = 1)
map_data <- read.delim("data/metadata.txt")

# --- Data Tidying Function ---
# This function tidies the raw PICRUSt output and merges it with metadata
tidy_picrust_output <- function(df, metadata) {
  df %>%
    rename(KO_Category = `#OTU ID`) %>%
    # Pivot to long format
    pivot_longer(
      cols = -c(KO_Category, KEGG_Pathways),
      names_to = "SampleID",
      values_to = "Abundance"
    ) %>%
    # Merge with sample metadata
    left_join(metadata, by = "SampleID")
}

l2_tidy <- tidy_picrust_output(l2_data, map_data)
l3_tidy <- tidy_picrust_output(l3_data, map_data)


# --- Statistical Analysis Function ---
# This function performs paired t-tests between two timepoints for a given treatment group.
# This avoids the massive code repetition in the original script.
run_paired_analysis <- function(df, treat_group, time1, time2) {
  df %>%
    filter(Treatment == treat_group, Timepoint %in% c(time1, time2)) %>%
    # Ensure each dog has data for both timepoints
    group_by(Dog, KO_Category) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    # Perform paired t-test for each KO category
    group_by(KO_Category, KEGG_Pathways) %>%
    summarise(
      p_value = t.test(Abundance ~ Timepoint, paired = TRUE, data = .)$p.value,
      mean_t1 = mean(Abundance[Timepoint == time1]),
      mean_t2 = mean(Abundance[Timepoint == time2]),
      .groups = 'drop'
    ) %>%
    mutate(
      FDR = p.adjust(p_value, method = "fdr"),
      Comparison = paste(treat_group, time1, "vs", time2),
      Change = ifelse(mean_t2 > mean_t1, "Increase", "Decrease")
    ) %>%
    select(Comparison, KO_Category, KEGG_Pathways, Change, FDR) %>%
    filter(FDR < 0.05) # Keep only significant results
}

# --- Run All Comparisons ---
# Define all comparisons to be made
comparisons <- list(
  c("Recipient", "T2", "T3"),
  c("Recipient", "T2", "T4"),
  c("Control", "T2", "T3"),
  c("Control", "T2", "T4")
)

# Apply the function to all comparisons for L2 and L3 data
all_l2_results <- purrr::map_dfr(comparisons, ~run_paired_analysis(l2_tidy, .x[1], .x[2], .x[3]))
all_l3_results <- purrr::map_dfr(comparisons, ~run_paired_analysis(l3_tidy, .x[1], .x[2], .x[3]))

# Save the summary tables
write.csv(all_l2_results, file.path(table_path, "picrust_L2_significant_results.csv"), row.names = FALSE)
write.csv(all_l3_results, file.path(table_path, "picrust_L3_significant_results.csv"), row.names = FALSE)

# --- Plotting Significant Results ---
# Create a plotting function to avoid repetition
plot_significant_kos <- function(tidy_data, sig_results, comparison_str, level_str) {
  # Get the KOs that were significant for this specific comparison
  sig_kos <- sig_results %>%
    filter(Comparison == comparison_str) %>%
    pull(KO_Category)
  
  if (length(sig_kos) == 0) {
    print(paste("No significant KOs for", comparison_str))
    return(NULL)
  }
  
  plot_df <- tidy_data %>%
    filter(KO_Category %in% sig_kos) %>%
    filter(grepl(strsplit(comparison_str, " ")[[1]][1], Treatment)) %>%
    filter(Timepoint %in% c(
      strsplit(comparison_str, " vs ")[[1]][1] %>% strsplit(" ") %>% .[[1]] %>% .[2],
      strsplit(comparison_str, " vs ")[[1]][2]
    ))
  
  p <- ggplot(plot_df, aes(x = Timepoint, y = Abundance, fill = Timepoint)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0) +
    facet_wrap(~KO_Category, scales = "free_y") +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Significantly Changed", level_str, "Pathways"),
      subtitle = comparison_str,
      x = "Timepoint",
      y = "Predicted Relative Abundance"
    ) +
    theme(legend.position = "none")
  
  ggsave(
    file.path(fig_path, paste0("picrust_", level_str, "_", gsub(" ", "_", comparison_str), ".png")),
    p, width = 12, height = 9
  )
}

# Generate plots for all significant comparisons
purrr::walk(unique(all_l2_results$Comparison), ~plot_significant_kos(l2_tidy, all_l2_results, .x, "L2"))
purrr::walk(unique(all_l3_results$Comparison), ~plot_significant_kos(l3_tidy, all_l3_results, .x, "L3"))

print("Functional analysis complete.")