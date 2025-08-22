#
# SCRIPT: 04_statistical_and_clinical_analysis.R
#
# PURPOSE: Perform advanced statistical analyses, including correlation with
#          clinical variables and predictive modeling with Random Forest.
#
# INPUTS:
#   - `results/phyloseq_object.rds`: The phyloseq object from script 02.
#   - `data/synthetic_clinical_data.csv`: Clinical data (BOP, PPD, Plaque).
#
# OUTPUTS:
#   - Correlation heatmaps saved to `results/figures/`.
#   - Random Forest variable importance plots saved to `results/figures/`.
#   - Statistical tables saved to `results/tables/`.
#

# --- 1. Load Libraries and Data ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(phyloseq, ggplot2, dplyr, tidyr, Hmisc, randomForest, pheatmap)

# Define output paths
fig_path <- "results/figures"
table_path <- "results/tables"

# Load phyloseq object and clinical data
ps <- readRDS("results/phyloseq_object.rds")
clinical_data <- read.csv("data/synthetic_clinical_data.csv")

# --- 2. Data Preparation: Merge Microbiome and Clinical Data ---

# For this analysis, we will agglomerate to the Genus level for stability
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

# Extract abundance table (relative) and metadata
genus_abun <- as.data.frame(otu_table(transform_sample_counts(ps_genus, function(x) x / sum(x))))
meta <- as(sample_data(ps_genus), "data.frame") %>%
  tibble::rownames_to_column("SampleID")

# Clean up clinical data and merge
# Note: The SampleID format must match between the two files.
# Here we assume clinical_data has a 'SampleID' column like 'MT-A-1'.
analysis_df <- meta %>%
  left_join(clinical_data, by = "SampleID") %>%
  # Remove samples with missing clinical data for correlation
  filter(!is.na(PPD) & !is.na(BOP) & !is.na(PlaqueIndex)) %>%
  # Join with abundance data
  left_join(
    t(genus_abun) %>% as.data.frame() %>% tibble::rownames_to_column("SampleID"),
    by = "SampleID"
  )

# --- 3. Correlation Analysis (Genus vs. Clinical) ---

# Select only numeric columns for correlation (clinical vars + genus abundances)
corr_matrix_input <- analysis_df %>%
  select(where(is.numeric)) %>%
  # We only want to correlate genera with clinical variables
  select(PPD, BOP, PlaqueIndex, starts_with("ASV"))

# Run correlation using Hmisc::rcorr (handles missing values and gives p-values)
corr_results <- rcorr(as.matrix(corr_matrix_input), type = "pearson")

# Function to flatten the rcorr object into a tidy table
flatten_rcorr <- function(rcorr_obj) {
  r <- rcorr_obj$r %>% as.data.frame() %>% tibble::rownames_to_column("var1") %>% gather(var2, r, -var1)
  p <- rcorr_obj$P %>% as.data.frame() %>% tibble::rownames_to_column("var1") %>% gather(var2, p, -var1)
  left_join(r, p, by = c("var1", "var2"))
}

corr_table <- flatten_rcorr(corr_results)

# Filter for significant correlations between clinical variables and genera
significant_correlations <- corr_table %>%
  filter(var1 %in% c("PPD", "BOP", "PlaqueIndex") & grepl("ASV", var2)) %>%
  filter(p < 0.05) %>%
  # Optional: filter for stronger correlations
  filter(abs(r) > 0.3) %>%
  # Add taxonomy back for plotting
  left_join(
    as.data.frame(tax_table(ps_genus)) %>% tibble::rownames_to_column("var2"),
    by = "var2"
  )

# Create a heatmap of significant correlations
if (nrow(significant_correlations) > 0) {
  heatmap_data <- significant_correlations %>%
    select(var1, Genus, r) %>%
    pivot_wider(names_from = var1, values_from = r, values_fill = 0) %>%
    tibble::column_to_rownames("Genus")
  
  pheatmap(
    heatmap_data,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    fontsize_number = 8,
    main = "Significant Pearson Correlations (p < 0.05)\nGenus vs. Clinical Parameters",
    filename = file.path(fig_path, "correlation_heatmap.png")
  )
  write.csv(significant_correlations, file.path(table_path, "significant_clinical_correlations.csv"))
}

# --- 4. Predictive Modeling (Random Forest) ---

# Goal: Predict Treatment group (Control vs. Recipient) at T4 based on T2 microbiome
rf_df <- meta %>%
  filter(Timepoint %in% c("T2", "T4"), Treatment %in% c("Control", "Recipient")) %>%
  select(SampleID, Dog, Treatment, Timepoint)

# We need a "baseline" (T2) and "outcome" (T4) for each dog
# Here, we simplify to predict Treatment status at T4
rf_data_t4 <- analysis_df %>%
  filter(Timepoint == "T4", Treatment %in% c("Control", "Recipient"))

# Prepare data for Random Forest
# Remove non-predictor columns and ensure outcome is a factor
rf_final_data <- rf_data_t4 %>%
  select(Treatment, starts_with("ASV")) %>%
  mutate(Treatment = factor(Treatment))

# Run Random Forest
set.seed(123) # for reproducibility
rf_model <- randomForest(
  formula = Treatment ~ .,
  data = rf_final_data,
  ntree = 1000,
  importance = TRUE
)

print(rf_model) # View OOB error rate

# Plot Variable Importance (Top 20 predictors)
imp <- importance(rf_model)
imp_df <- data.frame(MeanDecreaseGini = imp[, "MeanDecreaseGini"]) %>%
  tibble::rownames_to_column("ASV") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  # Add taxonomy
  left_join(
    as.data.frame(tax_table(ps_genus)) %>% tibble::rownames_to_column("ASV"),
    by = "ASV"
  ) %>%
  head(20)

# Create the plot
imp_plot <- ggplot(imp_df, aes(x = reorder(Genus, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 20 Most Important Genera for Predicting Treatment Group at T4",
    x = "Genus",
    y = "Mean Decrease in Gini Impurity"
  )

ggsave(file.path(fig_path, "random_forest_importance.png"), imp_plot, width = 10, height = 8)
write.csv(imp_df, file.path(table_path, "random_forest_importance.csv"))

print("Statistical and clinical analysis complete.")