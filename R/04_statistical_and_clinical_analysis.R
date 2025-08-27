#
# SCRIPT: 04_statistical_and_clinical_analysis.R
#
# PURPOSE: Perform advanced statistical analyses, including Principal Response
#          Curve (PRC), correlation with clinical variables, and predictive
#          modeling with Random Forest.
#
# INPUTS:
#   - `results/phyloseq_object.rds`: The phyloseq object from script 02.
#   - `data/synthetic_clinical_data.csv`: Clinical data (BOP, PPD, Plaque).
#
# OUTPUTS:
#   - PRC plots and tables saved to `results/`.
#   - Correlation heatmaps with FDR correction saved to `results/figures/`.
#   - Random Forest variable importance plots saved to `results/figures/`.
#

# --- 1. Load Libraries and Data ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(phyloseq, ggplot2, dplyr, tidyr, Hmisc, randomForest, pheatmap, vegan)

# Define output paths
fig_path <- "results/figures"
table_path <- "results/tables"
if (!dir.exists(fig_path)) dir.create(fig_path)
if (!dir.exists(table_path)) dir.create(table_path)

# Load phyloseq object and clinical data
ps <- readRDS("results/phyloseq_object.rds")
clinical_data <- read.csv("data/clinical_metadata.csv")

# --- 2. Data Preparation ---

# Create the initial metadata frame from the phyloseq object
meta <- as(sample_data(ps), "data.frame") %>%
  tibble::rownames_to_column("SampleID")

# --- 3. Principal Response Curve (PRC) Analysis ---

print("Performing Principal Response Curve (PRC) analysis...")

# Prepare data for PRC: Use raw counts at ASV level, not genus level.
# Filter for only Control and Recipient groups.
ps_prc <- subset_samples(ps, Treatment %in% c("Control", "Recipient"))

# Extract counts and log-transform as described in the document: log(x+1)
abun_table <- as.data.frame(otu_table(ps_prc))
abun_log <- decostand(abun_table, method = "log")

# Prepare the design matrix with treatment and time variables
design <- as(sample_data(ps_prc), "data.frame") %>%
  select(Treatment, Week) %>%
  mutate(Treatment = as.factor(Treatment),
         Week = as.factor(Week))


# Run the PRC analysis
# The formula specifies that we are looking at the effect of Treatment over Week.
prc_model <- prc(response = abun_log, treatment = design$Treatment, time = design$Week)

# --- Plotting PRC Results ---
# Plot the main PRC diagram (equivalent to Figure S5)
png(file.path(fig_path, "prc_main_plot.png"), width = 10, height = 7, units = "in", res = 300)
plot(prc_model, legpos = NA, main = "Principal Response Curve of OMT", xlab = "Weeks Post-Transplant")
legend("topright", legend = levels(design$Treatment), lty = "solid", col = c(1, 2), bty = "n")
dev.off()

# --- Extracting and Plotting Top ASVs ---
# Get the summary to find species scores
prc_summary <- summary(prc_model)

# Get the top 20 ASVs based on their PRC score (absolute value)
top_asvs <- prc_summary$sp %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  rename(prc_score = 2) %>%
  arrange(desc(abs(prc_score))) %>%
  head(20)

# Create a phyloseq object containing only Donor samples
ps_donor <- subset_samples(ps, Treatment == "Donor")

# Get the number of donor samples from this new object
num_donors <- nsamples(ps_donor)

# Get average abundance of these ASVs in the Donor samples
donor_abun <- ps_donor %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  prune_taxa(top_asvs$ASV, .) %>%
  taxa_sums() / num_donors # Use the pre-calculated variable

top_asvs <- top_asvs %>%
  left_join(
    data.frame(ASV = names(donor_abun), donor_avg_abun = donor_abun),
    by = "ASV"
  )

# Add taxonomy to the table (replicating Table S3)
tax_table_df <- as.data.frame(tax_table(ps)) %>%
  tibble::rownames_to_column("ASV")
prc_results_table <- top_asvs %>%
  left_join(tax_table_df, by = "ASV")

write.csv(prc_results_table, file.path(table_path, "prc_top20_asvs_summary.csv"), row.names = FALSE)

# Create abundance plot for the top 20 PRC scoring ASVs (equivalent to Figure S6)
ps_top_asvs <- prune_taxa(top_asvs$ASV, ps_prc) %>%
  transform_sample_counts(function(x) x / sum(x)) # Convert to relative abundance

plot_df <- psmelt(ps_top_asvs) %>%
  mutate(TaxonLabel = if_else(
    is.na(Genus),
    as.character(OTU), # If Genus is NA, just use the OTU/ASV name
    paste0(OTU, " (", Genus, ")") # Otherwise, create the full "ASV (Genus)" label
  ))

prc_abundance_plot <- ggplot(plot_df, aes(x = Week, y = Abundance, group = Dog, color = Treatment)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.7) +
  facet_wrap(~TaxonLabel, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 10) +
  labs(
    title = "Abundance of Top 20 PRC-Scoring ASVs Over Time",
    x = "Weeks Post-Transplant",
    y = "Relative Abundance"
  )

ggsave(file.path(fig_path, "prc_top20_asv_abundances.png"), prc_abundance_plot, width = 12, height = 10)

# --- 4. Genus-Level Data Preparation ---

print("Preparing data for Genus-level analysis...")

# Step 4.1: Agglomerate taxa to the Genus level using tax_glom.
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
#QC CHECK - Evaluate the effect of removing unclassified taxa.
print("--- QC Check: Evaluating abundance loss from tax_glom(NArm=TRUE) ---")
# Calculate sums before and after agglomeration
sums_before <- sample_sums(ps)
sums_after <- sample_sums(ps_genus)
# Combine into a data frame for comparison
qc_df <- data.frame(
  SampleID = names(sums_before),
  Abundance_Before = sums_before,
  Abundance_After = sums_after[names(sums_before)] # Ensure correct sample matching
) %>%
  mutate(Percent_Retained = (Abundance_After / Abundance_Before) * 100)
# Calculate average retention
average_retention <- mean(qc_df$Percent_Retained)
# Print results and recommendations
print("Percentage of abundance retained per sample after removing taxa unclassified at Genus level:")
print(qc_df %>% select(SampleID, Percent_Retained) %>% arrange(Percent_Retained))
print(sprintf("Average abundance retained across all samples: %.2f%%.", average_retention))
if (average_retention > 90) {
  print("QC PASSED: Average retention is high (>90%). It is safe to proceed.")
} else {
  print("QC WARNING: Average retention is low (<90%). A significant portion of abundance is from taxa unclassified at the Genus level.")
  print("RECOMMENDATION: Consider re-running with tax_glom(..., NArm = FALSE) to keep these taxa, or revisit the taxonomic classification pipeline to improve assignments.")
}
print("--------------------------------------------------------------------")
# Step 4.2: Extract the Genus-level abundance table(samples x taxa).
genus_abun <- as.data.frame(otu_table(transform_sample_counts(ps_genus, function(x) x / sum(x))))

# Step 4.3: Create a mapping to get clean, unique Genus names for column headers.
tax_key <- as.data.frame(tax_table(ps_genus)) %>%
  tibble::rownames_to_column("ASV_ID") %>%
  mutate(Genus_Name = ifelse(is.na(Genus) | Genus == "", paste0("Unassigned_Genus_", ASV_ID), as.character(Genus))) %>%
  mutate(Genus_Name_Unique = make.names(Genus_Name, unique = TRUE))

# Step 4.4: Replace the ASV_ID column names with the clean Genus names.
current_colnames <- colnames(genus_abun)
new_colnames <- tax_key$Genus_Name_Unique[match(current_colnames, tax_key$ASV_ID)]
colnames(genus_abun) <- new_colnames

# Step 4.5: Move the sample names from rownames to a column for joining.
genus_abun_wide <- genus_abun %>%
  tibble::rownames_to_column("SampleID")

# Step 4.6: Create the final, clean Genus-level analysis data frame.
analysis_df_genus <- meta %>%
  inner_join(clinical_data, by = "SampleID") %>%
  inner_join(genus_abun_wide, by = "SampleID")

# Check the resultant dataframe before proceeding
# Check is the rows are samples, columns are genera, content is the relative abundance of genus in a given sample.
head(analysis_df_genus)

# --- 5. Correlation Analysis (Genus Level) ---

print("Performing Genus-level correlation analysis with FDR correction...")

# Select the clinical variables and the new genus columns for correlation.
corr_matrix_input_genus <- analysis_df_genus %>%
  select(PPD_mean, BOP_pc, plaque_pc, all_of(new_colnames))

corr_results_genus <- rcorr(as.matrix(corr_matrix_input_genus), type = "pearson")

# Function to flatten the rcorr object
flatten_rcorr <- function(rcorr_obj) {
  r <- rcorr_obj$r %>% as.data.frame() %>% tibble::rownames_to_column("var1") %>% gather(var2, r, -var1)
  p <- rcorr_obj$P %>% as.data.frame() %>% tibble::rownames_to_column("var1") %>% gather(var2, p, -var1)
  left_join(r, p, by = c("var1", "var2"))
}

corr_table_genus <- flatten_rcorr(corr_results_genus)


# Filter for significant correlations, APPLYING FDR CORRECTION
significant_correlations_genus <- corr_table_genus %>%
  filter(var1 %in% c("PPD_mean", "BOP_pc", "plaque_pc") & var2 %in% new_colnames) %>%
  # Correct for multiple testing using Benjamini-Hochberg (FDR)
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  filter(p_adj < 0.05) %>%
  # Optional: filter for stronger correlations
  filter(abs(r) > 0.3) 

# Create a heatmap of significant correlations
if (nrow(significant_correlations_genus) > 0) {
  heatmap_data <- significant_correlations_genus %>%
    select(var1, var2, r) %>% # var1:clinical meta, var2: Genus 
    pivot_wider(names_from = var1, values_from = r, values_fill = 0) %>%
    tibble::column_to_rownames("var2") #r-values stored
  
  pheatmap(
    heatmap_data,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    fontsize_number = 8,
    labels_col = c("mean PPD", "BOP %", "Plaque %"),
    angle_col = ("0"),
    main = "Significant Pearson Correlations (FDR < 0.05)\nGenus vs. Clinical Parameters",
    filename = file.path(fig_path, "correlation_heatmap_fdr.png")
  )
  write.csv(significant_correlations, file.path(table_path, "significant_clinical_correlations_fdr.csv"))
} else {
  print("No significant correlations found after FDR correction.")
}


# --- 6. Predictive Modeling (Random Forest at Genus Level) ---
# Goal: Predict Treatment group (Control vs. Recipient) at T4 based on T4 microbiome

print("Performing Random Forest analysis at Genus level...")

rf_data_t4_genus <- analysis_df_genus %>%
  filter(Timepoint == "T4", Treatment %in% c("Control", "Recipient"))

# Prepare data: select outcome and predictors
rf_final_data_genus <- rf_data_t4_genus %>%
  select(Treatment, all_of(new_colnames)) %>%
  mutate(Treatment = factor(Treatment))

# Run Random Forest
if (nrow(rf_final_data_genus) > 1 && n_distinct(rf_final_data_genus$Treatment) > 1) {
  set.seed(123)
  rf_model_genus <- randomForest(
    formula = Treatment ~ .,
    data = rf_final_data_genus,
    ntree = 1000,
    importance = TRUE
  )
  
  print(rf_model_genus)
  
  imp <- importance(rf_model_genus)
  imp_df <- data.frame(MeanDecreaseGini = imp[, "MeanDecreaseGini"]) %>%
    tibble::rownames_to_column("Genus_Name_Unique") %>%
    arrange(desc(MeanDecreaseGini)) %>%
    left_join(tax_key, by = "Genus_Name_Unique") %>%
    head(20)
  
  imp_plot <- ggplot(imp_df, aes(x = reorder(Genus, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    coord_flip() +
    theme_bw(base_size = 14) +
    labs(
      title = "Top 20 Most Important Genera for Predicting Treatment Group at T4",
      x = "Genus",
      y = "Mean Decrease in Gini Impurity"
    )
  
  ggsave(file.path(fig_path, "random_forest_importance_genus.png"), imp_plot, width = 10, height = 8)
  write.csv(imp_df, file.path(table_path, "random_forest_importance_genus.csv"))
  
} else {
  print("Insufficient data for Genus-level Random Forest analysis.")
}

print("Genus-level analysis complete.")

# --- Side-Quest: Reproducing Figure 4 ---
# This section reproduces the Principal Response Curve (PRC) plot from the publication.
## If you have followed through all the previous analysis, skip S4.1 & S4.2, 
## those are already processed in the previous sessions.

# --- S4.1. Data Preparation for PRC ---
print("--- Side-Quest: Reproducing Figure 4 ---")
print("Preparing data for Principal Response Curve (PRC) analysis...")

ps_prc <- subset_samples(ps, Treatment %in% c("Control", "Recipient"))
abun_table <- as.data.frame(otu_table(ps_prc))
abun_log <- decostand(abun_table, method = "log")
design <- as(sample_data(ps_prc), "data.frame") %>%
  select(Treatment, Week) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Recipient")),
    Week = as.factor(Week)
  )

# --- S4.2. Run PRC and Prepare Labels for Plotting ---
prc_model <- prc(response = t(abun_log), treatment = design$Treatment, time = design$Week)
prc_summary <- summary(prc_model)

tax_table_df <- as.data.frame(tax_table(ps)) %>%
  tibble::rownames_to_column("ASV")

# Get the top 20 ASVs and create the descriptive labels.
top_asvs_labeled <- prc_summary$sp %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  rename(prc_score = 2) %>%
  arrange(desc(abs(prc_score))) %>%
  head(20) %>%
  left_join(tax_table_df, by = "ASV") %>%
  # **MODIFICATION**: Implement hierarchical fallback for taxon names.
  mutate(
    TaxonName = case_when(
      !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
      !is.na(Genus)                   ~ paste(Genus, "sp."),
      !is.na(Family)                  ~ paste(Family, "(Family)"),
      !is.na(Order)                   ~ paste(Order, "(Order)"),
      !is.na(Class)                   ~ paste(Class, "(Class)"),
      !is.na(Phylum)                  ~ paste(Phylum, "(Phylum)"),
      !is.na(Kingdom)                 ~ Kingdom,
      TRUE                            ~ ASV # Ultimate fallback
    ),
    FinalLabel = paste0(TaxonName, " (", ASV, ")")
  )

# --- S4.3. Generate the PRC Plot (Figure 4 Reproduction) ---
print("Generating PRC plot (Figure 4) with dynamic Y-axis and manual labels...")

# Extract the PRC line data
prc_line_data <- prc_summary$coefficients
# Manually create the numeric time points
numeric_time_points <- c(-2, 0, 2, 12)

# Dynamically determine the y-axis range
range_line <- range(prc_line_data)
range_species <- range(top_asvs_labeled$prc_score)
dynamic_ylim <- range(c(range_line, range_species))

# Add a 10% buffer for better visualization so points aren't on the edge
padding <- (dynamic_ylim[2] - dynamic_ylim[1]) * 0.1
dynamic_ylim_padded <- dynamic_ylim + c(-padding, padding)

png(file.path(fig_path, "Figure4_reproduction_prc_final.png"), width = 11, height = 8, units = "in", res = 300)
par(mar = c(5, 5, 4, 17) + 0.1)

# Step 1: Create the empty plot canvas using our new dynamic ylim
plot(
  NULL,
  xlim = c(-2,12),
  ylim = dynamic_ylim_padded, # Use the calculated, padded range
  type = "n",
  xlab = "Timepoint (Weeks)",
  ylab = "Effect Size (Difference from Control)",
  main = "Figure 4 Reproduction: Principal Response Curve (PRC)",
  xaxt = "n",
  yaxt = "n"
)

# Step 2: Draw the PRC line, points, and axes
lines(numeric_time_points, prc_line_data, col = "darkcyan", lty = "solid", lwd = 2)
points(numeric_time_points, prc_line_data, pch = 16, col = "darkcyan", cex = 1.5)
axis(side = 1, at = seq(-2, 12, by = 2)) 
axis(side = 2, las = 1)
abline(h = 0, lty = "dashed", col = "grey50", lwd = 1.5)

# Step 3: **NEW** - Algorithm to prevent label overlap
# Sort labels by score to process them from bottom to top
labels_df <- top_asvs_labeled %>% arrange(prc_score)
# Calculate minimum separation needed based on text height in plot coordinates
min_sep <- strheight("A", units = "user", cex = 0.9) * 1.5 
# Initialize a new column for the adjusted y-positions
labels_df$new_y <- labels_df$prc_score
# Loop from the second label upwards
if (nrow(labels_df) > 1) {
  for (i in 2:nrow(labels_df)) {
    # Calculate difference to the label below it
    diff <- labels_df$new_y[i] - labels_df$new_y[i-1]
    # If they are too close, nudge the current label up
    if (diff < min_sep) {
      labels_df$new_y[i] <- labels_df$new_y[i-1] + min_sep
    }
  }
}
# Step 4: - Add pointer lines and spaced labels
# Add the internal ticks on the right-hand axis first.
rug(labels_df$prc_score, side = 4, col = "darkred", ticksize = 0.02, lwd = 1.5)

# Define x-coordinates for the start and end of the pointer lines
x_axis_pos <- par("usr")[2] # The right edge of the plot area
x_label_start <- x_axis_pos + (par("usr")[2] - par("usr")[1]) * 0.02 # Start labels slightly off the axis

# Draw the pointer lines (segments)
segments(
  x0 = x_axis_pos, 
  y0 = labels_df$prc_score, 
  x1 = x_label_start, 
  y1 = labels_df$new_y, 
  col = "darkred",
  xpd = TRUE
)

# Draw the text labels at their new, non-overlapping positions
text(
  x = x_label_start,
  y = labels_df$new_y,
  labels = labels_df$FinalLabel,
  pos = 4, # Position text to the right of the coordinate
  xpd = TRUE, # Allow plotting into the margin
  col = "darkred",
  cex = 0.9
)

# Step 5: **MODIFICATION** - Create the complete legend
legend(
  "bottomleft",
  legend = c("OMT (Test Group)", "Control (y=0)"),
  col = c("darkcyan", "grey50"),
  lty = c("solid", "dashed"),
  lwd = c(2, 1.5),
  bty = "n"
)

dev.off()
par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset margins

# --- S4.4. Create and Save Supporting Data Table ---
print("Saving supporting data table for Figure 4...")

final_prc_table <- top_asvs_labeled %>%
  select(ASV, FinalLabel, prc_score, Kingdom, Phylum, Class, Order, Family, Genus, Species)

write.csv(final_prc_table, file.path(table_path, "Figure4_data_prc_top20_asvs_labeled.csv"), row.names = FALSE)

print(paste("Final Figure 4 reproduction saved to:", file.path(fig_path, "Figure4_reproduction_prc_final.png")))
print(paste("Supporting data table with labels saved to:", file.path(table_path, "Figure4_data_prc_top20_asvs_labeled.csv")))
print("------------------------------------------")
