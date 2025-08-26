#
# SCRIPT: 03_diversity_and_ordination_analysis.R
#
# PURPOSE: Perform core ecological analyses, including alpha and beta diversity.
#          Generates ordination plots (PCoA/MDS) and statistical tests (PERMANOVA).
#
# INPUTS:
#   - `results/phyloseq_object.rds`: The phyloseq object from script 02.
#
# OUTPUTS:
#   - Alpha and beta diversity plots saved to `results/figures/`.
#   - Statistical results tables saved to `results/tables/`.
#

# --- 1. Load Libraries and Data ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(phyloseq, ggplot2, vegan, dplyr, rstatix, ggpubr, devtools)
library(devtools)
devtools::install_github("gauravsk/ranacapa")
library(ranacapa)

# Define output paths
fig_path <- "results/figures"
table_path <- "results/tables"
if (!dir.exists(fig_path)) dir.create(fig_path)
if (!dir.exists(table_path)) dir.create(table_path)

# Load the phyloseq object
ps <- readRDS("results/phyloseq_object.rds")

# --- 2. Alpha Diversity Analysis ---

# Plot rarefaction curves for all samples
# Accessing sampling depth and coverage
p.rare <- ggrare(ps, step = 1000, color = "Treatment",label="Sample", se = FALSE)
p.rare + xlim(0, 25000)
p.all <- p.rare + facet_wrap(~Treatment)+ xlim(0, 25000)
ggsave(file.path(fig_path, "Rarefaction curves.png"), p.all, width = 10, height = 8)

# For reproducibility, rarefy to an even depth before alpha diversity calculation.
# The sample.size should be chosen based on the minimum library size of the dataset.
# Rarefied dataset is only used for alpha-diversity estimation.
# Check min library size
print(min(sample_sums(ps)))
ps.rf <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123, replace = FALSE)

# Plot alpha diversity (Shannon and Observed ASVs) comparing treatments at each timepoint
alpha_plot <- plot_richness(ps.rf, x = "Timepoint", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.7, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Alpha Diversity by Treatment at Each Timepoint")

ggsave(file.path(fig_path, "alpha_diversity_by_timepoint.png"), alpha_plot, width = 10, height = 8)

# Plot alpha diversity change over time, faceted by treatment
alpha_time_plot <- plot_richness(ps.rf, x = "Week", measures = "Shannon", color = "Dog") +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(aes(group = Dog), alpha = 0.5) +
  facet_wrap(~Treatment) +
  theme_bw(base_size = 14) +
  labs(title = "Shannon Diversity Trajectory Over Time", x = "Week", y = "Shannon Index")

ggsave(file.path(fig_path, "alpha_diversity_trajectory.png"), alpha_time_plot, width = 12, height = 6)


# ---Side-Quest: Reproducing Figure 2 ---
# --- S1.1. Data Preparation ---

# Calculate alpha diversity metrics (Shannon)
alpha_div <- estimate_richness(ps, measures = "Shannon")

# Combine diversity data with sample metadata
# The sample_data should contain 'Timepoint' and 'Treatment' columns
plot_data <- data.frame(alpha_div,sample_data(ps))

# Re-label the Timepoint and Treatment to match  Figure 2
plot_data$Timepoint <- factor(plot_data$Timepoint, 
                              levels = c("T1", "T2", "T3", "T4"),
                              labels = c("week –2", "baseline", "week 2", "week 12"))
plot_data$Treatment <- factor(plot_data$Treatment, 
                              levels = c( "Donor", "Recipient","Control"), 
                              labels = c( "Donor", "Test", "Control"))

# Define the colors used in the paper
treatment_colors <- c("Control" = "cadetblue1", "Test" = "chartreuse2", "Donor" = "coral1")

# --- S1.2. Create the Upper Panel Plot (Between-Group Comparisons) ---
# This performs pairwise Wilcoxon tests for each timepoint.
stat.test.upper <- plot_data %>%
  group_by(Timepoint) %>%
  wilcox_test(Shannon ~ Treatment, p.adjust.method = "fdr") %>%
  add_significance("p.adj") %>%
  filter(
    (group1 == "Donor" & group2 %in% c("Test", "Control")) |
      (group1 == "Test" & group2 == "Control")
  )
# This automatically calculates the correct y-coordinates for the brackets, avoiding errors.
stat.test.upper <- stat.test.upper %>%
  add_y_position(fun = "max", step.increase = 0.12) %>%
  add_x_position(x = "Timepoint", dodge = 0.8)

#Build the plot and add the p-values manually
upper_panel <- ggplot(plot_data, aes(x = Timepoint, y = Shannon, fill = Treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, alpha = 0.7) +
  
  # This function will now find the 'p.signif' column it needs in stat.test.upper
    stat_pvalue_manual(
    stat.test.upper, 
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = FALSE,
    inherit.aes = FALSE 
  ) +
  scale_fill_manual(values = treatment_colors) +
  # Adjust y-axis to ensure there is enough space for the automatically placed brackets
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "Comparison of Alpha Diversity Between Groups", x = "", y = "Shannon Entropy") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))


# --- S1.3. Create the Lower Panel Plot (Within-Group Comparisons) ---
# Pre-calculate the statistics for the lower panel
stat.test.lower <- plot_data %>%
  group_by(Treatment) %>%
  wilcox_test(Shannon ~ Timepoint, paired = FALSE, p.adjust.method = "fdr") %>%
  add_significance("p.adj")
# Add the xy positions for plotting. This now runs on the correctly calculated stats.
stat.test.lower <- stat.test.lower %>%
  add_xy_position(x = "Timepoint", fun = "max", step.increase = 0.12)

lower_panel <- ggplot(plot_data, aes(x = Timepoint, y = Shannon, fill = Treatment)) +
  geom_boxplot(color = "black", alpha = 0.8, show.legend = FALSE, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75), alpha = 0.4, show.legend = FALSE) +
  
  stat_pvalue_manual(
    stat.test.lower,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = FALSE, # Set to FALSE to see "ns" for non-significant results
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = treatment_colors) +
  facet_wrap(~Treatment, scales = "free_x") + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "Changes in Alpha Diversity Within Groups Over Time", x = "Timepoint", y = "Shannon Entropy") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# --- S1.4. Combine Plots and Save ---

# Arrange the two panels into a single figure
final_plot <- ggarrange(
  upper_panel, 
  lower_panel, 
  ncol = 1, 
  nrow = 2,
  common.legend = TRUE, # Use a single, shared legend
  legend = "right"
)

# Add a main title to the combined plot
final_plot_with_title <- annotate_figure(
  final_plot,
  top = text_grob("Figure 2 Reproduction: Changes in Oral Microbiota Alpha-Diversity", 
                  face = "bold", size = 16)
)

# Save the final figure 
# Figure 2 Changes in the alpha-diversities of the oral microbiota in the donor, test, and control groups throughout the study
ggsave(
  file.path(fig_path, "Figure2_reproduction_alpha_diversity.png"), 
  final_plot_with_title, 
  width = 10, 
  height = 12
)

print(final_plot_with_title)
print(paste("Figure saved to:", file.path(fig_path, "Figure2_reproduction_alpha_diversity.png")))

# --- 3. Beta Diversity Analysis ---

# Transform to relative abundance for Bray-Curtis and other ecological metrics
ps.ra <- transform_sample_counts(ps, function(x) { x / sum(x) })

# --- Ordination ---
# Calculate distances
dist_bc <- distance(ps.ra, method = "bray")
dist_wuf <- distance(ps, method = "wunifrac") # Use non-rarefied, non-relative for UniFrac

# Perform MDS/PCoA
ord_bc <- ordinate(ps.ra, "MDS", dist_bc)
ord_wuf <- ordinate(ps, "PCoA", dist_wuf)

# Plot Bray-Curtis ordination
pcoa_bc_plot <- plot_ordination(ps.ra, ord_bc, color = "Treatment", shape = "Timepoint") +
  geom_point(size = 5, alpha = 0.75) +
  stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity",
    color = "Treatment Group",
    shape = "Timepoint"
  )

ggsave(file.path(fig_path, "pcoa_bray_curtis.png"), pcoa_bc_plot, width = 8, height = 6)

# Plot Weighted UniFrac ordination
pcoa_wuf_plot <- plot_ordination(ps, ord_wuf, color = "Treatment", shape = "Timepoint") +
  geom_point(size = 5, alpha = 0.75) +
  stat_ellipse(aes(group = Treatment), type = "t", linetype = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCoA of Weighted UniFrac Distance",
    color = "Treatment Group",
    shape = "Timepoint"
  )

ggsave(file.path(fig_path, "pcoa_weighted_unifrac.png"), pcoa_wuf_plot, width = 8, height = 6)

# --- Statistical Testing (PERMANOVA) ---
metadata <- as(sample_data(ps), "data.frame")

# Test overall effect of Treatment and Timepoint
adonis_result <- adonis2(dist_bc ~ Treatment * Timepoint, data = metadata, permutations = 9999)
print(adonis_result)

# Save the result to a file
write.csv(as.data.frame(adonis_result), file.path(table_path, "permanova_bray_curtis_results.csv"))

print("Diversity and ordination analysis complete.")

# --- Side-Quest: Reproducing Figure 3 ---
# This section focuses on creating a publication-quality PCoA plot,
# faceted by timepoint, with appropriate PERMANOVA statistics displayed on each facet.

# --- S2.1. Data Preparation and Transformation ---
# We will use the relative abundance phyloseq object 'ps.ra' already created in the main script.
# For consistency, we ensure the factor levels and colors from the Figure 2 reproduction are used.
sample_data(ps.ra)$Timepoint <- factor(sample_data(ps.ra)$Timepoint, 
                                       levels = c("T1", "T2", "T3", "T4"),
                                       labels = c("week –2", "baseline", "week 2", "week 12"))
sample_data(ps.ra)$Treatment <- factor(sample_data(ps.ra)$Treatment, 
                                       levels = c("Donor", "Recipient", "Control"), 
                                       labels = c("Donor", "Test", "Control"))
treatment_colors <- c("Control" = "cadetblue1", "Test" = "chartreuse2", "Donor" = "coral1")

# --- S2.2. Beta Diversity Calculation (PCoA) ---
# Perform PCoA on the Bray-Curtis dissimilarity matrix calculated from relative abundance data.
# The 'ordinate' function is a wrapper that handles distance calculation and ordination.
pcoa_bc <- ordinate(ps.ra, method = "PCoA", distance = "bray")

# --- S2.3. PERMANOVA Statistical Analysis per Timepoint ---
# To get stats for each facet, we must run PERMANOVA on each timepoint subset.
# A single adonis2 call on the whole dataset would not give per-timepoint results.

# Get the unique timepoints from the data
timepoints <- levels(sample_data(ps.ra)$Timepoint)
permanova_results <- list()

for (tp in timepoints) {
  # Subset the relative abundance phyloseq object to the current timepoint
  ps_subset <- subset_samples(ps.ra, Timepoint == tp)
  
  # Check if the subset has enough data (e.g., multiple treatment groups) to test
  if (nsamples(ps_subset) > 1 && length(unique(sample_data(ps_subset)$Treatment)) > 1) {
    
    # Get sample data from the subset
    sample_df_subset <- data.frame(sample_data(ps_subset))
    
    # Calculate Bray-Curtis distance matrix for only the subset samples
    dist_matrix <- phyloseq::distance(ps_subset, method = "bray")
    
    # Run PERMANOVA using the vegan::adonis2 function
    permanova_fit <- adonis2(dist_matrix ~ Treatment, data = sample_df_subset, permutations = 9999)
    
    # Store the relevant results (R2 and p-value for the 'Treatment' term)
    permanova_results[[tp]] <- data.frame(
      Timepoint = tp,
      R2 = permanova_fit$R2[1],
      p_value = permanova_fit$`Pr(>F)`[1]
    )
  }
}

# Combine the list of results into a single, tidy data frame
permanova_stats_df <- do.call(rbind, permanova_results)

# Format the statistical results for clean plotting
permanova_stats_df <- permanova_stats_df %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"), # Adjust p-values for multiple comparisons (4 tests)
    p_label = if_else(p_adj < 0.001, "p < 0.001", paste("p =", round(p_adj, 3))),
    r2_label = paste("R² =", round(R2, 2)),
    full_label = paste(r2_label, p_label, sep = "\n") # Combine into one label for the plot
  )

# --- S2.4. Create the PCoA Plot ---
# Extract the variance explained by the first two axes to use in axis labels
axis1_var <- round(pcoa_bc$values$Relative_eig[1] * 100, 1)
axis2_var <- round(pcoa_bc$values$Relative_eig[2] * 100, 1)

# Build the plot layer by layer using ggplot2 for full control
fig3_plot <- plot_ordination(
  physeq = ps.ra,
  ordination = pcoa_bc,
  color = "Treatment"
) +
  geom_point(size = 4, alpha = 0.75) +
  # Facet by timepoint to create separate panels for each timepoint
  facet_wrap(~Timepoint, ncol=2, as.table =TRUE, dir="h") +
  # Add 95% confidence ellipses for each treatment group
  stat_ellipse(type = "t", level = 0.95, linetype = 2, aes(group = Treatment)) +
  # Apply our consistent color scheme
  scale_color_manual(values = treatment_colors) +
  # Add the PERMANOVA stats as text to each facet
  # We place the text at the top-left corner of each facet panel
  geom_text(
    data = permanova_stats_df, 
    aes(x = -Inf, y = Inf, label = full_label), 
    hjust = -0.1, vjust = 1.5, inherit.aes = FALSE, size = 4
  ) +
  # Set labels and theme for a publication-quality graphic
  labs(
    title = "Figure 3 Reproduction: Beta Diversity (PCoA) by Treatment Over Time",
    subtitle = "Bray-Curtis Dissimilarity with PERMANOVA Statistics per Timepoint",
    x = paste0("Axis 1 (", axis1_var, "%)"),
    y = paste0("Axis 2 (", axis2_var, "%)"),
    color = "Treatment Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90")
  )

# --- S2.5. Save and Display the Plot ---
ggsave(
  file.path(fig_path, "Figure3_reproduction_beta_diversity.png"), 
  fig3_plot, 
  width = 10, 
  height = 8
)

print(fig3_plot)
print(paste("Figure saved to:", file.path(fig_path, "Figure3_reproduction_beta_diversity.png")))
