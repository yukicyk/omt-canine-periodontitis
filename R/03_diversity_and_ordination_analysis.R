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
pacman::p_load(phyloseq, ggplot2, vegan, dplyr, ggsignif)

# Define output paths
fig_path <- "results/figures"
table_path <- "results/tables"
if (!dir.exists(fig_path)) dir.create(fig_path)
if (!dir.exists(table_path)) dir.create(table_path)

# Load the phyloseq object
ps <- readRDS("results/phyloseq_object.rds")

# --- 2. Alpha Diversity Analysis ---

# For reproducibility, rarefy to an even depth before alpha diversity calculation.
# The sample.size should be chosen based on the minimum library size of the dataset.
# Check min library size: min(sample_sums(ps))
ps.rf <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123, replace = FALSE)

# Plot alpha diversity (Shannon and Observed ASVs) comparing treatments at each timepoint
alpha_plot <- plot_richness(ps.rf, x = "Treatment", measures = c("Shannon", "Observed")) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.7, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~Timepoint, scales = "free_x") +
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