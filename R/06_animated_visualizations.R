#
# SCRIPT: 06_animated_visualizations.R
#
# PURPOSE: Create animated ordination plots to visualize microbial community
#          shifts over time in the canine OMT study using gganimate.
#
# INPUTS:
#   - `results/phyloseq_object.rds`: The phyloseq object from script 02.
#
# OUTPUTS:
#   - Animated GIF files saved to `results/figures/`.
#
# DEPENDENCIES: gganimate, magick, phyloseq, dplyr, ggplot2
#

# --- 1. Load Libraries and Data ---

if (!require("pacman")) install.packages("pacman")
pacman::p_load(phyloseq, ggplot2, dplyr, gganimate, magick)

# Define output path
fig_path <- "results/figures"
if (!dir.exists(fig_path)) dir.create(fig_path)

# Load the phyloseq object
ps <- readRDS("results/phyloseq_object.rds")

# --- 2. Data Preparation ---
# The original 'Week' variable is a factor with non-numeric levels ("baseline").
# gganimate's transition_reveal() works most reliably with a numeric variable.
# We will create a new, clean numeric column for time.

# Extract sample data
sample_df <- as(sample_data(ps), "data.frame")

# Create the new numeric week column
sample_df <- sample_df %>%
  mutate(
    WeekNumeric = case_when(
      Week == "week -2"  ~ -2,
      Week == "baseline" ~ 0,
      Week == "week 2"   ~ 2,
      Week == "week 12"  ~ 12,
      TRUE               ~ NA_real_ # Fallback for any unexpected values
    )
  )
# Replace the old sample data in the phyloseq object with our modified version
sample_data(ps) <- sample_df

# Transform counts to relative abundance for beta diversity analysis
ps.ra <- transform_sample_counts(ps, function(x) { x / sum(x) })

# Create a combined subset for Control and Recipient groups
# This allows for a single ordination space for better comparison
rc.ra <- subset_samples(ps.ra, Treatment %in% c("Control", "Recipient")) %>%
  prune_taxa(taxa_sums(.) > 0, .)

# --- 3. Ordination ---

# Calculate Bray-Curtis dissimilarity and perform MDS ordination on the combined data
mds_rc <- ordinate(rc.ra, "MDS", "bray")

# --- 4. Create Static Plot (Base for Animation) ---

# Create the base ggplot object, faceted by Treatment
p_rc <- plot_ordination(rc.ra, mds_rc, "samples", color = "Dog") +
  geom_point(size = 6, alpha = 0.7) +
  theme_bw(base_size = 16) +
  facet_wrap(~Treatment) +
  labs(
    x = "MDS Axis 1",
    y = "MDS Axis 2",
    title = 'Microbiome Trajectory (Weeks: {frame_along})',
    subtitle = 'Bray-Curtis Dissimilarity',
    color = "Dog ID"
  ) +
  # Add text labels for each dog for clarity
  geom_text(aes(label = Dog), vjust = -1.5, size = 4, show.legend = FALSE)

# --- 5. Add Animation Layers ---

# Add the gganimate layers to the plot
# `transition_reveal(Week)` animates the plot along the 'Week' variable.
# `shadow_wake` leaves a trail of past points.
anim_rc <- p_rc +
  transition_reveal(WeekNumeric) +
  shadow_wake(wake_length = 0.1, alpha = 0.5) +
  ease_aes('cubic-in-out') # Smooths the transition

# --- 6. Render and Save Animation ---

# Render the animation as a high-quality GIF file.
# magick_renderer is used instead of the default gifski
# The duration, fps, and resolution can be adjusted.
print("Rendering animation for Control and Recipient groups...")
animate(
  anim_rc,
  nframes = 150, # More frames for a smoother animation
  fps = 10,
  duration = 15,
  width = 1200,
  height = 700,
  renderer = magick_renderer(file.path(fig_path, "community_trajectory_animated.gif"))
)

print(paste("Animation saved to:", file.path(fig_path, "community_trajectory_animated.gif")))
