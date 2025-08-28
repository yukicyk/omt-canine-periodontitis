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
      Week == 0  ~ -2,
      Week == 2 ~ 0,
      Week == 4   ~ 2,
      Week == 14  ~ 12,
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
  renderer = magick_renderer())

print(paste("Animation saved to:", file.path(fig_path, "community_trajectory_animated.gif")))

### Just one more figure ###
# Include the donor in the same ordination space and see how the control and test samples move around
# Create a single data frame with ordination scores and metadata
# This makes plotting with ggplot2 more explicit and flexible.
mds_all <- ordinate(ps.ra, "MDS", "bray")
plot_data <- plot_ordination(ps.ra, mds_all , justDF = TRUE) %>%
  bind_cols(as(sample_data(ps.ra), "data.frame"))

# Separate data for the animated points (Control/Recipient) and static points (Donor)
animated_data <- plot_data %>% filter(Treatment...5 %in% c("Control", "Recipient"))
static_data <- plot_data %>% filter(Treatment...5 == "Donor")

# --- Create Static Plot (Base for Animation) ---
# Create the base ggplot object without faceting.
# We will use 'shape' to distinguish between Treatment groups.
p_unified <- ggplot(animated_data, aes(x = Axis.1, y = Axis.2, color = Dog...3)) +
  # Add the static Donor points as a fixed background layer.
  # They are larger and have a different shape to stand out.
  geom_point(
    data = static_data,
    aes(shape = Treatment...5), # Map shape to Treatment
    size = 8,
    alpha = 0.8,
    color = "royalblue4" # Make donors a distinct, constant color
  ) +
  # Add the points for Control and Recipient that will be animated.
  geom_point(aes(shape = Treatment...5), size = 6, alpha = 0.7) +
  theme_bw(base_size = 16) +
  labs(
    x = "MDS Axis 1",
    y = "MDS Axis 2",
    title = 'Microbiome Trajectory (Weeks: {frame_along})',
    subtitle = 'Bray-Curtis Dissimilarity | All Groups in One Space',
    color = "Dog ID",
    shape = "Group"
  ) +
  # Manually define shapes for clarity
  scale_shape_manual(values = c("Control" = 16, "Recipient" = 17, "Donor" = 18)) +
  # Add text labels for each dog
  geom_text(aes(label = Dog...3), vjust = -1.5, size = 4, show.legend = FALSE)

# --- 5. Add Animation Layers ---

# Add the gganimate layers to the plot
anim_unified <- p_unified +
  transition_reveal(WeekNumeric...12) +
  shadow_wake(wake_length = 0.1, alpha = 0.5) +
  ease_aes('cubic-in-out')

# --- 6. Render and Save Animation ---

print("Rendering unified animation for all groups...")
animate(
  anim_unified,
  nframes = 150,
  fps = 10,
  duration = 15,
  width = 1200,
  height = 800, # Increased height slightly for the single plot
  renderer = magick_renderer(),
  device = "ragg_png" # Use a high-quality rendering device
)

# Save the animation
anim_save(
  filename = file.path(fig_path, "community_trajectory_unified_animated.gif")
)

print(paste("Unified animation saved to:", file.path(fig_path, "community_trajectory_unified_animated.gif")))
