# ==============================================================================
# Step 1: Setup
# Next: Load required packages and set global options.
# ==============================================================================
library(MASS)
library(ggplot2)
library(ggpattern)

set.seed(1234)

# ==============================================================================
# Step 2: Configure output paths
# Next: Define where figures will be written.
# ==============================================================================
save_dir <- "./figures/"

# ==============================================================================
# Step 3: Define toy example data
# Next: Create the fixed 2D dataset used in the toy example.
# ==============================================================================
dat_toy_eg <- matrix(
  c(
    7, 5,
    7, 7,
    9, 4,
    5, 4,
    14, 9,
    0, 9,
    7, -3,
    19, 20
  ),
  ncol = 2, byrow = TRUE
)

colnames(dat_toy_eg) <- c("x", "y")

# ==============================================================================
# Step 4: Prepare plotting data
# Next: Convert to a data frame and create labels/styles for each point.
# ==============================================================================
# Convert matrix to data frame for ggplot
df_toy <- as.data.frame(dat_toy_eg)

# Add index and label expression for each point (for parse = TRUE)
df_toy$idx <- seq_len(nrow(df_toy))
df_toy$label <- paste0(
  "z[", df_toy$idx, ']*"=(',
  df_toy$x, ", ", df_toy$y, ')"'
)

blue_points <- df_toy[2:7, c("x", "y")]
blue_points_list <- split(blue_points, seq(nrow(blue_points)))
for (pt in blue_points_list) {
  df_toy$color[df_toy$x == pt$x & df_toy$y == pt$y] <- "blue"
}

df_toy$color[df_toy$x == df_toy$x[1] & df_toy$y == df_toy$y[1]] <- "red"

df_toy$color[df_toy$x == df_toy$x[nrow(df_toy)] & df_toy$y == df_toy$y[nrow(df_toy)]] <- "darkred"

df_toy$ptype <- c(
  "star",
  "hollow", "hollow", "hollow",
  "solid", "solid", "solid",
  "large"
)

# ==============================================================================
# Step 5: Configure plot geometry
# Next: Set axis limits so half-spaces and annotations are visible.
# ==============================================================================
# Define plot limits to ensure half-spaces are visible within a reasonable range
x_limits <- c(-4, 20)
y_limits <- c(-4, 20)

# ==============================================================================
# Step 6: Build the figure
# Next: Draw half-spaces, bag polygon, points, and text annotations.
# ==============================================================================
# Create the plot
p <- ggplot(df_toy, aes(x = x, y = y, color = color, shape = ptype)) +
  scale_color_identity(guide = "none") +

  # PURPLE half-space: y <= -0.55x+8.85
  geom_abline(intercept = 8.85, slope = -0.55, color = "red", linetype = "dashed", alpha = 0.3) +
  geom_polygon_pattern(
    data = data.frame(
      x = c(-10, 25, 25, -10),
      y = c(-0.55 * -10 + 8.85, -0.55 * 25 + 8.85, -16, -16)
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    fill = NA,
    color = NA,
    pattern = "stripe",
    pattern_fill = "red",
    pattern_colour = "red",
    pattern_angle = 135,
    pattern_spacing = 0.02,
    pattern_density = 0.05,
    pattern_size = 0.006,
    pattern_alpha = 0.3
  ) +

  # BLUE half-space: y <= 2x - 19
  geom_abline(intercept = -19, slope = 2, color = "blue", linetype = "dashed", alpha = 0.3) +
  geom_polygon_pattern(
    data = data.frame(
      x = c(-10, 25, 25, -10),
      y = c(2 * -10 - 19, 2 * 25 - 19, -16, -16)
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    fill = NA,
    color = NA,
    pattern = "stripe",
    pattern_fill = "blue",
    pattern_colour = "blue",
    pattern_angle = 0,
    pattern_spacing = 0.02,
    pattern_density = 0.05,
    pattern_size = 0.006,
    pattern_alpha = 0.3
  ) +

  # ORANGE half-space: y >= 1.2x - 1.4
  geom_abline(intercept = -1.4, slope = 1.2, color = "orange", linetype = "dashed", alpha = 0.3) +
  geom_polygon_pattern(
    data = data.frame(
      x = c(-10, 20, 20, -10),
      y = c(1.2 * -10 - 1.4, 1.2 * 20 - 1.4, 36, 36)
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    fill = NA,
    color = NA,
    pattern = "stripe",
    pattern_fill = "orange",
    pattern_colour = "orange",
    pattern_angle = 45,
    pattern_spacing = 0.02,
    pattern_density = 0.05,
    pattern_size = 0.006,
    pattern_alpha = 0.3
  ) +

  # polygon bag
  geom_polygon(
    data = data.frame(
      x = c(
        dat_toy_eg[2, 1], dat_toy_eg[3, 1], dat_toy_eg[4, 1]
      ),
      y = c(
        dat_toy_eg[2, 2], dat_toy_eg[3, 2], dat_toy_eg[4, 2]
      )
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +


  # Points last so they appear on top of polygons
  geom_point(size = 3, stroke = 1.2) +

  # Add LaTeX/markdown-style labels slightly below each point
  geom_text(
    aes(label = label),
    parse = TRUE,
    vjust = 1.8,
    nudge_x = c(-6.6, -0.65, -0.05, -1.85, -0.95, -0.85, -0.65, -0.9),
    nudge_y = c(2.2, 2.5, 0, 0, 0.2, 0.2, 2.4, 0),
    size = 6,
    color = "black",
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c("solid" = 1, "hollow" = 1, "star" = 8, "large" = 19),
    guide = "none"
  ) +
  scale_color_identity(guide = "none") +

  # Set limits and theme
  scale_x_continuous(breaks = seq(floor(x_limits[1]), ceiling(x_limits[2]), by = 3)) +
  scale_y_continuous(breaks = seq(floor(y_limits[1]), ceiling(y_limits[2]), by = 4)) +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dotted"),
    panel.grid.minor = element_line(color = "lightgray", linetype = "dotted"),
    axis.text = ggplot2::element_text(size = 18),
    axis.title = ggplot2::element_text(size = 18),
    plot.title = ggplot2::element_text(size = 18),
    legend.text = ggplot2::element_text(size = 18),
    legend.title = ggplot2::element_text(size = 18),
    strip.text = ggplot2::element_text(size = 18)
  )

# ==============================================================================
# Step 7: Add annotations
# Next: Add an arrow annotation and finalize background styling.
# ==============================================================================
# add arrow from (5,5) to (12,6)
p <- p +
  geom_segment(
    data = data.frame(x = 7, y = 5, xend = 3, yend = 5.8),
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth = 0.7,
    color = "black"
  )

# white background
p <- p + theme(
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white"),
  panel.border = element_blank()
)

# ==============================================================================
# Step 8: Export outputs
# Next: Save the figure and display it in the current session.
# ==============================================================================
ggsave(filename = paste0(save_dir, "toy_example_hdepth.png"), plot = p, width = 6, height = 6)

print(p)
