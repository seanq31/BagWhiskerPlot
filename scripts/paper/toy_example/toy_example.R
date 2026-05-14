# ==============================================================================
# Step 1: Setup
# Next: Load required packages and set global options.
# ==============================================================================
library(MASS)

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
# Step 4: Define experiment settings
# Next: Create the parameter grid to iterate over.
# ==============================================================================
dist_settings <- c("toy_eg")
inlier_settings <- c(FALSE)
outter_settings <- c(TRUE)
redefine_loop_settings <- c(TRUE)
conservative_lambda_settings <- c(1)
center_type_settings <- c("hdepth")

all_settings <- expand.grid(
  dist = dist_settings, normal_inlier = inlier_settings, normal_outter = outter_settings, conservative_lambda = conservative_lambda_settings,
  redefine_loop = redefine_loop_settings, center_type = center_type_settings
)

# ==============================================================================
# Step 5: Run the analysis for each setting
# Next: Loop over settings, generate plots, and save outputs.
# ==============================================================================
for (setting_idx in 1:nrow(all_settings)) {
  dist_type <- all_settings[setting_idx, "dist"]
  normal_inlier <- all_settings[setting_idx, "normal_inlier"]
  normal_outter <- all_settings[setting_idx, "normal_outter"]
  conservative_lambda <- all_settings[setting_idx, "conservative_lambda"]
  redefine_loop <- all_settings[setting_idx, "redefine_loop"]
  center_type <- all_settings[setting_idx, "center_type"]

  fig_title <- paste0(
    dist_type,
    ifelse(center_type == "hdepth", "_center_hdepth", "_center_mcd")
  )
  dat <- switch(as.character(dist_type),
    "toy_eg" = dat_toy_eg
  )
  bgplts <- list()

  # ==============================================================================
  # Step 5.1: Compute bag/whisker plots (BagWhiskerPlot)
  # Next: Run the bag/whisker routine under multiple adjustment types.
  # ==============================================================================
  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "unadjusted", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 1, pch = 1, precision = 3,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq", # debug.plots='all',
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = FALSE, naive_bag = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FWER", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 1, pch = 1, precision = 3,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq", # debug.plots='all',
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, naive_bag = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FDR", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 1, pch = 1, precision = 3,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq", # debug.plots='all',
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, naive_bag = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))


  old <- Sys.time()
  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "PFER", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 1, pch = 1, precision = 3,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq", # debug.plots='all',
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, naive_bag = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))
  print(Sys.time() - old)


  # ==============================================================================
  # Step 5.2: Prepare ggplot objects
  # Next: Extract the ggplot objects and compute shared axis limits.
  # ==============================================================================
  library(gridExtra)
  library(ggplot2)
  plts_of_bgplts <- lapply(bgplts, function(item) item$bgplt)

  # compute union xlim/ylim across all plots
  get_xy_limits <- function(p) {
    gb <- ggplot_build(p)
    xs <- unlist(lapply(gb$data, function(d) d[["x"]]), use.names = FALSE)
    ys <- unlist(lapply(gb$data, function(d) d[["y"]]), use.names = FALSE)
    xs <- xs[is.finite(xs)]
    ys <- ys[is.finite(ys)]
    list(
      x = if (length(xs)) range(xs, na.rm = TRUE) else c(NA_real_, NA_real_),
      y = if (length(ys)) range(ys, na.rm = TRUE) else c(NA_real_, NA_real_)
    )
  }
  lims_list <- lapply(plts_of_bgplts, get_xy_limits)
  xlim_union <- c(-8, 22)
  ylim_union <- c(-5, 25)

  # ==============================================================================
  # Step 5.3: Save reference plot (aplpack)
  # Next: Create an aplpack bagplot and write it to disk.
  # ==============================================================================
  dir.create(file.path(save_dir, paste0("simu_", dist_type)), recursive = TRUE)
  png(file.path(save_dir, paste0("simu_", dist_type), paste0(fig_title, "_aplpack.png")),
    width = 3, height = 3, units = "in", res = 300
  )
  par(mar = c(3.5, 3.5, 3.3, 1.5), mgp = c(1.5, 0.5, 0))
  bo_apl <- aplpack::bagplot(dat,
    factor = 3, type1 = "unadjusted", normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq",
    show.whiskers = FALSE, show.loophull = TRUE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE,
    show.fence_mag_bag = TRUE, xlim = xlim_union, ylim = ylim_union,
  )
  dev.off()

  # ==============================================================================
  # Step 5.4: Save BagWhiskerPlot subplots
  # Next: Apply shared limits and save each subplot as a PNG.
  # ==============================================================================
  # update each subplot with the union limits
  plts_of_bgplts <- lapply(plts_of_bgplts, function(p) {
    p + coord_cartesian(xlim = xlim_union, ylim = ylim_union) +
      theme(plot.margin = margin(t = 3.5, r = 1.3, b = 1.8, l = 1.8, unit = "lines"))
  })

  for (j in seq_along(plts_of_bgplts)) {
    ggsave(
      filename = file.path(
        save_dir,
        paste0("simu_", dist_type, "/", fig_title, "_sub", j, ".png")
      ),
      plot = plts_of_bgplts[[j]],
      width = 4,
      height = 4,
      bg = "white"
    )
  }

  # ==============================================================================
  # Step 5.5: Combine selected images into a single figure
  # Next: Read saved panels, add labels, and export the combined figure.
  # ==============================================================================
  library(cowplot)
  library(png)
  library(grid)

  subdir <- file.path(save_dir, paste0("simu_", dist_type))
  sub_paths <- file.path(subdir, paste0(fig_title, "_sub", 1:2, ".png"))

  # Read images
  read_img_grob <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  }
  grobs <- lapply(sub_paths, read_img_grob)

  # Wrap with labels
  labels <- c("(a) Unadjusted", "(b) Adjusted")
  labeled <- mapply(function(g, lb) {
    ggdraw() +
      draw_grob(g) +
      draw_label(lb,
        x = 0.02, y = 0.98, hjust = 0, vjust = 1,
        fontface = "bold", size = 15, color = "black"
      )
  }, grobs, labels, SIMPLIFY = FALSE)

  combined <- plot_grid(plotlist = labeled, ncol = 2)
  ggsave(
    filename = file.path(save_dir, paste0(fig_title, "_combined.png")),
    plot = combined,
    width = 8, height = 4, dpi = 300, bg = "white"
  )
}
