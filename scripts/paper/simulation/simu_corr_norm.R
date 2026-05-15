# install.packages("aplpack")

# ==============================================================================
# Step 1: Setup
# Next: Load required packages and set global options.
# ==============================================================================
library(MASS)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(png)
library(grid)
set.seed(1234)

# ==============================================================================
# Step 2: Configure output paths
# Next: Define where figures will be written.
# ==============================================================================
save_dir <- "./figures/"

# ==============================================================================
# Step 3: Simulate correlated mixture data
# Next: Construct a correlated two-component normal mixture dataset.
# ==============================================================================
corr_coef <- 0.3

p <- 0.05
n <- 5000
theta <- rbinom(n, 1, p)

mu0 <- c(100, 300)
mu1 <- c(103, 298)
dx <- mu1[1] - mu0[1]
dy <- mu1[2] - mu0[2]
between_cov <- p * (1 - p) * dx * dy
var_x <- 1 + p * (1 - p) * dx^2
var_y <- 1 + p * (1 - p) * dy^2
within_cov <- corr_coef * sqrt(var_x * var_y) - between_cov
within_cov <- max(min(within_cov, 0.999999), -0.999999)
Sigma <- matrix(c(1, within_cov, within_cov, 1), nrow = 2)
dat_normal <- matrix(NA_real_, nrow = n, ncol = 2)
idx0 <- which(theta == 0)
idx1 <- which(theta == 1)
if (length(idx0)) dat_normal[idx0, ] <- MASS::mvrnorm(length(idx0), mu = mu0, Sigma = Sigma)
if (length(idx1)) dat_normal[idx1, ] <- MASS::mvrnorm(length(idx1), mu = mu1, Sigma = Sigma)
colnames(dat_normal) <- c("x", "y")

# ==============================================================================
# Step 4: Define experiment settings
# Next: Create the parameter grid to iterate over.
# ==============================================================================
dist_settings <- c("corr_normal_mixture")
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
    "corr_normal_mixture" = dat_normal,
    "chisq" = dat_chisq
  )
  bgplts <- list()

  # ==============================================================================
  # Step 5.1: Compute bag/whisker plots (BagWhiskerPlot)
  # Next: Run the bag/whisker routine under multiple adjustment types.
  # ==============================================================================
  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FWER", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FDR", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "PFER", q = 0.5, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "chisq",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FWER", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "F",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "FDR", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "F",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))

  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    factor = 3, type1 = "PFER", q = 0.5, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "F",
    show.whiskers = TRUE, show.loophull = FALSE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = NULL, timing = TRUE
  )
  bgplts <- c(bgplts, list(tmp_plt))



  # ==============================================================================
  # Step 5.2: Prepare ggplot objects
  # Next: Extract ggplot objects and compute shared axis limits.
  # ==============================================================================
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
  xlim_union <- range(unlist(lapply(lims_list, function(l) l$x)), na.rm = TRUE)
  ylim_union <- range(unlist(lapply(lims_list, function(l) l$y)), na.rm = TRUE)

  # ==============================================================================
  # Step 5.3: Save reference plots (aplpack)
  # Next: Create aplpack bagplots and write them to disk.
  # ==============================================================================
  dir.create(file.path(save_dir, paste0("simu_", dist_type)), recursive = TRUE)
  png(file.path(save_dir, paste0("simu_", dist_type), paste0(fig_title, "_aplpack.png")),
    width = 3, height = 3, units = "in", res = 300
  )
  par(mar = c(3.5, 3.5, 3.3, 1.5), mgp = c(1.5, 0.5, 0))
  bo_apl <- aplpack::bagplot(dat[1:10, ],
    factor = 3, type1 = "unadjusted", normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "F",
    show.whiskers = FALSE, show.loophull = TRUE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE,
    show.fence_mag_bag = TRUE, xlim = xlim_union, ylim = ylim_union
  )
  dev.off()
  png(file.path(save_dir, paste0("simu_", dist_type), paste0(fig_title, "_aplpack.png")),
    width = 3, height = 3, units = "in", res = 300
  )
  par(mar = c(3.5, 3.5, 3.3, 1.5), mgp = c(1.5, 0.5, 0))
  bo_apl <- aplpack::bagplot(dat,
    factor = 3, type1 = "unadjusted", normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = TRUE, approx.limit = 10000, redefine_loop = redefine_loop,
    show.outlier = TRUE, show.looppoints = TRUE, whisker.end.prop = 0.7,
    show.bagpoints = TRUE, dkmethod = 1, asymp_dist_pv = "F",
    show.whiskers = FALSE, show.loophull = TRUE, center_type = center_type,
    show.baghull = TRUE, verbose = FALSE,
    show.fence_mag_bag = TRUE, xlim = xlim_union, ylim = ylim_union
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
}


# ==============================================================================
# Step 6: Combine panels into a single figure
# Next: Combine panels into a single figure.
# ==============================================================================
read_img_grob <- function(path) {
  img <- readPNG(path)
  rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
}

img_prefixes <- c(
  "./figures/simu_corr_normal_mixture/corr_normal_mixture_center_hdepth_"
)
img_names <- c(
  "corr_normal_mixture_center_hdepth_combined.png"
)

for (i in 1:length(img_prefixes)) {
  img_prefix <- img_prefixes[i]
  img_paths <- c(
    paste0(img_prefix, "aplpack.png"),
    paste0(img_prefix, "sub", 1:6, ".png")
  )

  grobs <- c(lapply(img_paths, read_img_grob))

  labels <- c(
    "(a) Unadjusted\n      (Rousseeuw et al. 1999)",
    "(b) Chi-square approximation\n      FWER, q=0.1",
    "(c) Chi-square approximation\n      FDR, q=0.1",
    "(d) Chi-square approximation\n      PFER, q=0.5",
    "(e) F approximation\n      FWER, q=0.1",
    "(f) F approximation\n     FDR, q=0.1",
    "(g) F approximation\n      PFER, q=0.5"
  )

  grobs <- append(grobs, list(grid::nullGrob()), after = 4)
  labels <- append(labels, "", after = 4)

  labeled <- mapply(function(g, lb) {
    ggdraw() +
      draw_grob(g, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label(lb,
        x = 0, y = 0.98, hjust = 0, vjust = 1,
        fontface = "bold", size = 26, color = "black"
      )
  }, grobs, labels, SIMPLIFY = FALSE)

  combined <- plot_grid(plotlist = labeled, ncol = 4)
  ggsave(
    filename = file.path(save_dir, img_names[i]),
    plot = combined,
    width = 24, height = 12, dpi = 300, bg = "white"
  )
}
