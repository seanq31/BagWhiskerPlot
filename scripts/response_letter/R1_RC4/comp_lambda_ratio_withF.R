# ==============================================================================
# Step 1: Setup
# Next: Set the random seed and load required packages.
# ==============================================================================
seed_base <- 1234
set.seed(seed_base)

n_reps <- 500

library(foreach)      # makes %dopar% available
library(doParallel)   # optional, but consistent with registerDoParallel()
library(doSNOW)   # <- add
library(ggplot2)
library(cowplot)
library(grid)
library(jpeg)


# ==============================================================================
# Step 2: Configure output paths
# Next: Create the output directory for tables and figures.
# ==============================================================================
save_dir <- "./tables/"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Step 3: Define simulation grid
# Next: Specify the parameter values used in the sensitivity study.
# ==============================================================================
# Grid for lambda_data/lambda_stat sensitivity
p_values <- c(0.01, 0.02, 0.03, 0.04, 0.05)
n_values <- c(20, 50, 100, 500, 1000)

# ==============================================================================
# Step 4: Define data-generation helper
# Next: Provide a single function that simulates data for each distribution.
# ==============================================================================
generate_data <- function(dist_type, n, p, corr_coef = 0.3) {
  dist_type <- as.character(dist_type)
  if (dist_type == "normal_mixture") {
    theta <- rbinom(n, 1, p)
    x <- (1 - theta) * rnorm(n, 100) + theta * rnorm(n, 103)
    y <- (1 - theta) * rnorm(n, 300) + theta * rnorm(n, 298)
    return(cbind(x, y))
  }
  if (dist_type == "corr_normal_mixture") {
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
    return(dat_normal)
  }
  if (dist_type == "chisq") {
    x <- rchisq(n, 10)
    y <- rchisq(n, 10)
    return(cbind(x, y))
  }
  if (dist_type == "t_dist") {
    x <- rt(n, df = 10)
    y <- rt(n, df = 10)
    return(cbind(x, y))
  }
  if (dist_type == "log_normal") {
    x <- rlnorm(n, meanlog = 0, sdlog = 0.5)
    y <- rlnorm(n, meanlog = 0, sdlog = 0.5)
    return(cbind(x, y))
  }
  if (dist_type == "cauchy") {
    x <- rcauchy(n, location = 0, scale = 1)
    y <- rcauchy(n, location = 0, scale = 1)
    return(cbind(x, y))
  }
  stop("Unknown dist_type: ", dist_type)
}

# ==============================================================================
# Step 5: Define experiment settings
# Next: Create the full settings grid used in the parallel experiment.
# ==============================================================================
dist_settings <- c("log_normal", "normal_mixture", "corr_normal_mixture")
inlier_settings <- c(FALSE)
outter_settings <- c(TRUE)
redefine_loop_settings <- c(TRUE)
conservative_lambda_settings <- c(1)
center_type_settings <- c("hdepth")
corr_coef_values <- c(0.3)
asymp_dist_pv_values <- c("chisq", "F")

all_settings <- expand.grid(
  dist = dist_settings,
  n = n_values,
  p = p_values,
  corr_coef = corr_coef_values,
  normal_inlier = inlier_settings,
  normal_outter = outter_settings,
  conservative_lambda = conservative_lambda_settings,
  redefine_loop = redefine_loop_settings,
  center_type = center_type_settings,
  asymp_dist_pv = asymp_dist_pv_values
)

# ==============================================================================
# Step 6: Create the parallel task grid
# Next: Expand settings by repetitions to form independent tasks.
# ==============================================================================
# Run in parallel over (setting_idx, rep_idx)
task_grid <- expand.grid(
  setting_idx = seq_len(nrow(all_settings)),
  rep_idx = seq_len(n_reps)
)


n_tasks <- nrow(task_grid)

# ==============================================================================
# Step 7: Initialize parallel backend and progress reporting
# Next: Start a cluster, register the backend, and configure a progress bar.
# ==============================================================================
cl <- parallel::makeCluster(80)
doSNOW::registerDoSNOW(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
}, add = TRUE)

pb <- txtProgressBar(min = 0, max = n_tasks, style = 3)
on.exit(close(pb), add = TRUE)

progress_fun <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress_fun)

# ==============================================================================
# Step 8: Run simulations in parallel
# Next: For each task, simulate data, run the method, and collect summary metrics.
# ==============================================================================
results <- foreach::foreach(
  task_idx = seq_len(n_tasks),
  .combine = rbind,
  .packages = c("BagWhiskerPlot", "MASS", "parallel"),
  .options.snow = opts          # <- add
) %dopar% {
  setting_idx <- task_grid$setting_idx[[task_idx]]
  rep_idx <- task_grid$rep_idx[[task_idx]]

  dist_type <- all_settings[setting_idx, "dist"]
  n <- all_settings[setting_idx, "n"]
  p <- all_settings[setting_idx, "p"]
  corr_coef <- all_settings[setting_idx, "corr_coef"]
  normal_inlier <- all_settings[setting_idx, "normal_inlier"]
  normal_outter <- all_settings[setting_idx, "normal_outter"]
  conservative_lambda <- all_settings[setting_idx, "conservative_lambda"]
  redefine_loop <- all_settings[setting_idx, "redefine_loop"]
  center_type <- all_settings[setting_idx, "center_type"]
  asymp_dist_pv <- all_settings[setting_idx, "asymp_dist_pv"]

  seed <- seed_base + (setting_idx - 1) * n_reps + rep_idx
  set.seed(seed)
  dat <- generate_data(dist_type = dist_type, n = n, p = p, corr_coef = corr_coef)

  old <- Sys.time()
  tmp_plt <- BagWhiskerPlot::bag_whisker(dat,
    type1 = "FDR", q = 0.1, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, create.plot = FALSE, approx.limit = 10000, redefine_loop = redefine_loop, cex = 0.6, pch = 1,
    show.outlier = TRUE, show.loophull = FALSE,
    show.bagpoints = TRUE, dkmethod = 1, center_type = center_type,
    show.whiskers = TRUE, asymp_dist_pv = asymp_dist_pv, whisker.end.prop = 0.7,
    show.baghull = TRUE, verbose = FALSE, col.baghull = "#D3D3D3",
    show.fence_mag_bag = TRUE, n_cores = 1
  )
  elapsed_secs <- as.numeric(difftime(Sys.time(), old, units = "secs"))

  bo <- tmp_plt$bo
  data.frame(
    dist = as.character(dist_type),
    n = as.integer(n),
    p = as.numeric(p),
    corr_coef = as.numeric(corr_coef),
    center_type = as.character(center_type),
    normal_inlier = as.logical(normal_inlier),
    normal_outter = as.logical(normal_outter),
    conservative_lambda = as.numeric(conservative_lambda),
    redefine_loop = as.logical(redefine_loop),
    asymp_dist_pv = as.character(asymp_dist_pv),
    rep = as.integer(rep_idx),
    seed = as.integer(seed),
    lambda_data = as.numeric(bo$lambda_data),
    lambda_stat = as.numeric(bo$lambda_stat),
    lambda_data_over_stat = as.numeric(bo$lambda_data_over_stat),
    elapsed_secs = elapsed_secs
  )
}

# ==============================================================================
# Step 9: Export results
# Next: Save the long-format results and per-distribution summaries/heatmaps.
# ==============================================================================
outfile_long <- file.path(save_dir, "lambda_data_over_stat_long.csv")
write.csv(results, outfile_long, row.names = FALSE)

# summary_settings is the product of dist_settings and asymp_dist_pv_values
summary_settings <- expand.grid(
  dist = dist_settings,
  asymp_dist_pv = asymp_dist_pv_values
)

for (summary_setting_idx in seq_len(nrow(summary_settings))) {
  summary_setting <- summary_settings[summary_setting_idx, ]
  dist_results <- subset(results, dist == summary_setting$dist & asymp_dist_pv == summary_setting$asymp_dist_pv)
  dist_summary <- aggregate(lambda_data_over_stat ~ n + p, data = dist_results, FUN = mean)
  dist_summary_wide <- reshape(
    dist_summary,
    idvar = "n",
    timevar = "p",
    direction = "wide"
  )
  dist_summary_wide <- dist_summary_wide[order(dist_summary_wide$n), ]
  outfile_dist <- file.path(save_dir, paste0("lambda_data_over_stat_summary_", summary_setting$dist, "_", summary_setting$asymp_dist_pv, ".csv"))
  write.csv(dist_summary_wide, outfile_dist, row.names = FALSE)

  dist_summary$p <- as.factor(dist_summary$p)
  p_heatmap <- ggplot(dist_summary, aes(x = p, y = factor(n), fill = lambda_data_over_stat)) +
    geom_tile() +
    scale_fill_gradient2(low = "#3B5B92", mid = "#D0D0D0", high = "#9C3A3A", midpoint = 0.85, limits = c(0.69, 1.0)) +
    labs(x = "p", y = "n", fill = expression(lambda[plain(data)]/lambda[plain(stat)])) +
    theme_minimal()
  outfile_heatmap <- file.path(save_dir, paste0("lambda_data_over_stat_heatmap_", summary_setting$dist, "_", summary_setting$asymp_dist_pv, ".jpeg"))
  ggsave(outfile_heatmap, p_heatmap, width = 6, height = 4)
}


read_img_grob <- function(path) {
  img <- readJPEG(path)
  rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
}

for (dist_setting in dist_settings) {
  img_paths <- c(
    paste0(save_dir, "lambda_data_over_stat_heatmap_", dist_setting, "_chisq.jpeg"),
    paste0(save_dir, "lambda_data_over_stat_heatmap_", dist_setting, "_F.jpeg")
  )

  grobs <- c(lapply(img_paths, read_img_grob))
  labels <- c(
      "(a) Chi-square approximation, FDR, q=0.1",
      "(b) F approximation, FDR, q=0.1"
  )
  labeled <- mapply(function(g, lb) {
  ggdraw() +
      draw_grob(g, x = 0, y = 0, width = 1, height = 0.92) +
      draw_label(lb,
      x = 0.02, y = 0.98, hjust = 0, vjust = 1,
      fontface = "bold", size = 18, color = "black"
      )
  }, grobs, labels, SIMPLIFY = FALSE)

  combined <- plot_grid(plotlist = labeled, ncol = 2)
  ggsave(
  filename = file.path(save_dir, paste0("lambda_data_over_stat_heatmap_", dist_setting, ".jpeg")),
  plot = combined,
  width = 12, height = 5, dpi = 300, bg = "white"
  )
}

