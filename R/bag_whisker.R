#' Bag-and-whisker plot implementation, with outlier detection 
#' under different type-I error controls
#'
#' @description
#' Compute and draw a bag-and-whisker plot of bivariate data that
#' highlights outliers under different type-I error controls. This is
#' a thin user-facing wrapper around the internal computation and
#' ggplot2-based plotting helpers in this package.
#'
#' @param x,y Numeric vectors giving the coordinates of the bivariate sample.
#'   You can either supply both `x` and `y`, or a two-column matrix / data
#'   frame via `x` and leave `y` missing.
#' @param type1 Character string specifying the type-I error control used in
#'   the multiple-testing step. One of `"unadjusted"`, `"FWER"`, `"FDR"`,
#'   or `"PFER"`, default `"unadjusted"`.
#' @param q Numeric, the control level for the multiple testing procedure, default `0.1`.
#' @param normal_inlier Logical; if `TRUE`, the bag is calculated based on a
#'   normal reference distribution instead of a depth-based bag, default `FALSE`.
#' @param normal_outter Logical; if `TRUE`, the fence is calculated based on a
#'   normal reference distribution, otherwise it is obtained by magnifying
#'   the bag, default `FALSE`.
#' @param asymp_dist_pv Character string; asymptotic distribution used for
#'   p-value calculation, default `"chisq"`.
#' @param center_type Character string; center definition, default `"hdepth"`.
#' @param factor Numeric; the factor lambda for `"unadjusted"` that controls
#'   how far the bag is magnified to obtain the fence, default `3`.
#' @param na.rm Logical; if `TRUE`, rows with missing values are removed,
#'   otherwise medians are used to impute, default `FALSE`.
#' @param approx.limit Integer; threshold above which a subsample is used for
#'   approximating the bag-plot computation, default `300`.
#' @param dkmethod Integer in `1:2`; depth kernel method, default `1`.
#' @param precision Numeric; controls precision of hull expansion, default `1`.
#' @param n_cores Integer or NULL. If greater than 1 or NULL, uses a parallel
#'   backend (via [parallel::makeCluster()]) for some calculation intensive
#'   loops. If `NULL`, uses `parallel::detectCores() - 1` cores
#'   (minimum 1). Default `1`.
#' @param verbose Logical; if `TRUE`, progress messages from the computational
#'   engine are printed, default `FALSE`.
#' @param timing Logical; if `TRUE`, prints a simple timing summary and adds a
#'   `timings` component (named numeric vector, elapsed seconds) to the returned
#'   object, default `FALSE`.
#' @param debug.plots Character string; controls generation of additional
#'   diagnostic plots for debugging, default `no`.
#' @param show.outlier Logical; if `TRUE`, identified outliers are shown in
#'   the plot, default `TRUE`.
#' @param show.whiskers Logical; if `TRUE`, whisker segments are added to the
#'   plot, default `TRUE`.
#' @param show.looppoints Logical; if `TRUE`, data points classified as
#'   belonging to the loop are shown, default `TRUE`.
#' @param show.bagpoints Logical; if `TRUE`, data points classified as inside
#'   the bag are shown, default `TRUE`.
#' @param show.loophull Logical; if `TRUE`, the convex hull of the loop is
#'   drawn, default `FALSE`.
#' @param show.baghull Logical; if `TRUE`, the convex hull of the bag is
#'   drawn, default `TRUE`.
#' @param create.plot Logical; if `TRUE`, a plot is created for the
#'   computed bag-and-whisker representation, default `TRUE`.
#' @param add Logical; if `TRUE`, graphical elements are added to the current
#'   plot; otherwise a new plot is started, default `FALSE`.
#' @param pch,cex Graphical parameters forwarded to the plotting method to
#'   control point character and expansion.
#' @param transparency Logical; if `TRUE`, semi-transparent fills are used for
#'   the bag and loop hulls, default `FALSE`.
#' @param col.loophull Fill colour used for the
#'   loop hull.
#' @param col.looppoints Point colour used for the
#'   loop points.
#' @param col.baghull Fill colour used for the bag
#'   hull.
#' @param col.bagpoints Point colour used for the bag
#'   points.
#' @param show.center Logical; if `TRUE`, the chosen center of the data is
#'   highlighted in the plot, default `TRUE`.
#' @param show.fence_mag_bag Logical; if `TRUE`, the fence by magnifying
#'   the bag is visualised, default `TRUE`.
#' @param ... Passed on to the S3 method [plot.bagWhiskerPlot()] when
#'   `create.plot = TRUE`.
#'
#' @details
#' This function is the only exported user-facing entry point of the
#' package. It calls the internal computational engine
#' `compute.bagWhiskerPlot()` and then, by default, plots the resulting
#' object using a ggplot2-based implementation of the S3 method
#' [plot.bagWhiskerPlot()].
#'
#' @return
#' An object of class `"bagWhiskerPlot"` containing the bag-plot
#' decomposition and outlier information. When `show = TRUE`, the object
#' is returned invisibly after drawing the plot.
#'
#' @export
bag_whisker <- function(
    x, y,
    type1 = "unadjusted", # type 1 error notion, other options are FDR, FWER, PFER
    q = 0.1, # the control level for the type1 error notion
    normal_inlier = FALSE,
    normal_outter = FALSE, # shape of the outerfence, according to normal or multiple of inner bag
    asymp_dist_pv = "chisq",
    center_type = "hdepth",
    factor = 3, # expanding factor for bag to get the loop
    na.rm = FALSE, # should 'NAs' values be removed or exchanged
    approx.limit = 300, # limit
    show.outlier = TRUE, # if TRUE outlier are shown
    show.whiskers = TRUE, # if TRUE whiskers are shown
    show.looppoints = TRUE, # if TRUE points in loop are shown
    show.bagpoints = TRUE, # if TRUE points in bag are shown
    show.loophull = FALSE, # if TRUE loop is shown
    show.baghull = TRUE, # if TRUE bag is shown
    create.plot = TRUE, # if TRUE a plot is created
    add = FALSE, # if TRUE graphical elements are added to actual plot
    pch = 1, cex = 0.6, # some graphical parameters
    dkmethod = 1, # in 1:2; there are two methods for approximating the bag
    precision = 1, # controls precision of computation
    n_cores = 1,
    verbose = FALSE, debug.plots = "no", # tools for debugging
    timing = FALSE,
    col.loophull = "#aaccff", # Alternatives: #ccffaa, #ffaacc
    col.looppoints = "#3355ff", # Alternatives: #55ff33, #ff3355
    col.baghull = "#D3D3D3", # Alternatives: #99ff77, #ff7799
    col.bagpoints = "#000088", # Alternatives: #008800, #880000
    transparency = FALSE,
    show.center = TRUE, # if TRUE center is shown
    show.fence_mag_bag = TRUE,
    ... # to define further parameters of plot
    ) {
  if (missing(x)) {
    return()
  }
  if ((is.data.frame(x) || is.matrix(x)) && ncol(x) == 2) {
    y <- x[, 2]
    x <- x[, 1]
  } # 180308
  if (missing(y)) {
    y <- x
    x <- seq(along = y)
  } # 180308
  bo <- compute.bagWhiskerPlot(
    x = x, y = y, normal_inlier = normal_inlier, normal_outter = normal_outter, type1 = type1, q = q, factor = factor, na.rm = na.rm, asymp_dist_pv = asymp_dist_pv,
    approx.limit = approx.limit, dkmethod = dkmethod, center_type = center_type,
    precision = precision, n_cores = n_cores, verbose = verbose, debug.plots = debug.plots, timing = timing
  )
  if (create.plot) {
    tmp_plt <- plot(bo,
      show.outlier = show.outlier,
      show.whiskers = show.whiskers,
      show.looppoints = show.looppoints,
      show.bagpoints = show.bagpoints,
      show.loophull = show.loophull,
      show.baghull = show.baghull,
      show.fence_mag_bag = show.fence_mag_bag,
      add = add, pch = pch, cex = cex,
      verbose = verbose,
      col.loophull = col.loophull,
      col.looppoints = col.looppoints,
      col.baghull = col.baghull,
      col.bagpoints = col.bagpoints,
      transparency = transparency,
      show.center = show.center, ...
    )
    invisible(list(bo = bo, bgplt = tmp_plt))
  } else {
    invisible(list(bo = bo, add_args = list(
      show.outlier = show.outlier,
      show.whiskers = show.whiskers,
      show.looppoints = show.looppoints,
      show.bagpoints = show.bagpoints,
      show.loophull = show.loophull,
      show.baghull = show.baghull,
      add = add, pch = pch, cex = cex,
      verbose = verbose,
      col.loophull = col.loophull,
      col.looppoints = col.looppoints,
      col.baghull = col.baghull,
      col.bagpoints = col.bagpoints,
      transparency = transparency,
      show.center = show.center,
      ...
    )))
  }
}
