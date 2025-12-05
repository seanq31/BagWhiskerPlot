#' Bag-Whisker Plot with Multiple-Testing Based Outlier Control
#'
#' @description
#' Compute a bag-plot based representation of bivariate data and draw
#' whisker-like segments that highlight outliers under different type-I
#' error controls. This is a thin user-facing wrapper around the internal
#' computation and ggplot2-based plotting helpers in this package.
#'
#' @param x,y Numeric vectors giving the coordinates of the bivariate sample.
#'   You can either supply both `x` and `y`, or a two-column matrix / data
#'   frame via `x` and leave `y` missing.
#' @param type1 Character string specifying the type-I error control used in
#'   the multiple-testing step. One of `"unadjusted"`, `"FWER"`, `"FDR"`,
#'   or `"PFER"`.
#' @param q Numeric, the control level for the type-I error measure.
#' @param factor Numeric, expansion factor used to construct the outer loop
#'   from the central bag.
#' @param normal_inlier,normal_outter Logical flags for using normal-theory
#'   based calibration for inner / outer regions (advanced; see paper).
#' @param conservative_lambda Numeric tuning parameter used in the
#'   multiple-testing adjustment.
#' @param redefine_loop Logical; if `TRUE`, the loop (outer hull) may be
#'   recomputed after multiple-testing.
#' @param asymp_dist_pv Character string; asymptotic distribution used for
#'   p-value calculation, default `"chisq"`.
#' @param center_type Character string; center definition, default `"hdepth"`.
#' @param na.rm Logical; if `TRUE`, rows with missing values are removed,
#'   otherwise medians are used to impute.
#' @param approx.limit Integer; threshold above which a subsample is used for
#'   approximating the bag-plot computation.
#' @param dkmethod Integer in `1:2`; depth kernel method, with `2` recommended.
#' @param precision Numeric; controls precision of hull expansion.
#' @param show Logical; if `TRUE` (default), a plot is produced via
#'   [plot.bagWhiskerPlot()]. If `FALSE`, only the computed object is
#'   returned.
#' @param ... Passed on to [plot.bagWhiskerPlot()] when `show = TRUE`.
#'
#' @details
#' This function is the only exported user-facing entry point of the
#' package. It calls the internal computational engine
#' [compute.bagWhiskerPlot()] and then, by default, plots the resulting
#' object using a ggplot2-based implementation of
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
    conservative_lambda = 0,
    asymp_dist_pv = 'chisq',
    center_type = 'hdepth',
    redefine_loop = TRUE, # if TRUE loop is redefined by fence_mag_bag
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
    verbose = FALSE, debug.plots = "no", # tools for debugging
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
    x = x, y = y, normal_inlier = normal_inlier, normal_outter = normal_outter, conservative_lambda = conservative_lambda, redefine_loop = redefine_loop, type1 = type1, q = q, factor = factor, na.rm = na.rm, asymp_dist_pv = asymp_dist_pv,
    approx.limit = approx.limit, dkmethod = dkmethod, center_type = center_type,
    precision = precision, verbose = verbose, debug.plots = debug.plots
  )
  if (create.plot) {
    tmp_plt = plot(bo,
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