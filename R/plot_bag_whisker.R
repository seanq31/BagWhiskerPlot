# Declare global variables to silence R CMD check notes from ggplot2's NSE
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "y", "xend", "yend", "group", "alpha"))
}

# Internal: compute polar angle
.bp_win <- function(dx, dy) {
  atan2(y = dy, x = dx)
}

# Internal: compute intersection of ray z->(infty) with polygon edge (p1, p2)
.bp_cut_z_pg <- function(zx, zy, p1x, p1y, p2x, p2y, verbose = FALSE, debug.plots = "no") {
  a2 <- (p2y - p1y) / (p2x - p1x)
  a1 <- zy / zx
  sx <- (p1y - a2 * p1x) / (a1 - a2)
  sy <- a1 * sx
  sxy <- cbind(sx, sy)
  h <- any(is.nan(sxy)) || any(is.na(sxy)) || any(Inf == abs(sxy))
  if (h) {
    if (verbose) cat("special")
    h0 <- 0 == zx
    sx <- ifelse(h0, zx, sx)
    sy <- ifelse(h0, p1y - a2 * p1x, sy)
    a1 <- ifelse(abs(a1) == Inf, sign(a1) * 123456789 * 1E10, a1)
    a2 <- ifelse(abs(a2) == Inf, sign(a2) * 123456789 * 1E10, a2)
    h <- 0 == (a1 - a2) & sign(zx) == sign(p1x)
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, p1y, sy)
    h <- 0 == (a1 - a2) & sign(zx) != sign(p1x)
    sx <- ifelse(h, p2x, sx); sy <- ifelse(h, p2y, sy)
    h <- p1x == p2x & zx != p1x & p1x != 0
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, zy * p1x / zx, sy)
    h <- p1x == p2x & zx != p1x & p1x == 0
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, 0, sy)
    h <- p1x == p2x & zx == p1x & p1x != 0
    sx <- ifelse(h, zx, sx);  sy <- ifelse(h, zy, sy)
    h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) == sign(p1y)
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, p1y, sy)
    h <- p1x == p2x & zx == p1x & p1x == 0 & sign(zy) != sign(p1y)
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, p2y, sy)
    h <- zx == p1x & zy == p1y
    sx <- ifelse(h, p1x, sx); sy <- ifelse(h, p1y, sy)
    h <- zx == p2x & zy == p2y
    sx <- ifelse(h, p2x, sx); sy <- ifelse(h, p2y, sy)
    h <- zx == 0 & zy == 0
    sx <- ifelse(h, 0, sx);   sy <- ifelse(h, 0, sy)
    sxy <- cbind(sx, sy)
  }
  if (debug.plots == "all") {
    segments(sxy[, 1], sxy[, 2], zx, zy, col = "red")
    segments(0, 0, sxy[, 1], sxy[, 2], col = "green", lty = 2)
    points(sxy, col = "red")
  }
  sxy
}

# Internal: for each z, find intersection with polygon pg (centered around center)
.bp_find_cut_z_pg <- function(z, pg, center = c(0, 0), debug.plots = "no") {
  if (!is.matrix(z)) z <- rbind(z)
  if (is.null(pg) || is.vector(pg) || 1 == nrow(pg)) {
    return(matrix(center, nrow(z), 2, TRUE))
  }
  n.pg <- nrow(pg)
  zc <- cbind(z[, 1] - center[1], z[, 2] - center[2])
  pgc <- cbind(pg[, 1] - center[1], pg[, 2] - center[2])
  apg <- .bp_win(pgc[, 1], pgc[, 2]); apg[is.nan(apg)] <- 0
  ord <- order(apg)
  apg <- apg[ord]; pgc <- pgc[ord, ]
  az <- .bp_win(zc[, 1], zc[, 2])
  segm.no <- apply(outer(apg, az, "<"), 2, sum)
  segm.no <- ifelse(segm.no == 0, n.pg, segm.no)
  next.no <- 1 + (segm.no %% length(apg))
  cuts <- .bp_cut_z_pg(
    zc[, 1], zc[, 2], pgc[segm.no, 1], pgc[segm.no, 2],
    pgc[next.no, 1], pgc[next.no, 2]
  )
  cuts <- cbind(cuts[, 1] + center[1], cuts[, 2] + center[2])
  cuts
}

# Internal: data.frame for polygon
.bp_poly_df <- function(poly, group) {
  if (is.null(poly) || length(poly) == 0 || is.vector(poly)) return(NULL)
  data.frame(
    x = poly[, 1],
    y = poly[, 2],
    group = group,
    stringsAsFactors = FALSE
  )
}

# Internal: data.frame for segments
.bp_segments_df <- function(x, y, xend, yend) {
  data.frame(x = x, y = y, xend = xend, yend = yend)
}

.bp_fade_segments_df <- function(x, y, xend, yend,
                                n = 20,
                                alpha_start = 0.9,
                                alpha_end = 0.05,
                                end_prop = 0.5,
                                fade_towards = c("end", "start")) {
  fade_towards <- match.arg(fade_towards)
  if (length(x) == 0) return(data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0), alpha = numeric(0)))
  n <- max(1L, as.integer(n))
  end_prop <- max(0, min(1, end_prop))
  out <- vector("list", length(x) * n)
  idx <- 1L
  for (i in seq_along(x)) {
    x0 <- x[i]; y0 <- y[i]; x1 <- xend[i]; y1 <- yend[i]
    for (k in 0:(n - 1L)) {
      t0 <- k / n; t1 <- (k + 1L) / n
      xs <- x0 + (x1 - x0) * t0
      ys <- y0 + (y1 - y0) * t0
      xe <- x0 + (x1 - x0) * t1
      ye <- y0 + (y1 - y0) * t1
      tm <- (t0 + t1) / 2
      if (tm <= end_prop && end_prop > 0) {
        g <- tm / end_prop
        if (fade_towards == "end") {
          a <- alpha_start + (alpha_end - alpha_start) * g
        } else {
          a <- alpha_end + (alpha_start - alpha_end) * g
        }
      } else {
        a <- alpha_end
      }
      out[[idx]] <- data.frame(x = xs, y = ys, xend = xe, yend = ye, alpha = a)
      idx <- idx + 1L
    }
  }
  do.call(rbind, out)
}

# Internal: hex alpha helper. If color is hex like "#RRGGBB", append alphaHex.
.bp_apply_alpha_hex <- function(col, alphaHex = "99") {
  ifelse(grepl("^#", col), paste0(col, alphaHex),
         grDevices::adjustcolor(col, alpha = as.integer(paste0("0x", alphaHex)) / 255))
}


.bp_points_outside_poly <- function(pts, poly, tol = 1e-8) {
  # Robust point-in-polygon test (ray casting), treating boundary points as INSIDE
  # and guarding against numerical problems on (almost) horizontal edges.

  if (is.null(poly) || length(poly) == 0) {
    # No polygon => nothing is considered outside for sizing purposes
    if (!is.matrix(pts)) pts <- rbind(pts)
    return(rep(FALSE, nrow(pts)))
  }

  if (!is.matrix(pts)) pts <- as.matrix(pts)
  if (!is.matrix(poly)) poly <- as.matrix(poly)

  x <- pts[, 1]; y <- pts[, 2]
  xv <- poly[, 1]; yv <- poly[, 2]
  n <- length(xv)
  inside <- logical(length(x))

  # Helper: distance of point to segment (for boundary detection)
  point_on_segment <- function(px, py, x1, y1, x2, y2, tol) {
    # Bounding box quick rejection (with tolerance)
    if (px < min(x1, x2) - tol || px > max(x1, x2) + tol ||
        py < min(y1, y2) - tol || py > max(y1, y2) + tol) {
      return(FALSE)
    }
    vx <- x2 - x1; vy <- y2 - y1
    wx <- px - x1; wy <- py - y1
    # Cross-product magnitude (area ~ 0 for collinear)
    cross <- vx * wy - vy * wx
    if (abs(cross) > tol * (abs(vx) + abs(vy) + 1)) return(FALSE)
    # Now check projection lies within segment
    dot <- vx * wx + vy * wy
    if (dot < -tol) return(FALSE)
    if (dot > (vx * vx + vy * vy) + tol) return(FALSE)
    TRUE
  }

  j <- n
  for (i in seq_len(n)) {
    xi <- xv[i]; yi <- yv[i]
    xj <- xv[j]; yj <- yv[j]

    # First: mark points that lie (numerically) on this edge as inside
    on_edge <- mapply(
      point_on_segment,
      px = x, py = y,
      MoreArgs = list(x1 = xi, y1 = yi, x2 = xj, y2 = yj, tol = tol)
    )
    inside[on_edge] <- TRUE

    # Skip ray-crossing contribution for horizontal edges (avoid division by ~0)
    dy_edge <- yj - yi
    if (abs(dy_edge) < tol) {
      j <- i
      next
    }

    # Standard ray-casting test for edges that cross horizontal ray to the right
    # Use half-open interval on upper endpoint to avoid double-counting vertices.
    y_low <- pmin(yi, yj)
    y_high <- pmax(yi, yj)

    crosses <- (y > y_low) & (y <= y_high)
    # Intersection x-coordinate
    x_int <- xi + (y - yi) * (xj - xi) / dy_edge

    intersect <- crosses & (x <= x_int + tol)
    inside <- xor(inside, intersect)

    j <- i
  }

  # outside = NOT inside
  !inside
}

# Build layers only (for advanced composition)
bp_build_layers <- function(
  x,
  show.outlier = TRUE,
  show.whiskers = TRUE,
  show.looppoints = TRUE,
  show.bagpoints = TRUE,
  show.loophull = TRUE,
  show.baghull = TRUE,
  show.fence_mag_bag = TRUE,
  pch = 16, cex = .4,
  col.loophull = "#aaccff",
  col.looppoints = "#3355ff",
  col.baghull = "#7799ff",
  col.bagpoints = "#000088",
  col.fence_mag_bag = "#CC33CC",
  transparency = FALSE,
  show.center = TRUE,
  whisker.fade = TRUE,
  whisker.n = 10,
  whisker.alpha.start = 0.4,
  whisker.alpha.end = 0.0,
  whisker.end.prop = 0.7
) {
  # unpack bagplot object (same as original code)
  center <- hull.center <- hull.bag <- hull.loop <- fence_mag_bag <- pxy.bag <- pxy.outer <- pxy.outlier <- NULL
  hdepths <- is.one.dim <- prdata <- xy <- xydata <- exp.dk <- exp.dk.1 <- hdepth <- NULL
  tphdepth <- tp <- NULL
  bagplotobj <- x
  for (i in seq(along = bagplotobj)) {
    eval(parse(text = paste(names(bagplotobj)[i], "<-bagplotobj[[", i, "]]")))
  }

  # colors with transparency if requested
  if (isTRUE(transparency)) {
    col.loophull <- .bp_apply_alpha_hex(col.loophull, "99")
    col.baghull <- .bp_apply_alpha_hex(col.baghull, "99")
  }

  # We'll collect whiskers separately so they can be forced to the bottom layer
  whisker_layers <- list()
  layers <- list()

  # One-dimensional special case: emulate original 1D boxplot logic using ggplot layers
  if (isTRUE(is.one.dim)) {
    ROT <- round(prdata[[2]], digits = 5)
    IROT <- round(solve(ROT), digits = 5)

    # compute 5-number summary in projected axis and map back
    usrX <- range(xydata[, 1], na.rm = TRUE)
    usrY <- range(xydata[, 2], na.rm = TRUE)

    if (ROT[1, 1] == 0) {
      xy2 <- cbind(mean(usrX), xydata[, 2])
      bpr <- boxplot(xy2[, 2], plot = FALSE)
      five <- cbind(mean(usrX), bpr$stat)
      dx <- 0.1 * diff(usrX); dy <- 0
      idx.out <- if (length(bpr$out)) match(bpr$out, xy2[, 2]) else integer(0)
      out.df <- if (length(idx.out)) data.frame(x = xy2[idx.out, 1], y = xy2[idx.out, 2]) else NULL
    } else if (ROT[1, 2] == 0) {
      xy2 <- cbind(xydata[, 1], mean(usrY))
      bpr <- boxplot(xy2[, 1], plot = FALSE)
      five <- cbind(bpr$stat, mean(usrY))
      dx <- 0; dy <- 0.1 * diff(usrY)
      idx.out <- if (length(bpr$out)) match(bpr$out, xy2[, 1]) else integer(0)
      out.df <- if (length(idx.out)) data.frame(x = xy2[idx.out, 1], y = xy2[idx.out, 2]) else NULL
    } else {
      xytr <- xydata %*% ROT
      bpr <- boxplot(xytr[, 1], plot = FALSE)
      five <- cbind(bpr$stat, xytr[1, 2]) %*% IROT
      vec <- five[5, ] - five[1, ]
      vec.ortho <- c(vec[2], -vec[1])
      # scale orthogonal vector roughly to 15% of smaller data range (approx original)
      rng <- c(diff(usrX), diff(usrY)); s <- 0.15 * min(rng) / sqrt(sum(vec.ortho^2))
      xy.delta <- vec.ortho * s
      dx <- xy.delta[1]; dy <- xy.delta[2]
      out.df <- NULL # original couldn't reliably back-map outliers here
    }

    # Whiskers segments (min-max to hinges) — add FIRST (bottom-most)
    wh.df <- rbind(
      .bp_segments_df(five[1, 1], five[1, 2], five[2, 1], five[2, 2]),
      .bp_segments_df(five[5, 1], five[5, 2], five[4, 1], five[4, 2])
    )
    wh_layer <- ggplot2::geom_segment(data = wh.df,
                                      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                                      color = "red", linewidth = 0.6)
    attr(wh_layer, "bp_role") <- "whisker"
    whisker_layers <- c(whisker_layers, list(wh_layer))

    # Box polygon (Q1-Q3)
    box.poly <- rbind(
      c(five[2, 1] + dx, five[2, 2] + dy),
      c(five[4, 1] + dx, five[4, 2] + dy),
      c(five[4, 1] - dx, five[4, 2] - dy),
      c(five[2, 1] - dx, five[2, 2] - dy)
    )
    box.df <- .bp_poly_df(box.poly, "bag_1d")
    layers <- c(layers, list(
      ggplot2::geom_polygon(data = box.df, ggplot2::aes(x = x, y = y, group = group),
                            fill = col.baghull, color = "black", linewidth = 0.5)
    ))
    # Endpoints at whiskers
    layers <- c(layers, list(
      ggplot2::geom_point(data = data.frame(x = five[c(1, 5), 1], y = five[c(1, 5), 2]),
                          ggplot2::aes(x = x, y = y),
                          color = col.looppoints, shape = 16, size = 1.2)
    ))
    # Median segment
    layers <- c(layers, list(
      ggplot2::geom_segment(
        data = .bp_segments_df(five[3, 1] + dx, five[3, 2] + dy, five[3, 1] - dx, five[3, 2] - dy),
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "red", linewidth = 1.2
      )
    ))
    # Outliers (if available)
    if (!is.null(out.df) && nrow(out.df)) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = out.df, ggplot2::aes(x = x, y = y),
                            color = "red", shape = 16, size = 1.2)
      ))
    }
    # Center point (median along main axis not defined; keep geometric center)
    if (show.center) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = data.frame(x = center[1], y = center[2]),
                            ggplot2::aes(x = x, y = y),
                            color = "red", shape = 8, size = 2)
      ))
    }

    # Ensure whiskers are at the very bottom
    return(c(whisker_layers, layers))
  }

  # 2D standard case
  # Loop hull
  if (show.loophull && !is.null(hull.loop) && length(hull.loop) > 0) {
    df.loop <- .bp_poly_df(hull.loop, "loop")
    layers <- c(layers, list(
      ggplot2::geom_polygon(data = df.loop, ggplot2::aes(x = x, y = y, group = group),
                            fill = col.loophull, color = "black", linewidth = 0.5)
    ))
  }
  # Loop points (outer, non-outliers)
  if (show.looppoints && !is.null(pxy.outer) && length(pxy.outer) > 0) {
    df.lpts <- data.frame(x = pxy.outer[, 1], y = pxy.outer[, 2])
    outside_mask <- .bp_points_outside_poly(as.matrix(df.lpts), fence_mag_bag)
    if (any(!outside_mask)) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = df.lpts[!outside_mask, , drop = FALSE], ggplot2::aes(x = x, y = y),
                            color = col.looppoints, shape = pch, size = cex)
      ))
    }
    if (any(outside_mask)) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = df.lpts[outside_mask, , drop = FALSE], ggplot2::aes(x = x, y = y),
                            color = col.looppoints, shape = pch, size = cex)
      ))
    }
  }
  # Bag hull
  if (show.baghull && !is.null(hull.bag) && length(hull.bag) > 0) {
    df.bag <- .bp_poly_df(hull.bag, "bag")
    layers <- c(layers, list(
      ggplot2::geom_polygon(data = df.bag, ggplot2::aes(x = x, y = y, group = group),
                            fill = col.baghull, color = "black", linewidth = 0.5)
    ))
  }
  # Data-adaptive magnified fence from multiple-testing threshold
  if (show.fence_mag_bag && !is.null(fence_mag_bag) && length(fence_mag_bag) > 0) {
    df.fence <- .bp_poly_df(fence_mag_bag, "fence_mag_bag")
    layers <- c(layers, list(
      ggplot2::geom_polygon(data = df.fence, ggplot2::aes(x = x, y = y, group = group),
                            fill = NA, color = col.fence_mag_bag, linewidth = 0.8, linetype = "dashed")
    ))
  }
  # Bag points (inside bag)
  if (show.bagpoints && !is.null(pxy.bag) && length(pxy.bag) > 0) {
    df.bpts <- data.frame(x = pxy.bag[, 1], y = pxy.bag[, 2])
    layers <- c(layers, list(
      ggplot2::geom_point(data = df.bpts, ggplot2::aes(x = x, y = y),
                          color = col.bagpoints, shape = pch, size = cex)
    ))
  }
  # Whiskers — collect into whisker_layers so they render underneath everything
  if (show.whiskers) {
    if (!is.null(pxy.outer) && length(pxy.outer) > 0 && !is.null(hull.bag) && length(hull.bag) > 0) {
      pkt.cut <- .bp_find_cut_z_pg(pxy.outer, hull.bag, center = center)
      if (isTRUE(whisker.fade)) {
        df.wh.fade <- .bp_fade_segments_df(
          x = pxy.outer[, 1], y = pxy.outer[, 2],
          xend = pkt.cut[, 1], yend = pkt.cut[, 2],
          n = whisker.n, alpha_start = whisker.alpha.start, alpha_end = whisker.alpha.end,
          fade_towards = "end", end_prop = whisker.end.prop
        )
        wh_layer <- ggplot2::geom_segment(
          data = df.wh.fade,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = alpha),
          color = "red", linewidth = cex, show.legend = FALSE
        )
        attr(wh_layer, "bp_role") <- "whisker"
        whisker_layers <- c(whisker_layers, list(wh_layer))
      } else {
        df.wh <- .bp_segments_df(
          x = pxy.outer[, 1], y = pxy.outer[, 2],
          xend = pkt.cut[, 1], yend = pkt.cut[, 2]
        )
        wh_layer <- ggplot2::geom_segment(data = df.wh,
                              ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                              color = "red", linewidth = 0.4)
        attr(wh_layer, "bp_role") <- "whisker"
        whisker_layers <- c(whisker_layers, list(wh_layer))
      }
    } else if (!is.null(pxy.outer) && length(pxy.outer) > 0) {
      if (isTRUE(whisker.fade)) {
        df.wh2.fade <- .bp_fade_segments_df(
          x = pxy.outer[, 1], y = pxy.outer[, 2],
          xend = rep(center[1], nrow(pxy.outer)), yend = rep(center[2], nrow(pxy.outer)),
          n = whisker.n, alpha_start = whisker.alpha.start, alpha_end = whisker.alpha.end,
          fade_towards = "end"
        )
        wh_layer <- ggplot2::geom_segment(
          data = df.wh2.fade,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = alpha),
          color = "red", linewidth = cex, show.legend = FALSE
        )
        attr(wh_layer, "bp_role") <- "whisker"
        whisker_layers <- c(whisker_layers, list(wh_layer))
      } else {
        df.wh2 <- .bp_segments_df(pxy.outer[, 1], pxy.outer[, 2], rep(center[1], nrow(pxy.outer)), rep(center[2], nrow(pxy.outer)))
        wh_layer <- ggplot2::geom_segment(data = df.wh2,
                              ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                              color = "red", linewidth = 0.4)
        attr(wh_layer, "bp_role") <- "whisker"
        whisker_layers <- c(whisker_layers, list(wh_layer))
      }
    }
  }
  # Outliers (red)
  if (show.outlier && !is.null(pxy.outlier) && length(pxy.outlier) > 0) {
    df.ol <- data.frame(x = pxy.outlier[, 1], y = pxy.outlier[, 2])
    outside_mask_ol <- .bp_points_outside_poly(as.matrix(df.ol), fence_mag_bag)
    if (any(!outside_mask_ol)) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = df.ol[!outside_mask_ol, , drop = FALSE], ggplot2::aes(x = x, y = y),
                            color = "blue", shape = pch, size = cex)
      ))
    }
    if (any(outside_mask_ol)) {
      layers <- c(layers, list(
        ggplot2::geom_point(data = df.ol[outside_mask_ol, , drop = FALSE], ggplot2::aes(x = x, y = y),
                            color = "darkred", shape = 16, size = cex * 3)
      ))
    }
  }
  # # Hull around center (orange)
  # if (!is.null(hull.center) && length(hull.center) > 2 && show.center) {
  #   df.ch <- .bp_poly_df(hull.center, "center_hull")
  #   layers <- c(layers, list(
  #     ggplot2::geom_polygon(data = df.ch, ggplot2::aes(x = x, y = y, group = group),
  #                           fill = "orange", color = "black", linewidth = 0.5)
  #   ))
  # }
  # Center point
  if (show.center) {
    layers <- c(layers, list(
      ggplot2::geom_point(data = data.frame(x = center[1], y = center[2]),
                          ggplot2::aes(x = x, y = y),
                          color = "red", shape = 8, size = 2)
    ))
  }

  # Return with whiskers first to guarantee they are drawn underneath
  c(whisker_layers, layers)
}

# -----------------------------
# ggplot2-based plot.bagWhiskerPlot
# -----------------------------
plot.bagWhiskerPlot <- function(
    x,
    show.outlier = TRUE,
    show.whiskers = TRUE,
    show.looppoints = TRUE,
    show.bagpoints = TRUE,
    show.loophull = TRUE,
    show.baghull = TRUE,
  show.fence_mag_bag = TRUE,
    add = FALSE,
    pch = 16, cex = .4,
    verbose = FALSE,
    col.loophull = "#aaccff",
    col.looppoints = "#3355ff",
    col.baghull = "#7799ff",
    col.bagpoints = "#000088",
  col.fence_mag_bag = "#CC33CC",
    transparency = FALSE,
    show.center = TRUE,
    whisker.fade = TRUE,
    whisker.n = 10,
    whisker.alpha.start = 0.4,
    whisker.alpha.end = 0.0,
    whisker.end.prop = 0.7,
    main = NULL,
    ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot.bagplot (ggplot2): package 'ggplot2' is required.")
  }

  # Unpack bagplot object to access xydata for scales if needed
  center <- hull.center <- hull.bag <- hull.loop <- pxy.bag <- pxy.outer <- pxy.outlier <- NULL
  hdepths <- is.one.dim <- prdata <- xy <- xydata <- exp.dk <- exp.dk.1 <- hdepth <- NULL
  tphdepth <- tp <- NULL
  bagplotobj <- x
  for (i in seq(along = bagplotobj)) {
    eval(parse(text = paste(names(bagplotobj)[i], "<-bagplotobj[[", i, "]]")))
  }

  # Build layers according to original logic
  layers <- bp_build_layers(
    x = x,
    show.outlier = show.outlier,
    show.whiskers = show.whiskers,
    show.looppoints = show.looppoints,
    show.bagpoints = show.bagpoints,
    show.loophull = show.loophull,
    show.baghull = show.baghull,
    show.fence_mag_bag = show.fence_mag_bag,
    pch = pch, cex = cex,
    col.loophull = col.loophull,
    col.looppoints = col.looppoints,
    col.baghull = col.baghull,
    col.bagpoints = col.bagpoints,
    col.fence_mag_bag = col.fence_mag_bag,
    transparency = transparency,
    show.center = show.center,
    whisker.fade = whisker.fade,
    whisker.n = whisker.n,
    whisker.alpha.start = whisker.alpha.start,
    whisker.alpha.end = whisker.alpha.end,
    whisker.end.prop = whisker.end.prop
  )
  
  # Split whisker layers and others so whiskers can be forced to the very bottom
  is_whisker <- function(ly) identical(attr(ly, "bp_role"), "whisker")
  whisker_layers <- layers[vapply(layers, is_whisker, logical(1))]
  other_layers <- layers[!vapply(layers, is_whisker, logical(1))]
  
  # Utility to support %||% without rlang
  if (!exists("%||%", mode = "function")) {
    `%||%` <- function(a, b) if (!is.null(a)) a else b
  }
  dots <- list(...)
  base_gg <- dots$gg %||% dots$base_gg %||% NULL

  # Prepare base ggplot
  # If one-dimensional, we still use xydata as background extent
  df_all <- if (!is.null(xydata)) data.frame(x = xydata[, 1], y = xydata[, 2]) else data.frame(x = numeric(0), y = numeric(0))
  if (isTRUE(add) && inherits(base_gg, "ggplot")) {
    # Prepend whiskers to existing layers so they are truly at the bottom
    p <- base_gg
    p$layers <- c(whisker_layers, p$layers, other_layers)
  } else {
    p <- ggplot2::ggplot(df_all, ggplot2::aes(x = x, y = y)) +
      whisker_layers +
      other_layers +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.x = ggplot2::element_text(size = 16),
        axis.title.y = ggplot2::element_text(size = 16, angle = 90),
        plot.title = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16),
        strip.text = ggplot2::element_text(size = 16)
      )
  }

  p <- p + ggplot2::scale_alpha_identity(guide = "none")

  # Optional xlim/ylim from ... (if provided)
  if (!is.null(dots$xlim)) {
    p <- p + ggplot2::scale_x_continuous(
      limits = dots$xlim,
      n.breaks = 6,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    )
  }
  if (!is.null(dots$ylim)) {
    p <- p + ggplot2::scale_y_continuous(
      limits = dots$ylim,
      n.breaks = 6,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    )
  }
  # If user did not specify explicit limits, still ensure non-overlapping tick labels
  if (is.null(dots$xlim)) {
    p <- p + ggplot2::scale_x_continuous(
      n.breaks = 6,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    )
  }
  if (is.null(dots$ylim)) {
    p <- p + ggplot2::scale_y_continuous(
      n.breaks = 6,
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    )
  }

  # if main title provided
  if (!is.null(main)) {
    p <- p + ggplot2::ggtitle(main)
  }

  if (!isTRUE(add) || !inherits(base_gg, "ggplot")) {
    print(p)
  }

  invisible(p)
}
