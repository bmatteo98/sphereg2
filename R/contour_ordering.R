# Ordering of unordered contour pixels into a curve
#
# These helpers are used by the two applications (see
# `applications/speech inversion/data/speaker1_at1/regenerate_vt_shapes.R` and
# `applications/speech inversion/01_extract_tongue_shapes.R`) to turn the set
# of pixels labelled as an articulator in a segmented MRI frame into an ordered
# polyline that can be smoothed and converted to an SRVF: `order.contour()`
# orders the pixels of one articulator (minimum spanning tree + diameter path),
# `orient_segments()` and `bind_with_connector()` join several articulators
# into one curve, `find.contours()` and `load.frames()` read the frames.
#
# None of these functions is exported: use `sphereg2:::order.contour()` or
# `devtools::load_all()`.

## ---------------------------------------------------------------------------
## Current implementation
## ---------------------------------------------------------------------------

#' Order the pixels of a contour along the shape
#'
#' Given the coordinates of the pixels of a thin, open, possibly noisy contour
#' (in arbitrary order), return them ordered from one end of the contour to the
#' other, so that `lines()` draws the shape instead of a tangle of segments.
#'
#' The points are treated as the vertices of a complete graph with Euclidean
#' edge weights; the minimum spanning tree (MST) of this graph is computed with
#' [igraph::mst()] and its diameter, i.e. the longest (weighted) simple path in
#' the tree, is returned by [igraph::get_diameter()]. Because the MST always
#' connects all points, gaps of any size in the pixel contour are bridged
#' without a distance tolerance to tune; stray pixels, small side branches and
#' places where the contour is more than one pixel thick become short spurs of
#' the tree that the diameter path does not visit. The diameter is a simple path
#' from one end of the curve to the other and hence already an ordering.
#'
#' Points not on the returned path are dropped (their indices are returned in
#' `$dropped`). Isolated pixels, whose nearest neighbour is further away than
#' `iso.tol`, are removed before building the tree so that they cannot extend
#' the path at its ends.
#'
#' Tracing the contour by repeatedly jumping to the closest point, the
#' approach used in an earlier version of the data preparation, frequently
#' returned only part of the shape when the next point was more than a couple
#' of pixels away; the spanning tree bridges such gaps.
#'
#' @param pts numeric matrix (or data frame) with one row per point and the
#'   x and y coordinates in the first two columns, in any order.
#' @param iso.tol points whose nearest neighbour is further than `iso.tol` are
#'   dropped beforehand as isolated noise; use `Inf` to keep every point.
#' @param min.frac if the ordered path contains less than `min.frac` times the
#'   number of (non-isolated) input points, the contour is flagged by setting
#'   `tr = tractlab` in the output.
#' @param tractlab label returned in `$tr` when the contour is flagged
#'   (the label of the articulator, in the applications).
#' @param orient logical; if `TRUE` the path is oriented so that its first point
#'   has the smaller x coordinate (same convention as [orient_segments()]).
#' @param plo logical; if `TRUE` the points and the ordered path are plotted.
#'
#' @return A list with components
#' \describe{
#'   \item{`cont`}{matrix of the ordered points (a subset of the rows of
#'     `pts`, same columns).}
#'   \item{`tr`}{`tractlab` if the contour is flagged, `0` otherwise.}
#'   \item{`frac`}{fraction of the input points that lie on the path.}
#'   \item{`idx`}{row indices of `pts` used, in path order.}
#'   \item{`dropped`}{row indices of `pts` not used.}
#' }
#'
#' @details Requires the suggested package \pkg{igraph}. The complete graph has
#'   `nrow(pts)^2` edges, which is immaterial for the few dozen to few hundred
#'   pixels of a segmented articulator but would become slow for many thousands
#'   of points.
#'
#' @seealso [orient_segments()], [bind_with_connector()] for combining several
#'   ordered contours; [find.contours()] and [load.frames()] for reading the
#'   segmented frames.
#'
#' @examples
#' \dontrun{
#' theta <- sort(runif(60, 0, 3))
#' pts <- cbind(cos(theta), sin(theta))[sample(60), ]
#' plot(pts); lines(pts)                       # unordered: a mess
#' res <- sphereg2:::order.contour(pts)
#' lines(res$cont, col = "red")                # ordered
#' }
#' @keywords internal
order.contour <- function(pts, iso.tol = 3, min.frac = 0.8, tractlab = NaN,
                          orient = TRUE, plo = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("order.contour() needs the 'igraph' package")
  pts <- as.matrix(pts)
  n <- nrow(pts)
  if (is.null(n) || n <= 2)
    return(list(cont = pts, tr = 0, frac = 1, idx = seq_len(max(n, 0)), dropped = integer(0)))
  storage.mode(pts) <- "double"
  D <- as.matrix(dist(pts[, 1:2]))

  # 1. drop isolated pixels (nearest neighbour further than iso.tol)
  nn <- apply(D + diag(Inf, n), 1, min)
  keep <- which(nn <= iso.tol)
  if (length(keep) <= 2) keep <- seq_len(n)
  Dk <- D[keep, keep, drop = FALSE]

  # 2. MST of the complete graph, then its (weighted) diameter path
  g <- igraph::graph_from_adjacency_matrix(Dk, mode = "undirected", weighted = TRUE, diag = FALSE)
  tr_g <- igraph::mst(g, weights = igraph::E(g)$weight)
  path <- as.integer(igraph::get_diameter(tr_g, weights = igraph::E(tr_g)$weight))
  idx <- keep[path]

  # 3. orientation: start from the end with the smaller x
  if (orient && pts[idx[1], 1] > pts[idx[length(idx)], 1]) idx <- rev(idx)

  cont <- pts[idx, , drop = FALSE]
  frac <- length(idx) / n
  tr <- if (frac < min.frac) tractlab else 0

  if (plo) {
    plot(pts, pch = 16, col = "grey60", asp = 1,
         main = sprintf("%d/%d pts used, tr = %s", length(idx), n, tr))
    lines(cont, col = "red")
    points(cont[1, 1], cont[1, 2], col = "blue", pch = 15)
  }
  list(cont = cont, tr = tr, frac = frac, idx = idx, dropped = setdiff(seq_len(n), idx))
}

#' Order the pixels of one articulator and keep only the contour and its flag
#'
#' Thin wrapper around [order.contour()] used by [find.contours()].
#'
#' @inheritParams order.contour
#' @param ... further arguments passed to [order.contour()] (`iso.tol`,
#'   `min.frac`, `orient`).
#' @return A list with components `cont` (ordered points) and `tr` (`tractlab`
#'   if the contour is flagged, `0` otherwise).
#' @seealso [order.contour()]
#' @keywords internal
draw.cont2 <- function(pts, plo = FALSE, tractlab = NaN, ...) {
  res <- order.contour(pts, tractlab = tractlab, plo = plo, ...)
  list(cont = res$cont, tr = res$tr)
}

## ---------------------------------------------------------------------------
## Combining several ordered contours (vocal tract: several articulators)
## ---------------------------------------------------------------------------

#' Squared Euclidean distance between two points
#'
#' @param p,q numeric vectors of the same length.
#' @return `sum((p - q)^2)`.
#' @keywords internal
dist2 <- function(p, q) sum((p - q)^2)

#' Orient a list of ordered contours consistently
#'
#' Given a list of ordered contours (matrices with x and y in the first two
#' columns), reverse them where needed so that the first one runs from smaller
#' to larger x and each following one starts close to where the previous one
#' ends. This is used to chain the contours of several articulators into a
#' single vocal-tract curve.
#'
#' @param segments list of matrices / data frames, each with the coordinates of
#'   one ordered contour in its first two columns. `NULL` elements are skipped.
#' @return The same list, with some elements reversed row-wise.
#' @seealso [bind_with_connector()]
#' @keywords internal
orient_segments <- function(segments) {
  n <- length(segments)
  segs <- segments

  if (n == 0) return(segs)

  # orient the 1st segment:
  # first point should have smaller x than last point
  first_seg <- segs[[1]]
  if (!is.null(first_seg)) {
    if (!is.matrix(first_seg) && !is.data.frame(first_seg)) {
      first_seg <- as.matrix(first_seg)
    }

    if (ncol(first_seg) >= 2 && nrow(first_seg) >= 1) {
      first_xy <- first_seg[, 1:2, drop = FALSE]
      x_first <- as.numeric(first_xy[1, 1])
      x_last  <- as.numeric(first_xy[nrow(first_xy), 1])

      if (x_first > x_last) {
        segs[[1]] <- segs[[1]][nrow(first_xy):1, , drop = FALSE]
      }
    }
  }

  # orient the following segments relative to the previous one
  for (k in seq_len(n - 1)) {
    cur <- segs[[k]]
    nxt <- segs[[k + 1]]

    if (is.null(cur) || is.null(nxt)) next

    if (!is.matrix(cur) && !is.data.frame(cur)) {
      cur <- as.matrix(cur)
    }
    if (!is.matrix(nxt) && !is.data.frame(nxt)) {
      nxt <- as.matrix(nxt)
    }

    if (ncol(cur) < 2 || ncol(nxt) < 2) next

    cur_xy <- cur[, 1:2, drop = FALSE]
    nxt_xy <- nxt[, 1:2, drop = FALSE]

    if (nrow(cur_xy) < 1 || nrow(nxt_xy) < 1) next

    p_end   <- as.numeric(cur_xy[nrow(cur_xy), ])
    q_start <- as.numeric(nxt_xy[1, ])
    q_end   <- as.numeric(nxt_xy[nrow(nxt_xy), ])

    d_to_start <- dist2(p_end, q_start)
    d_to_end   <- dist2(p_end, q_end)

    if (d_to_end < d_to_start) {
      segs[[k + 1]] <- segs[[k + 1]][nrow(nxt_xy):1, , drop = FALSE]
    }
  }

  segs
}

#' Check that every part of a frame has enough points
#'
#' @param frame_parts list of contours (matrices / data frames), one per
#'   articulator; `NULL` elements count as 0 points.
#' @param min_points numeric vector of required minimum numbers of points, of
#'   the same length as `frame_parts`.
#' @return `TRUE` if every part has at least the required number of rows,
#'   `FALSE` otherwise (also when the inputs are not a list or their lengths
#'   differ).
#' @keywords internal
is_well_recovered <- function(frame_parts, min_points) {
  # frame_parts: list of parts (data frames / matrices)
  # min_points: vector of required minima, same length

  if (!is.list(frame_parts)) return(FALSE)
  if (length(frame_parts) != length(min_points)) return(FALSE)

  # lengths of each part
  frame_lengths <- sapply(frame_parts, function(p) {
    if (is.null(p)) return(0L)
    nrow(p)
  })

  all(frame_lengths >= min_points)
}

#' Straight-line connector between two points
#'
#' @param p1,p2 numeric vectors `c(x, y)`.
#' @param n number of equally spaced points on the segment (end points
#'   included).
#' @return An `n x 2` matrix with columns `x` and `y`.
#' @keywords internal
make_connector <- function(p1, p2, n = 10) {
  # p1, p2 are numeric vectors of length 2: c(x, y)
  cbind(
    x = seq(p1[1], p2[1], length.out = n),
    y = seq(p1[2], p2[2], length.out = n)
  )
}

#' Bind oriented contours into one curve, bridging parts 3 and 4
#'
#' Row-binds a list of consistently oriented contours (see
#' [orient_segments()]) into a single matrix, inserting a straight-line
#' connector of `n_conn` points between the end of the third and the start of
#' the fourth part. This is specific to the vocal-tract application, where the
#' third and fourth articulators are not adjacent in the image.
#'
#' @param segments_fixed list of at least four matrices / data frames with the
#'   coordinates in their first two columns, already oriented.
#' @param n_conn number of points of the connector, see [make_connector()].
#' @return A matrix with two columns containing parts 1-3, the connector and
#'   parts 4 onwards.
#' @keywords internal
bind_with_connector <- function(segments_fixed, n_conn = 10) {
  # keep only first two columns if segments are matrices/data frames
  segs <- lapply(segments_fixed, function(s) as.matrix(s[, 1:2, drop = FALSE]))

  # last point of part 3
  p3_end <- as.numeric(segs[[3]][nrow(segs[[3]]), ])

  # first point of part 4
  p4_start <- as.numeric(segs[[4]][1, ])

  # interpolated connector
  conn34 <- make_connector(p3_end, p4_start, n = n_conn)

  # bind in the desired order:
  # parts 1,2,3, connector, 4,5,...
  out <- do.call(rbind, c(segs[1:3], list(conn34), segs[4:length(segs)]))
  out
}

## ---------------------------------------------------------------------------
## Extracting the contours from segmented MRI frames
## ---------------------------------------------------------------------------

#' Extract and order the articulator contours of one segmented frame
#'
#' For each label in `labels`, collects the pixels of `im` carrying that label
#' and orders them along the contour with [draw.cont2()].
#'
#' `im` is expected in the orientation produced by [load.frames()], i.e. an
#' `nx x ny` matrix (or an `nx x ny x 1 x 1` array as used with \pkg{imager})
#' whose first index is the horizontal pixel coordinate x and whose second
#' index is the vertical coordinate y, y increasing upwards. Pixel coordinates
#' are then `which(im == i, arr.ind = TRUE)`, which coincides with
#' `imager::where(im == i)`.
#'
#' @param im integer matrix / array of segmentation labels, see Details.
#' @param labels integer vector of the labels to extract (default `2`, the
#'   tongue in the speech-inversion data).
#' @param plt logical; if `TRUE` the ordered contours are drawn (one colour per
#'   label).
#' @param ... further arguments passed to [draw.cont2()] / [order.contour()]
#'   (`iso.tol`, `min.frac`, `orient`).
#' @return A list with components `beta` (list of ordered contours, one
#'   matrix per label, in the order of `labels`) and `labtr` (vector of the
#'   `tr` flags returned by [draw.cont2()], i.e. the label if that contour is
#'   suspicious and `0` otherwise).
#' @seealso [load.frames()], [order.contour()]
#' @keywords internal
find.contours <- function(im, labels = 2, plt = FALSE, ...) {
  if (length(dim(im)) > 2) im <- im[, , 1, 1]
  beta <- list()
  labtr <- c()
  mycol <- grDevices::colorRampPalette(c("green", "blue", "red"))(max(labels))
  if (plt) plot(-10:85, ty = "n", xlab = "", ylab = "")
  for (i in labels) {
    label <- which(im == i, arr.ind = TRUE)
    colnames(label) <- c("x", "y")
    seg <- draw.cont2(label, plo = FALSE, tractlab = i, ...)
    if (plt) lines(seg$cont, col = mycol[i], pch = 16)
    beta <- append(beta, list(seg$cont))
    labtr <- c(labtr, seg$tr)
  }
  list(beta = beta, labtr = labtr)
}

#' Load a folder of segmented MRI frames and extract their contours
#'
#' Reads every `.npy` file of `path1` (one segmented frame each, as produced by
#' the neural-network contour extraction), flips it to the x/y orientation
#' expected by [find.contours()] and extracts the ordered contours of the
#' requested labels.
#'
#' The `.npy` files are read with \pkg{reticulate} through NumPy
#' (`numpy.load`), which must therefore be available in the Python
#' configuration used by \pkg{reticulate}.
#'
#' @param path1 path to the folder containing the `.npy` frames.
#' @param labels labels to extract, passed to [find.contours()].
#' @param np optional NumPy module as returned by
#'   `reticulate::import("numpy")`; imported if `NULL`.
#' @param ... further arguments passed to [find.contours()].
#' @return A list with components `video` (list with one element per frame,
#'   each the `beta` list of [find.contours()]) and `trlab` (matrix with one row
#'   per frame and one column per label containing the `tr` flags).
#' @seealso [find.contours()]
#' @keywords internal
load.frames <- function(path1, labels = 2, np = NULL, ...) {
  if (is.null(np)) {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("load.frames() needs the 'reticulate' package (or pass `np`)")
    np <- reticulate::import("numpy")
  }
  trlab <- c()
  label.pixels <- list.files(path1, pattern = "\\.npy$")
  video <- list()
  for (imag in label.pixels) {
    curr_frame <- np$load(file.path(path1, imag))
    # rows of the .npy are image rows from the top: flip so that the first
    # index is x and the second is y (increasing upwards)
    im <- t(curr_frame[nrow(curr_frame):1, ])
    cont <- find.contours(im, labels = labels, ...)
    video <- append(video, list(cont$beta))
    trlab <- rbind(trlab, cont$labtr)
  }
  list(video = video, trlab = trlab)
}
