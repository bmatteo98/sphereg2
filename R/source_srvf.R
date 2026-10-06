# SRVF functions

#' Cumulative trapezoidal integration of a single vector
#'
#' Drop-in replacement for `pracma::cumtrapz(x, y)` restricted to vector
#' input, computing the identical values without the `repmat`/`kronecker`
#' overhead of the general matrix version. Called once per iteration in
#' the innermost alignment loops, where that overhead dominated runtime.
#'
#' @param y numeric vector of function values.
#' @param x numeric vector of evaluation points, same length as `y`; may be
#'   omitted when `dxh` is supplied.
#' @param dxh optional precomputed half-steps `diff(x) / 2`, so callers
#'   integrating repeatedly on the same grid can skip recomputing them.
#' @return An `length(y) x 1` matrix, as returned by [pracma::cumtrapz()].
#' @noRd
cumtrapz_vec <- function(y, x = NULL, dxh = NULL) {
  m <- length(y)
  if (is.null(dxh)) dxh <- diff(x) / 2
  ct <- c(0, cumsum(dxh * (y[1:(m - 1)] + y[2:m])))
  dim(ct) <- c(m, 1L)
  ct
}

# cache for the trapezoidal weights diff(seq(0, 1, length.out = n)) used by
# L2prod on the default equispaced grid; a handful of grid lengths recur, so
# they are kept in parallel vectors and found with match() (a character-keyed
# environment turned out measurably slower in this hot path)
.L2prod_dx_cache <- new.env(parent = emptyenv())
.L2prod_dx_cache$ns <- integer(0)
.L2prod_dx_cache$dxs <- list()

# curves and their square-root velocity representation ---------------------

#' Smooth and resample a planar curve
#'
#' Apply a univariate smoothing spline to each coordinate of a discretised
#' planar curve and re-evaluate it on a regular grid of `lo` points. This
#' is the standard pre-processing step applied to hand-drawn or
#' landmark-digitised outlines before passing them to [fdasrvf::curve_to_q()]: raw
#' digitisations are typically noisy and irregularly sampled, whereas the
#' SRVF transform involves a numerical derivative and hence amplifies
#' high-frequency noise considerably.
#'
#' Both coordinates are smoothed *independently as functions of their own
#' index*, so the operation is a smoothing of the parametrised curve
#' \eqn{t \mapsto \beta(t)} rather than of its image. The orientation is
#' normalised at the end: if the smoothed \eqn{x}-coordinate is decreasing
#' on average (negative slope of a simple linear fit against the index),
#' both coordinates are reversed, so that curves traversed in opposite
#' directions become comparable.
#'
#' @param betaj matrix (or object coercible to one) with two columns
#'   holding the \eqn{x}- and \eqn{y}-coordinates of the curve.
#' @param lo number of points of the output grid; defaults to `1000`.
#' @param spar smoothing parameter in \eqn{[0, 1]} forwarded to
#'   [stats::smooth.spline()]. Larger values give smoother curves.
#'
#' @return A `lo x 2` matrix with the smoothed curve, possibly reversed so
#'   that the \eqn{x}-coordinate increases on average.
#'
#' @seealso [fdasrvf::curve_to_q()] for the subsequent SRVF transform.
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 60)
#' beta <- cbind(t, sin(2 * pi * t)) + matrix(rnorm(120, sd = .02), ncol = 2)
#' smooth_beta <- smooth.curve(beta, lo = 200)
smooth.curve <- function (betaj, lo= 1000, spar=0.3){
  betaj = as.matrix(betaj)
  b.sx = smooth.spline(as.numeric(betaj[,1]), spar=spar)
  b.sy = smooth.spline(as.numeric(betaj[,2]), spar=spar)
  pr.b.sx = predict(b.sx, c(seq(range(b.sx$x)[1],range(b.sx$x)[2] ,length.out=lo)))$y
  pr.b.sy = predict(b.sy, c(seq(range(b.sy$x)[1],range(b.sy$x)[2] ,length.out=lo)))$y
  if (lm(pr.b.sx~c(1:length(pr.b.sx)))$coef[2]<0) return(cbind(rev(pr.b.sx),rev(pr.b.sy)))
  return(cbind(pr.b.sx,pr.b.sy))
}

# rotation alignment -------------------------------------------------------

#' Optimal rotation between two square-root velocity functions
#'
#' Solve the Procrustes problem of rotating one SRVF onto another: find
#' \deqn{O^\ast = \arg\min_{O \in SO(2)} \|q_1 - O q_2\|_{L^2}
#'   = \arg\max_{O \in SO(2)} \mathrm{tr}\bigl(O^\top A\bigr),
#'   \qquad A = \int_0^1 q_1(t)\, q_2(t)^\top\, dt.}
#' Writing \eqn{A = U \Sigma V^\top} for the singular value
#' decomposition, the maximiser is \eqn{O^\ast = U V^\top}; when
#' \eqn{\det A < 0} the second singular direction is flipped
#' (\eqn{O^\ast = U \, \mathrm{diag}(1, -1)\, V^\top}) so that the
#' result is a rotation rather than a reflection, i.e. shapes are *not*
#' identified with their mirror images.
#'
#' The rotation group is one of the shape-preserving transformations that
#' the elastic distance quotients out; the other is reparametrisation,
#' handled by [opt.reparam.curve()]. [align_curves()] applies both.
#'
#' @param q1 matrix with two columns: the reference SRVF.
#' @param q2 matrix with two columns of the same dimension as `q1`: the
#'   SRVF to be rotated.
#'
#' @return A `2 x 2` rotation matrix in \eqn{SO(2)}.
#'
#' @seealso [align_curves()], [opt.reparam.curve()].
#' @export
opt.O <- function (q1, q2){
  I = 0
  for ( i in 1: nrow(q1)){
    I = I + (q1[i,])%*% t(q2[i,])/ (nrow(q1)-1)
  }
  C = I
  svdC= svd(C)
  if (det(C)>=0){
    O = svdC$u%*%t(svdC$v)
  }
  if (det(C)<0){
    H = diag(1, 2,2)
    H[2,2]=-1
    O = svdC$u%*%H%*%t(svdC$v)
  }
  return(O)
}

#' Align one SRVF curve to another over rotations and reparametrisations
#'
#' Bring `q2` into optimal correspondence with the reference `q1` by
#' composing the two group actions the elastic shape distance quotients
#' out: reparametrisation and rotation. Concretely, the routine
#' \enumerate{
#'   \item computes the optimal warping of `q2` towards `q1` with
#'     [opt.reparam.curve()] and applies it,
#'   \item computes the optimal rotation of the warped curve with
#'     [opt.O()] and applies it,
#'   \item renormalises the result to unit norm with respect to `S`.
#' }
#' The renormalisation matters because the discretised warping action is
#' only approximately norm-preserving, whereas the downstream spherical
#' machinery ([sphere_regression()], [karcher_mean()]) requires its
#' inputs to lie exactly on the unit sphere.
#'
#' The distance between `q1` and the returned `aligned_q2`, measured on
#' the sphere, is the (discretised) elastic shape distance between the
#' two curves. This is exactly how [sphere_regression()] with
#' `elastic = TRUE` re-aligns the responses to the current fit at every
#' iteration.
#'
#' @param q1 matrix with two columns: the reference SRVF, assumed to be
#'   of unit norm.
#' @param q2 matrix with two columns of the same dimension as `q1`: the
#'   SRVF to be aligned.
#' @param S a sphere geometry as returned by [initialize_sphere()],
#'   supplying the norm used for the final renormalisation.
#' @param eps,maxit convergence threshold and maximum number of iterations
#'   of the warping optimisation; only used when `optimizer = "sd"`.
#' @param optimizer solver for the warping, passed on to [opt.gamj()]:
#'   `"lbfgs"` (the default), `"sd"` or `"rlbfgs"`.
#'
#' @return A list with elements
#'   \describe{
#'     \item{`aligned_q2`}{Matrix of the same dimension as `q2`: the
#'       rotated, reparametrised and renormalised SRVF.}
#'     \item{`opt_gam`}{The warping as returned by [opt.gamj()], a list
#'       with the warping `path` and its derivative `dpath`.}
#'   }
#'
#' @seealso [align_vectorized_curves()] for the vectorised interface used
#'   throughout the package, [align_multiple_curves()] for aligning a
#'   whole sample, [opt.O()], [opt.reparam.curve()].
#' @export
align_curves = function (q1, q2, S, eps = 0.05, maxit = 50,
                         optimizer = c("lbfgs", "sd", "rlbfgs")){

  opt_rep = opt.reparam.curve(q1, q2, eps = eps, maxit = maxit,
                              optimizer = match.arg(optimizer))
  q2_ = opt_rep$rep_q2
  O = opt.O(q1, q2_)
  q2_ =  t(O %*% t(q2_))
  aligned_q2 = q2_/S$norm.u(q2_)
  return (list(aligned_q2 = aligned_q2 , opt_gam = opt_rep$otp_gam))
}

#' Align two SRVF curves given in stacked-vector form
#'
#' Vectorised front end to [align_curves()]. Throughout the regression
#' code a planar SRVF sampled at \eqn{m} points is stored as a single
#' numeric vector of length \eqn{2m},
#' \eqn{v = (q_{1,x}, \ldots, q_{m,x}, q_{1,y}, \ldots, q_{m,y})}, so
#' that a sample of curves is an ordinary matrix with one curve per
#' column and the generic sphere routines apply unchanged. This function
#' unstacks the two inputs, aligns them, and returns the result in the
#' same stacked form.
#'
#' The inner product of `S` — which is defined on the stacked
#' representation — is lifted to matrix-valued arguments internally, so
#' that the renormalisation inside [align_curves()] uses the geometry the
#' caller intended. Note that the low-dimensional coordinate
#' representation used by [sphere_regression()] (see [subspace_setup()])
#' is not a valid input here: alignment has to be carried out in the
#' curves' native representation, which is why the elastic fitting loop
#' transforms back and forth via `iso$trafo()`.
#'
#' @param v1 numeric vector of length \eqn{2m}: the stacked reference SRVF.
#' @param v2 numeric vector of length \eqn{2m}: the stacked SRVF to be aligned.
#' @param S a sphere geometry as returned by [initialize_sphere()], whose
#'   inner product refers to the stacked representation. Defaults to the
#'   Euclidean geometry.
#' @param eps,maxit convergence threshold and maximum number of iterations
#'   of the warping optimisation; only used when `optimizer = "sd"`.
#' @param optimizer solver for the warping, passed on to [opt.gamj()]:
#'   `"lbfgs"` (the default), `"sd"` or `"rlbfgs"`.
#'
#' @return A list with elements `aligned_q2` (numeric vector of length
#'   \eqn{2m}) and `opt_gam`, as described in [align_curves()].
#'
#' @seealso [align_curves()], [align_multiple_curves()],
#'   [sphere_regression()].
#'
#' @references
#' Huang, W., Gallivan, K. A., Srivastava, A., and Absil, P.-A. (2015).
#' Riemannian optimization for registration of curves in elastic shape analysis.
#'
#' @export
align_vectorized_curves <- function(v1, v2, S = initialize_sphere(), eps = 0.05, maxit = 50,
                                    optimizer = c("lbfgs", "sd", "rlbfgs")) {
  mm <- length(v1) / 2
  q1 <- cbind(v1[1:mm], v1[(mm + 1):(2 * mm)])
  q2 <- cbind(v2[1:mm], v2[(mm + 1):(2 * mm)])

  pprod <- function(q1, q2) S$prod.u(c(q1), c(q2))
  S_intrn <- initialize_sphere(prod.u = pprod)

  return(align_curves(q1, q2, S_intrn, eps = eps, maxit = maxit,
                      optimizer = match.arg(optimizer)))
}

#' Trapezoidal inner product and norm on \eqn{L^2([0,1])}
#'
#' `L2prod()` approximates \eqn{\langle v, l\rangle_{L^2} =
#' \int v(t)\, l(t)\, dt} by the trapezoidal rule from function values
#' given on a (by default equispaced) grid; `L2norm()` is the induced norm.
#'
#' These are the primitives of the Hilbert-sphere geometry in which the
#' warping functions live: a warping \eqn{\gamma} is represented by
#' \eqn{\psi = \sqrt{\gamma'}}, which satisfies \eqn{\|\psi\|_{L^2} = 1}
#' because \eqn{\int_0^1 \gamma'(t)\, dt = \gamma(1) - \gamma(0) = 1}.
#' Optimisation over warpings therefore becomes optimisation over the
#' unit sphere of \eqn{L^2([0,1])}, with [proj_tg_l()] and [exp_map_l()]
#' as the corresponding projection and exponential map.
#'
#' Unlike [trap()], which assumes a regular grid, `L2prod()` accepts an
#' arbitrary increasing `grid`.
#'
#' @param v numeric vector of function values.
#' @param l numeric vector of function values, of the same length as `v`.
#' @param grid optional numeric vector of evaluation points in
#'   \eqn{[0,1]}; if `NULL` (the default), an equispaced grid is assumed.
#'
#' @return Numeric scalar: the approximate inner product, respectively norm.
#'
#' @seealso [trap()], [proj_tg_l()], [exp_map_l()].
#' @keywords internal
L2prod <- function(v, l, grid = NULL) {
  if (length(v) != length(l)) stop("length(v) == length(l) is not TRUE")

  if (is.null(grid)) {
    # equispaced default grid: reuse the cached weights diff(seq(0, 1, ...))
    # instead of rebuilding grid and differences on every call
    n <- length(v)
    i <- match(n, .L2prod_dx_cache$ns)
    if (is.na(i)) {
      .L2prod_dx_cache$ns <- c(.L2prod_dx_cache$ns, n)
      i <- length(.L2prod_dx_cache$ns)
      .L2prod_dx_cache$dxs[[i]] <- diff(seq(0, 1, length.out = n))
    }
    dx <- .L2prod_dx_cache$dxs[[i]]
  } else {
    dx <- diff(grid)
  }

  f <- v * l
  integral <- sum((f[-1] + f[-length(f)]) / 2 * dx)

  return(integral)
}

#' @rdname L2prod
#' @keywords internal
L2norm <- function(v) sqrt(L2prod(v, v))

#' Tangent projection and exponential map on the sphere of warpings
#'
#' Geometric primitives for the unit sphere of \eqn{L^2([0,1])}, on which
#' the square-root representations \eqn{\psi = \sqrt{\gamma'}} of the
#' warping functions live. They mirror the sphere operations used in the
#' elastic registration literature.
#'
#' `proj_tg_l()` removes from `v` its component along `l`,
#' \deqn{P_{T_l}(v) = v - \frac{\langle v, l\rangle}{\langle l, l\rangle}\, l,}
#' which is needed because the Euclidean gradient of the alignment cost
#' does not lie in the tangent space of the manifold.
#' `exp_map_l()` moves from `l` along the geodesic in direction `v`,
#' \deqn{\mathrm{Exp}_l(v) = \cos(\|v\|)\, l + \sin(\|v\|)\, v/\|v\|,}
#' returning a point on the sphere again.
#'
#' @param l numeric vector: a point on the unit sphere of
#'   \eqn{L^2([0,1])}.
#' @param v numeric vector of the same length as `l`.
#'
#' @return Numeric vector of the same length as `l`.
#'
#' @seealso [L2prod()], [steep_desc()], [initialize_sphere()].
#' @keywords internal
proj_tg_l <- function(l, v) {
  normaliz <- L2prod(v, l) / L2prod(l, l)
  lnor <- l * normaliz
  return(v - lnor)
}

#' @rdname L2prod
#' @keywords internal
exp_map_l <- function(l, v) {
  n_v <- L2norm(v)
  return(cos(n_v) * l + v * (sin(n_v) / n_v))
}

#' Gradient of the elastic alignment cost with respect to the warping
#'
#' Riemannian gradient of the alignment criterion
#' \deqn{E(\gamma) = \bigl\| q_1 O -
#'   (q_2 \circ \gamma)\sqrt{\gamma'} \bigr\|_{L^2}^2 + w\, B(\psi)}
#' with respect to the square-root representation
#' \eqn{\psi = \sqrt{\gamma'}} of the warping, projected onto the tangent
#' space \eqn{T_\psi} by [proj_tg_l()].
#'
#' The composition \eqn{q_2 \circ \gamma} and its derivative are obtained
#' by evaluating natural cubic spline interpolants of the two coordinates
#' of `q2` at `gam`; `gam` itself is interpolated by a monotonicity-
#' preserving spline so that the numerical derivative \eqn{\gamma'}
#' remains non-negative.
#'
#' @param gam numeric vector of length \eqn{n}: the current warping
#'   \eqn{\gamma} evaluated on the regular grid of \eqn{[0,1]}.
#' @param q1 matrix with two columns and \eqn{n} rows: the reference SRVF.
#' @param q2 matrix of the same dimension as `q1`: the SRVF whose
#'   parametrisation is being matched to that of `q1`.
#' @param O `2 x 2` rotation matrix in \eqn{SO(2)} applied to `q1`.
#' @param w non-negative weight of the barrier term.
#' @param spq2x,spq2y optional natural-spline interpolants of the two
#'   coordinates of `q2` (as returned by [stats::splinefun()]); if `NULL`,
#'   they are constructed internally. Callers evaluating the gradient
#'   repeatedly for the same `q2` should construct them once and pass them
#'   in.
#'
#' @return A list with elements `g`, `l`, and `lossg`.
#'
#' @seealso [loss()], [armijo()], [steep_desc()], [opt.gamj()].
#' @references
#' Huang, W., Gallivan, K. A., Srivastava, A., and Absil, P.-A. (2015).
#' Riemannian optimization for registration of curves in elastic shape analysis.
#' @keywords internal
grad <- function(gam, q1, q2, O = diag(rep(1, 2)), w = 1e-4, spq2x = NULL, spq2y = NULL,
                 l = NULL) {
  n <- nrow(q1)
  ts <- seq(0, 1, length.out = n)
  # callers that already hold psi should pass it in: recovering l from gam by
  # differentiating a monotone spline is not the inverse of gam = int psi^2, and
  # the round trip costs 0.6-2% in relative L2 norm on typical iterates, which
  # makes the returned gradient that of a slightly different function than the
  # one loss() evaluates
  if (is.null(l)) {
    spgam <- splinefun(ts, gam, method = "monoH.FC")
    l <- sqrt(spgam(ts, deriv = 1))
  }
  # the q2 interpolants only depend on q2, so callers iterating over warpings
  # (steep_desc) construct them once and pass them in
  if (is.null(spq2x)) {
    spq2x <- splinefun(ts, q2[, 1], method = "natural")
    spq2y <- splinefun(ts, q2[, 2], method = "natural")
  }

  q2rox <- spq2x(gam); q2roy <- spq2y(gam)
  q2ro1x <- spq2x(gam, 1); q2ro1y <- spq2y(gam, 1)

  # multiplication by the default identity rotation is exact, so it is skipped
  q1L <- if (missing(O)) q1 else q1 %*% O
  q1O <- if (missing(O)) q1 else q1 %*% t(O) # rows hold (O %*% q1[t, ])^T
  # stacking the columns reproduces the column-major flattening that the
  # former matrix expression produced inside L2norm
  lossg <- L2norm(c(q1L[, 1] - q2rox * l, q1L[, 2] - q2roy * l))
  # explicit two-term dot products rather than rowSums: keeps the rounding
  # order of the original per-row BLAS products, so results are bit-identical
  x <- q1O[, 1] * q2rox + q1O[, 2] * q2roy
  y1 <- q1O[, 1] * (2 * l * q2ro1x) + q1O[, 2] * (2 * l * q2ro1y)
  z <- w * l * (2 - 1 / (l^4)) * sqrt(1 + l^4)
  # the barrier gradient diverges where gamma' = 0 (monoH.FC clamps flat
  # regions to exactly zero); drop those contributions to keep the descent finite
  z[l == 0] <- 0
  y <- cumtrapz_vec(y1, ts)
  g <- proj_tg_l(l, 2 * (y * l) - 2 * x + 2 * z)
  return(list(g = g, l = l, lossg = lossg))
}

#' Elastic alignment cost for a given warping
#'
#' Evaluate the alignment criterion
#' \deqn{\bigl\| q_1 O - (q_2 \circ \gamma)\, \psi \bigr\|_{L^2}}
#' whose square is minimised by [steep_desc()].
#'
#' @param q1 matrix with two columns: the reference SRVF.
#' @param q2 matrix of the same dimension as `q1`: the SRVF being warped.
#' @param gam numeric vector: the warping \eqn{\gamma} on the regular
#'   grid of \eqn{[0,1]}.
#' @param l numeric vector: the corresponding \eqn{\psi = \sqrt{\gamma'}}.
#' @param O `2 x 2` rotation matrix applied to `q1`.
#' @param spq2x,spq2y optional natural-spline interpolants of the two
#'   coordinates of `q2`; if `NULL`, they are constructed internally. See
#'   [grad()].
#'
#' @return Numeric scalar: the \eqn{L^2} distance between the rotated
#'   reference and the warped curve.
#'
#' @seealso [grad()], [armijo()], [steep_desc()].
#' @references
#' Huang, W., Gallivan, K. A., Srivastava, A., and Absil, P.-A. (2015).
#' Riemannian optimization for registration of curves in elastic shape analysis.
#' @keywords internal
loss <- function(q1, q2, gam, l, O = diag(rep(1, 2)), spq2x = NULL, spq2y = NULL) {
  if (is.null(spq2x)) {
    n <- length(gam)
    spq2x <- splinefun(seq(0, 1, length.out = n), q2[, 1], method = "natural")
    spq2y <- splinefun(seq(0, 1, length.out = n), q2[, 2], method = "natural")
  }
  # skip the exact multiplication by the default identity rotation and
  # reproduce the former matrix expression's column-major flattening directly
  q1L <- if (missing(O)) q1 else q1 %*% O
  return(L2norm(c(q1L[, 1] - spq2x(gam) * l, q1L[, 2] - spq2y(gam) * l)))
}

#' Armijo backtracking line search on the sphere of warpings
#'
#' Choose a step size \eqn{\alpha} for the gradient step
#' \eqn{\psi \mapsto \mathrm{Exp}_\psi(-\alpha g)} that achieves
#' sufficient decrease of the alignment cost,
#' \deqn{E(\psi) - E(\mathrm{Exp}_\psi(-\alpha g)) \;\geq\;
#'   \sigma\, \alpha\, \|g\|^2,}
#' by starting from `al` and repeatedly multiplying it by `beta` until
#' the condition holds. Each candidate step is mapped back onto the
#' sphere by [exp_map_l()] and re-integrated to a warping via
#' \eqn{\gamma(t) = \int_0^t \psi^2}, so that only admissible warpings
#' are ever evaluated.
#'
#' The loop terminates early if both sides of the condition fall below
#' `1e-20`, which happens once the iterate is essentially at a stationary
#' point and the comparison becomes numerical noise.
#'
#' @param q1 matrix with two columns: the reference SRVF.
#' @param q2 matrix of the same dimension as `q1`: the SRVF being warped.
#' @param lossg numeric scalar: the cost at the current iterate, as
#'   returned by [grad()].
#' @param loss the cost function to evaluate candidate steps with,
#'   normally [loss()]; passed as an argument so that the line search
#'   stays independent of the particular criterion.
#' @param gam numeric vector: the current warping.
#' @param g numeric vector: the current projected gradient.
#' @param al initial (largest) step size.
#' @param sigma sufficient-decrease constant \eqn{\sigma} of the Armijo
#'   condition.
#' @param beta backtracking factor in \eqn{(0,1)} by which the step size
#'   is shrunk after each failed trial.
#' @param l optional numeric vector \eqn{\psi = \sqrt{\gamma'}} matching
#'   `gam`, as returned by [grad()]; if `NULL`, it is recomputed from
#'   `gam`.
#' @param dxh optional precomputed trapezoidal half-steps of the regular
#'   grid, `diff(seq(0, 1, length.out = length(gam))) / 2`; if `NULL`,
#'   they are recomputed.
#'
#' @return Numeric scalar: the accepted step size.
#'
#' @seealso [grad()], [loss()], [steep_desc()].
#' @keywords internal
armijo <- function (q1, q2, lossg, loss, gam, g, al = 1, sigma = 1e-5, beta = 0.5, l = NULL, dxh = NULL){
  n = length(gam)
  if (is.null(dxh)) dxh = diff(seq(0,1,length.out=n)) / 2
  if (is.null(l)) {
    # identical to gradlist$l of grad(); callers that already hold it pass it in
    ts = seq(0,1,length.out=n)
    spgam = splinefun(ts, gam, method = 'monoH.FC')
    l = sqrt(spgam(ts, 1))
  }
  m = 0
  al0 = al   # backtracking is scaled from the starting step, see the note below
  newl = exp_map_l(l, - al * g)
  newg =  cumtrapz_vec(newl**2, dxh = dxh)
  sgg = sigma * t(g)%*% g
  c1 = (lossg - loss(q1, q2, newg, newl))
  c2 = al * sgg
  # a NaN in either side (degenerate warping) makes the comparison error out;
  # signal "no usable step" with a zero step size, handled by the caller
  if (!is.finite(c1) || !is.finite(c2)) return (0)
  while (c1 < c2){
    m =m+ 1
    al = al0 * beta ^ m
    newl = exp_map_l(l, - al * g)
    newg =  cumtrapz_vec(newl**2, dxh = dxh)
    c1 = (lossg - loss(q1, q2, newg, newl))
    c2 = al * sgg
    if (!is.finite(c1) || !is.finite(c2)) return (0)
    if ((c1 < 1e-20) &&(c2< 1e-20)) return (al)
  }
  return (al)
}

#' Steepest descent for the optimal warping
#'
#' Minimise the elastic alignment cost over warping functions by
#' Riemannian steepest descent on the unit sphere of \eqn{L^2([0,1])}:
#' starting from `gam`, repeatedly compute the projected gradient with
#' [grad()], pick a step size with [armijo()], move along the geodesic
#' with [exp_map_l()], and re-integrate the resulting
#' \eqn{\psi = \sqrt{\gamma'}} to a warping
#' \eqn{\gamma(t) = \int_0^t \psi^2}. Iteration stops as soon as the
#' gradient norm drops below `eps`, or after `maxit` steps.
#'
#' The return value is \eqn{\psi}, not \eqn{\gamma}: callers such as
#' [opt.gamj()] recover the warping by numerical integration, which
#' guarantees monotonicity of the result by construction.
#'
#' @param gam numeric vector: starting warping, evaluated on the regular
#'   grid of \eqn{[0,1]}. [opt.gamj()] starts from the identity
#'   \eqn{\gamma(t) = t}.
#' @param q1 matrix with two columns: the reference SRVF.
#' @param q2 matrix of the same dimension as `q1`: the SRVF being warped.
#' @param O `2 x 2` rotation matrix applied to `q1`; defaults to the
#'   identity.
#' @param w weight of the barrier term, see [grad()].
#' @param eps convergence threshold on the \eqn{L^2} norm of the
#'   projected gradient.
#' @param maxit maximum number of iterations; the current iterate is
#'   returned without warning if it is reached.
#'
#' @return Numeric vector: the final \eqn{\psi = \sqrt{\gamma'}}.
#'
#' @seealso [opt.gamj()], [grad()], [armijo()], [exp_map_l()].
#' @references
#' Huang, W., Gallivan, K. A., Srivastava, A., and Absil, P.-A. (2015).
#' Riemannian optimization for registration of curves in elastic shape analysis.
#' @keywords internal
steep_desc <- function (gam, q1, q2, O = diag(rep(1, 2)), w = 1e-4, eps = 0.05, maxit = 50){
  # l is a vector; len(l) = q1.shape[0] = q2.shape[0]; represents the reparametriz
  # grad is the gradient function
  # stop when L2norm(grad_k) < eps

  # q2 is fixed throughout the descent: construct its interpolants once and
  # share them with every grad()/loss() evaluation of the line search
  n <- nrow(q2)
  ts <- seq(0, 1, length.out = n)
  dxh <- diff(ts) / 2
  spq2x <- splinefun(ts, q2[, 1], method = "natural")
  spq2y <- splinefun(ts, q2[, 2], method = "natural")
  # forward the missing-ness of O so grad()/loss() can skip the exact
  # multiplication by the default identity rotation
  Oid <- missing(O)
  loss_ <- function(q1, q2, gam, l, O)
    if (missing(O)) loss(q1, q2, gam, l, spq2x = spq2x, spq2y = spq2y) else
      loss(q1, q2, gam, l, O, spq2x = spq2x, spq2y = spq2y)

  for (i in 1:maxit){
    gradlist = if (Oid) grad(gam, q1, q2, w = w, spq2x = spq2x, spq2y = spq2y) else
      grad(gam, q1, q2, O, w, spq2x = spq2x, spq2y = spq2y) # compute the gradient
    g = gradlist$g
    l = gradlist$l
    lossg = gradlist$lossg
    normg = L2norm(g)

    # near-degenerate warpings (gamma' touching 0) can produce NaN in the
    # gradient; no further descent is possible from there, so return the
    # current point rather than crash in the line search
    if (!is.finite(normg)) return(l)
    if (normg < eps) return(l) # stop if the norm of the gradient is close to 0

    al = armijo (q1, q2, lossg, loss_, gam, g, l = l, dxh = dxh)
    if (al <= 0) return(l) # line search found no usable step (see armijo)
    l = exp_map_l(l, - al * g) # take a step towards the negative direction of the gradient and project on the manifold
    gam = cumtrapz_vec(l^2, dxh = dxh)
  }
  #print('Out for iterations')
  return (l)
}

#' Quasi-Newton solver for the optimal warping
#'
#' Minimise the same alignment cost as [steep_desc()], but with L-BFGS-B in
#' place of the Riemannian steepest descent, reusing [grad()] unchanged.
#'
#' The unknown is \eqn{\psi = \sqrt{\gamma'}} as an *unconstrained* vector in
#' \eqn{R^n}. Dividing \eqn{\psi} by its norm inside the objective makes the
#' cost scale invariant, so the unit-sphere constraint disappears and the
#' Euclidean gradient is automatically tangent to the sphere --- which is what
#' [grad()] returns, since it projects onto the tangent space. Two conversions
#' are needed before [grad()] can be handed to [stats::optim()], both fixed by
#' comparison with central finite differences (correlation 0.9995, ratio
#' 1.000): [grad()] returns the \eqn{L^2} gradient whereas `optim()` expects
#' the Euclidean one, hence the multiplication by the trapezoidal quadrature
#' weights; and it carries a factor 2 relative to the derivative of
#' \eqn{\mathrm{loss}^2}. The barrier of [grad()] is switched off (`w = 0`):
#' \eqn{\gamma' = \psi^2 \geq 0} holds by construction, so it is not needed.
#'
#' On the vocal tract contours this is roughly an order of magnitude more
#' accurate than the steepest descent at the same cost, and its runtime is
#' almost independent of the grid size; see
#' `simulations/simulation reparametrization/`.
#'
#' @param q1j matrix with two columns: the reference SRVF.
#' @param q2j matrix of the same dimension as `q1j`: the SRVF being warped.
#' @param maxit,factr,pgtol control parameters of `optim(method = "L-BFGS-B")`.
#'
#' @return A list with elements `path` and `dpath`, exactly as [opt.gamj()].
#'
#' @seealso [opt.gamj()], [grad()], [steep_desc()].
#' @keywords internal
opt.gamj.bfgs <- function(q1j, q2j, maxit = 100, factr = 1e7, pgtol = 1e-6) {
  n   <- nrow(q1j)
  ts  <- seq(0, 1, length.out = n)
  dxh <- diff(ts) / 2
  h   <- 1 / (n - 1)
  wts <- c(h / 2, rep(h, n - 2), h / 2)   # trapezoidal quadrature weights
  spq2x <- splinefun(ts, q2j[, 1], method = "natural")
  spq2y <- splinefun(ts, q2j[, 2], method = "natural")

  prep <- function(psi) {
    nrm <- sqrt(L2prod(psi, psi))
    l <- psi / nrm
    list(l = l, nrm = nrm, gam = as.numeric(cumtrapz_vec(l^2, dxh = dxh)))
  }
  fn <- function(psi) {
    p <- prep(psi)
    loss(q1j, q2j, p$gam, p$l, spq2x = spq2x, spq2y = spq2y)^2
  }
  gr <- function(psi) {
    p <- prep(psi)
    wts * grad(p$gam, q1j, q2j, w = 0, spq2x = spq2x, spq2y = spq2y,
               l = p$l)$g / (2 * p$nrm)
  }
  o <- stats::optim(par = rep(1, n), fn = fn, gr = gr, method = "L-BFGS-B",
                    control = list(maxit = maxit, factr = factr, pgtol = pgtol))
  p <- prep(o$par)
  gam <- p$gam
  gam <- (gam - gam[1]) / (gam[n] - gam[1])
  list(path = cbind(ts, gam), dpath = p$l^2)
}


#' Parallel translation on the sphere of warpings
#'
#' Translate the tangent vector `xi` at `l` along the geodesic leaving `l`
#' with velocity `v`, i.e. the one followed by [exp_map_l()]. The component
#' of `xi` orthogonal to `v` is unchanged; the component along `v` rotates
#' with the geodesic. Parallel translation is an isometry for the
#' \eqn{L^2} metric, so it preserves the curvature pairs stored by
#' [opt.gamj.rlbfgs()].
#'
#' @param l numeric vector: a point of the unit sphere of \eqn{L^2([0,1])}.
#' @param v numeric vector: the velocity of the geodesic, tangent at `l`.
#' @param xi numeric vector: the tangent vector to translate.
#' @return The translated vector, tangent at `exp_map_l(l, v)`.
#' @seealso [exp_map_l()], [opt.gamj.rlbfgs()].
#' @keywords internal
par_transp_l <- function(l, v, xi) {
  nv <- L2norm(v)
  if (!is.finite(nv) || nv < 1e-14) return(xi)
  vh <- v / nv
  c1 <- L2prod(xi, vh)
  xi - c1 * vh + c1 * (-sin(nv) * l + cos(nv) * vh)
}

#' Riemannian quasi-Newton solver for the optimal warping
#'
#' Minimise the same alignment cost as [steep_desc()] by a Riemannian
#' limited-memory BFGS on the unit sphere of \eqn{L^2([0,1])} with the
#' trapezoidal metric --- the representation of Huang et al. (2016). The
#' exponential map [exp_map_l()] is used as retraction, the stored
#' curvature pairs are carried between tangent spaces by
#' [par_transp_l()], the update is skipped when the curvature condition
#' fails, and the step size comes from Armijo backtracking along the
#' geodesic. The Riemannian gradient is the one [grad()] returns, which is
#' already projected onto the tangent space.
#'
#' This is the manifold counterpart of [opt.gamj.bfgs()]. The two reach the
#' same solution to within a percent; the Riemannian one takes fewer
#' gradient evaluations but is slower in wall-clock time, being plain
#' \R code against the compiled `optim()`. See
#' `simulations/simulation reparametrization/`.
#'
#' @param q1j matrix with two columns: the reference SRVF.
#' @param q2j matrix of the same dimension as `q1j`: the SRVF being warped.
#' @param mem number of curvature pairs kept.
#' @param maxit maximum number of iterations.
#' @param tol stopping threshold on the norm of the Riemannian gradient.
#' @param sigma,beta_ls sufficient-decrease constant and contraction factor
#'   of the Armijo backtracking.
#'
#' @return A list with elements `path` and `dpath`, exactly as [opt.gamj()],
#'   plus the counts `nf` and `ng` of cost and gradient evaluations.
#'
#' @seealso [opt.gamj()], [opt.gamj.bfgs()], [par_transp_l()], [grad()].
#' @references
#' Huang, W., Gallivan, K. A., Srivastava, A., and Absil, P.-A. (2015).
#' Riemannian optimization for registration of curves in elastic shape analysis.
#' @keywords internal
opt.gamj.rlbfgs <- function(q1j, q2j, mem = 10, maxit = 200, tol = 1e-7,
                            sigma = 1e-4, beta_ls = 0.5) {
  n   <- nrow(q1j)
  ts  <- seq(0, 1, length.out = n)
  dxh <- diff(ts) / 2
  spq2x <- splinefun(ts, q2j[, 1], method = "natural")
  spq2y <- splinefun(ts, q2j[, 2], method = "natural")

  gam_of <- function(l) as.numeric(cumtrapz_vec(l^2, dxh = dxh))
  Cost <- function(l) loss(q1j, q2j, gam_of(l), l, spq2x = spq2x, spq2y = spq2y)^2
  Grad <- function(l) grad(gam_of(l), q1j, q2j, w = 0,
                           spq2x = spq2x, spq2y = spq2y, l = l)$g / 2

  l <- rep(1, n); l <- l / L2norm(l)
  g <- Grad(l); f <- Cost(l)
  S <- list(); Y <- list(); RHO <- numeric(0)
  nf <- 1L; ng <- 1L

  for (k in seq_len(maxit)) {
    if (!is.finite(L2norm(g)) || L2norm(g) < tol) break
    ## two-loop recursion, every inner product in the L2 metric
    qv <- g; alp <- numeric(length(S))
    if (length(S)) for (i in rev(seq_along(S))) {
      alp[i] <- RHO[i] * L2prod(S[[i]], qv); qv <- qv - alp[i] * Y[[i]] }
    j <- length(S)
    r <- (if (j) L2prod(S[[j]], Y[[j]]) / L2prod(Y[[j]], Y[[j]]) else 1) * qv
    if (length(S)) for (i in seq_along(S)) {
      bq <- RHO[i] * L2prod(Y[[i]], r); r <- r + S[[i]] * (alp[i] - bq) }
    d <- -r
    if (L2prod(d, g) >= 0) d <- -g          # safeguard: steepest descent

    a <- 1; okstep <- FALSE
    for (bt in 0:50) {
      lnew <- exp_map_l(l, a * d); fnew <- Cost(lnew); nf <- nf + 1L
      if (is.finite(fnew) && (f - fnew) >= -sigma * a * L2prod(g, d)) { okstep <- TRUE; break }
      a <- a * beta_ls }
    if (!okstep) break

    gnew <- Grad(lnew); ng <- ng + 1L; step <- a * d
    sk <- par_transp_l(l, step, step)
    yk <- gnew - par_transp_l(l, step, g)
    if (length(S)) for (i in seq_along(S)) {
      S[[i]] <- par_transp_l(l, step, S[[i]]); Y[[i]] <- par_transp_l(l, step, Y[[i]]) }
    sy <- L2prod(sk, yk)
    if (is.finite(sy) && sy > 1e-10 * L2norm(sk) * L2norm(yk)) {   # cautious update
      S[[length(S) + 1L]] <- sk; Y[[length(Y) + 1L]] <- yk; RHO <- c(RHO, 1 / sy)
      if (length(S) > mem) { S <- S[-1]; Y <- Y[-1]; RHO <- RHO[-1] } }
    l <- lnew; f <- fnew; g <- gnew
  }
  gam <- gam_of(l)
  gam <- (gam - gam[1]) / (gam[n] - gam[1])
  # `path` and `dpath` are the contract of opt.gamj(); the evaluation counts are
  # extra, and are what the comparison script reports
  list(path = cbind(ts, gam), dpath = l^2, nf = nf, ng = ng)
}


#' Optimal warping between two SRVF curves
#'
#' Minimise the elastic alignment cost over warpings, starting from the
#' identity, and return the result in the form expected by
#' [reparam.curvej()]: the graph \eqn{\{(t, \gamma(t))\}} together with the
#' derivative \eqn{\gamma'}. The warping is obtained as
#' \eqn{\gamma(t) = \int_0^t \psi^2} from the optimal \eqn{\psi}, hence is
#' automatically non-decreasing with \eqn{\gamma(0) = 0}.
#'
#' Two solvers are available. `"lbfgs"` (the default) hands the gradient to
#' L-BFGS-B via [opt.gamj.bfgs()]; `"sd"` is the original Riemannian steepest
#' descent of [steep_desc()]. On the vocal tract contours the quasi-Newton
#' solver is about an order of magnitude more accurate at the same cost, which
#' is why it is the default; `"sd"` is kept for reproducing earlier results.
#'
#' @param q1j matrix with two columns: the reference SRVF.
#' @param q2j matrix of the same dimension as `q1j`: the SRVF being
#'   warped.
#' @param eps,maxit convergence threshold and maximum number of iterations.
#'   With `optimizer = "sd"` both are passed on to [steep_desc()]; with
#'   `"lbfgs"` `eps` is ignored (L-BFGS-B has its own tolerances, see
#'   [opt.gamj.bfgs()]) and `maxit` caps its number of iterations. A
#'   small cap (1 or 2) yields a deliberately shallow alignment: the
#'   alternating elastic fits are only stable with such shallow steps
#'   (the steepest descent owes its stability to the same effect), whereas
#'   a fully converged warping makes the tangent-space alternation drift.
#'   See `align_maxit` in [sphere_regression()].
#' @param optimizer either `"lbfgs"` (default) or `"sd"`, see Details.
#'
#' @return A list with elements
#'   \describe{
#'     \item{`path`}{Matrix with two columns holding the identity grid
#'       \eqn{t} and the optimal warping \eqn{\gamma(t)}.}
#'     \item{`dpath`}{Numeric vector with the derivative
#'       \eqn{\gamma' = \psi^2}.}
#'   }
#'
#' @seealso [opt.reparam.curve()], [reparam.curvej()], [opt.gamj.bfgs()],
#'   [steep_desc()].
#' @keywords internal
opt.gamj <- function(q1j, q2j, eps = 0.05, maxit = 50,
                     optimizer = c("lbfgs", "sd", "rlbfgs")){
  optimizer <- match.arg(optimizer)
  if (optimizer == "lbfgs")  return(opt.gamj.bfgs(q1j, q2j, maxit = maxit))
  if (optimizer == "rlbfgs") return(opt.gamj.rlbfgs(q1j, q2j))

  nn = nrow(q2j)
  x  = seq(0,1,length.out=nn)
  l = steep_desc(x, q1j, q2j, eps = eps, maxit = maxit)
  gam = cumtrapz_vec(l^2, x)


  path = cbind(x, gam)
  dpath = l^2
  return(list (path=path, dpath = dpath))
}

#' Apply a warping to an SRVF curve
#'
#' Evaluate the action of a warping \eqn{\gamma} on a square-root
#' velocity function,
#' \deqn{(q \ast \gamma)(t) = q(\gamma(t))\, \sqrt{\gamma'(t)},}
#' which is the group action under which the \eqn{L^2} norm --- and hence
#' the sphere the SRVFs live on --- is preserved. The composition
#' \eqn{q \circ \gamma} is computed by evaluating a smoothing-spline
#' interpolant of each coordinate of `q` at the warped grid `path[, 2]`.
#'
#' Because the discretised action is only approximately norm-preserving,
#' the vectorised branch (`vec = TRUE`) renormalises the result with
#' respect to `S` before returning it; the matrix branch does not, and
#' callers that need an exact sphere point should normalise themselves.
#'
#' @param qj the SRVF to be warped: either a matrix with two columns, or
#'   --- if `vec = TRUE` --- a numeric vector of length \eqn{2m} in the
#'   stacked representation described in [align_vectorized_curves()].
#' @param path matrix with two columns as returned by [opt.gamj()]; only
#'   the second column, holding \eqn{\gamma(t)}, is used.
#' @param dpath numeric vector with the derivative \eqn{\gamma'}.
#' @param S a sphere geometry as returned by [initialize_sphere()], used
#'   for the renormalisation when `vec = TRUE`.
#' @param vec logical; if `TRUE`, input and output are stacked numeric
#'   vectors and the result is renormalised to unit norm.
#'
#' @return The warped SRVF, in the same representation as `qj`.
#'
#' @seealso [opt.reparam.curve()], [align_curves()].
#' @export
reparam.curvej <- function (qj, path, dpath, S = initialize_sphere(), vec = F){
  if (vec) {
    m = length(qj)/2
    qj = cbind(qj[1:m], qj[(m+1):(2*m) ] )}
  qjx = smooth.spline(seq(0,1, length.out = nrow(qj)), as.numeric(qj[,1]))
  qjy = smooth.spline(seq(0,1, length.out = nrow(qj)), as.numeric(qj[,2]))
  qj.sx = predict(qjx, path[,2])
  qj.sy = predict(qjy, path[,2])
  if (vec)  {
    curv = c(qj.sx$y*sqrt(dpath),qj.sy$y*sqrt(dpath))
    return(curv/S$norm.u(curv))
  }
  return(cbind(qj.sx$y*sqrt(dpath),qj.sy$y*sqrt(dpath)))
  #return(cbind(qj.sx$y*sqrt(path[,2]),qj.sy$y*sqrt(path[,2])))
}




#' Optimally reparametrise one SRVF curve onto another
#'
#' Convenience wrapper combining [opt.gamj()] and [reparam.curvej()]:
#' estimate the warping that best matches the parametrisation of `q2` to
#' that of `q1`, and return the warped curve alongside the warping
#' itself. This is the reparametrisation half of [align_curves()]; the
#' rotation half is [opt.O()].
#'
#' @param q1 matrix with two columns: the reference SRVF.
#' @param q2 matrix of the same dimension as `q1`: the SRVF to be
#'   reparametrised.
#' @param eps,maxit convergence threshold and maximum number of iterations
#'   of the warping optimisation; only used when `optimizer = "sd"`.
#' @param optimizer solver for the warping, passed on to [opt.gamj()]:
#'   `"lbfgs"` (the default), `"sd"` or `"rlbfgs"`.
#'
#' @return A list with elements
#'   \describe{
#'     \item{`rep_q2`}{Matrix of the same dimension as `q2`: the warped
#'       SRVF.}
#'     \item{`otp_gam`}{The warping, as returned by [opt.gamj()].}
#'   }
#'
#' @seealso [align_curves()], [opt.gamj()], [reparam.curvej()].
#' @export
opt.reparam.curve <- function (q1, q2, eps = 0.05, maxit = 50,
                               optimizer = c("lbfgs", "sd", "rlbfgs")){
  optgam = opt.gamj(q1, q2, eps = eps, maxit = maxit,
                    optimizer = match.arg(optimizer))
  rep.q2 = reparam.curvej(q2, optgam$path, optgam$dpath)
  return(list (rep_q2= rep.q2, otp_gam = optgam ))
}


#' Reconstruct and centre a shape from its stacked representation
#'
#' Turn a stacked vector back into a planar curve centred at the origin.
#' With `srvf = TRUE` (the default) the input is understood as a
#' square-root velocity function and integrated back to a curve via
#' `fdasrvf::q_to_curve()`; with `srvf = FALSE` it is understood as the
#' curve itself and merely unstacked. In both cases the column means are
#' subtracted, removing the translation that the SRVF representation does
#' not determine anyway.
#'
#' This is the standard last step before plotting or before comparing
#' fitted and observed outlines: predictions of [sphere_regression()] on
#' shape data are SRVFs, and `translate_shape()` puts them back into the
#' space where they can be drawn.
#'
#' @param q numeric vector of length \eqn{2m} in the stacked
#'   representation described in [align_vectorized_curves()].
#' @param srvf logical; if `TRUE`, `q` holds a square-root velocity
#'   function and is integrated back to a curve, otherwise it already
#'   holds curve coordinates.
#'
#' @return An `m x 2` matrix with the centred curve.
#'
#' @seealso [fdasrvf::q_to_curve()] for the package's own inverse SRVF transform,
#'   [align_vectorized_curves()].
#' @export
translate_shape = function (q, srvf = TRUE){
  m = length(q)/2
  if (srvf) {
    shape = t(q_to_curve(t(cbind(q[1:m], q[(m+1):(2*m) ] ))))} else {shape = cbind(q[1:m], q[(m+1):(2*m) ] )}

  shape - matrix(rep(colMeans(shape), each = nrow(shape)), ncol = 2)
}



#' Align a sample of curves to their Karcher mean
#'
#' Register a whole sample of SRVF curves to a common reference by
#' aligning each of them --- over rotation and reparametrisation --- to
#' the Karcher mean of the sample, as computed by [karcher_mean()]. The
#' result is the elastic analogue of centring a sample: distances to the
#' mean become elastic shape distances, and the aligned curves can be
#' mapped to a common tangent space for subsequent linear analysis (this
#' is what `cross_val_tg_space_regression()` does with
#' `elastic = TRUE`).
#'
#' Note that the mean is computed once, from the *unaligned* sample, and
#' is not iterated with the alignment.
#'
#' Only alignment to the Karcher mean is implemented; any other value of
#' `method` raises an error. In the present implementation only
#' `vectorized = TRUE` is supported.
#'
#' @param curves matrix with one stacked curve per column, of dimension
#'   \eqn{2m \times n}; see [align_vectorized_curves()] for the stacked
#'   representation.
#' @param S a sphere geometry as returned by [initialize_sphere()],
#'   referring to the stacked representation.
#' @param vectorized logical; must be `TRUE`, indicating that `curves`
#'   holds stacked vectors.
#' @param method alignment target; currently only `"karcher mean"`.
#' @param ... further arguments passed on to [align_vectorized_curves()],
#'   such as `optimizer` and `maxit` of the warping optimisation
#'   (`optimizer = "sd"` reproduces the alignment used before the
#'   L-BFGS-B solver became the default).
#'
#' @return A matrix of the same dimension as `curves`, holding the
#'   aligned curves in its columns.
#'
#' @seealso [align_vectorized_curves()], [align_curves()],
#'   [karcher_mean()].
#' @export
align_multiple_curves = function (curves, S = initialize_sphere(), vectorized = TRUE, method = c('karcher mean'), ...){
  if (method != c('karcher mean')) stop('Sorry, for now only multiple alignment to the karcher mean is implemented.')
  m = nrow(curves)/2
  curves_aligned = matrix(NA, nrow = 2*m, ncol = ncol(curves))
  km = karcher_mean(curves, S = S)
  if (!vectorized) {q1 = cbind(km[1:m], km[(m+1):(2*m)])} else {q1 = km}
  for (i in 1:ncol(curves)){
    if (!vectorized) {
      q2 = cbind(curves[1:m,i], curves[(m+1):(2*m),i])
      q2_aligned = align_curves(q1, q2, S = S, ...)$aligned_q2
      curves_aligned[,i] = c(q2_aligned[1:m], q2_aligned[(m+1):(2*m)])
    } else {
        q2 = curves[,i]
        q2_aligned =  align_vectorized_curves (q1, q2, S = S, ...)$aligned_q2
        curves_aligned[,i] = q2_aligned
        }
  }
  curves_aligned
}

#' Align every observation of a sample to its predicted counterpart
#'
#' Internal workhorse of the elastic fitting loops: align each column of
#' `Y` to the corresponding column of `y_pred` with
#' [align_vectorized_curves()], optionally distributing the (independent)
#' per-curve alignments over a cluster.
#'
#' @param y_pred matrix of stacked reference curves, one per column.
#' @param Y matrix of the same dimension holding the curves to be aligned.
#' @param S sphere geometry of the stacked representation.
#' @param eps,maxit passed on to [align_vectorized_curves()].
#' @param optimizer solver for the warping, passed on to [opt.gamj()]:
#'   `"lbfgs"` (the default), `"sd"` or `"rlbfgs"`.
#' @param cl optional cluster as returned by [make_align_cluster()]; if
#'   `NULL` (the default), the alignments run serially.
#' @return A list with the aligned sample `Y_aligned` (same dimension as
#'   `Y`) and the list `gammas` of the individual warpings.
#' @noRd
align_sample <- function(y_pred, Y, S, eps = 0.05, maxit = 50, cl = NULL,
                         optimizer = c("lbfgs", "sd", "rlbfgs")) {
  optimizer <- match.arg(optimizer)
  # everything the workers need is passed as arguments (so it is evaluated and
  # serialised as values, not as promises pointing into this frame), and the
  # function is enclosed by the package namespace so that only a reference to
  # it, not this frame, is shipped
  align_one <- function(j, y_pred, Y, S, eps, maxit, optimizer)
    align_vectorized_curves(y_pred[, j], Y[, j], S = S, eps = eps, maxit = maxit,
                            optimizer = optimizer)
  environment(align_one) <- topenv(environment())
  res <- if (is.null(cl))
    lapply(seq_len(ncol(y_pred)), align_one, y_pred = y_pred, Y = Y, S = S,
           eps = eps, maxit = maxit, optimizer = optimizer) else
      parallel::parLapply(cl, seq_len(ncol(y_pred)), align_one, y_pred = y_pred, Y = Y, S = S,
                          eps = eps, maxit = maxit, optimizer = optimizer)
  list(Y_aligned = vapply(res, function(r) as.vector(r$aligned_q2), numeric(nrow(Y))),
       gammas = lapply(res, `[[`, "opt_gam"))
}

#' Set up a cluster for parallel elastic alignment
#'
#' Create a socket cluster (the mechanism that also works on Windows,
#' where forking is unavailable) whose workers are equipped with the
#' alignment routines of this package and its dependencies, ready to be
#' passed as `cl` to [sphere_regression()],
#' [cross_val_sphere_regression()] or [tg_space_regression_model()].
#' The per-curve alignments of the elastic fitting loops are independent
#' and are then distributed over the workers.
#'
#' Each worker loads the same code as the master: the source tree when
#' the package was loaded with `devtools::load_all()`, the installed
#' package otherwise, searched in the library paths of the master. This matters because the closures shipped to the
#' workers reference the package namespace, which a worker resolves by
#' loading the package itself; exporting the functions would not help,
#' as an installed but outdated copy would still take precedence.
#'
#' Remember to release the workers with `parallel::stopCluster(cl)` when
#' done. Parallelisation pays off for the elastic fits, where each
#' alignment costs tens of milliseconds; for small samples the
#' communication overhead can outweigh the gain.
#'
#' @param ncores number of worker processes.
#' @return A cluster object as returned by [parallel::makeCluster()].
#' @export
make_align_cluster <- function(ncores) {
  ns <- asNamespace("sphereg2")
  dev <- requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("sphereg2")
  path <- getNamespaceInfo(ns, "path")   # source tree under load_all(), library path otherwise
  libs <- .libPaths()                     # so that the workers find the same installed packages
  cl <- parallel::makeCluster(ncores)
  ok <- tryCatch({
    parallel::clusterCall(cl, function(path, dev, libs) {
      .libPaths(libs)
      suppressMessages({
        if (dev) pkgload::load_all(path, quiet = TRUE) else loadNamespace("sphereg2")
      })
      NULL
    }, path, dev, libs)
    TRUE
  }, error = function(e) e)
  if (!isTRUE(ok)) {
    parallel::stopCluster(cl)
    stop("could not load the package on the workers: ", conditionMessage(ok))
  }
  cl
}
