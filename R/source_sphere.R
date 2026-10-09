
# geometric tools ---------------------------------------------------------

#' Construct a spherical Riemannian geometry
#'
#' Build a self-contained collection of geometric primitives for the unit
#' sphere \eqn{\mathbb{S} = \{y \in \mathcal{Y} : \|y\| = 1\}} of a real
#' separable Hilbert space \eqn{\mathcal{Y}}, parametrised by a user-supplied
#' inner product `prod.u`. The resulting object decouples the rest of the
#' package from any specific representation of \eqn{\mathcal{Y}}, so the same
#' code base handles, for instance, square-root density responses in
#' \eqn{L^2([0,1])} (using a trapezoidal inner product) just as readily as
#' compositional or directional data in \eqn{\mathbb{R}^{D+1}} (using the
#' Euclidean inner product).
#'
#' For points \eqn{p, q \in \mathbb{S}} and tangent vectors
#' \eqn{v \in T_p\mathbb{S} = \{v \in \mathcal{Y} : \langle v, p\rangle = 0\}},
#' the returned list contains the following functions:
#' \describe{
#'   \item{`exp.map(p, v)`}{Riemannian exponential map
#'     \eqn{\mathrm{Exp}_p(v) = \cos\|v\|\, p + \sin\|v\|\, v/\|v\|},
#'     extended by continuity at `v = 0`.}
#'   \item{`log.map(p, q)`}{Riemannian logarithmic map
#'     \eqn{\mathrm{Log}_p(q)} on \eqn{\mathbb{S}\setminus\{-p\}}; inverse of
#'     `exp.map` restricted to the open ball of radius \eqn{\pi} in
#'     \eqn{T_p\mathbb{S}}. Returns \eqn{0} when \eqn{p} and \eqn{q} are
#'     numerically equal.}
#'   \item{`proj_tg(p, v)`}{Orthogonal projection
#'     \eqn{\mathcal{Y} \to T_p\mathbb{S}}, i.e.\ removing the component
#'     along \eqn{p}.}
#'   \item{`parall.transp(p, q, v)`}{Parallel transport of
#'     \eqn{v \in T_p\mathbb{S}} along the unique minimising geodesic from
#'     \eqn{p} to \eqn{q}, yielding an element of \eqn{T_q\mathbb{S}}.}
#' }
#'
#' Numerically sensitive scalar quantities (norms, the argument of `acos`,
#' inner products of nearly-aligned vectors, ...) are routed through
#' `precise`, which by default is the identity and can be replaced by a
#' higher-precision arithmetic if needed.
#'
#' Each routine performs sanity checks (e.g.\ that the supplied point lies
#' on the sphere up to a small tolerance, that the supplied tangent vector
#' is orthogonal to the base point) and returns invisibly with an
#' informative message rather than throwing if the checks fail; callers
#' should therefore inspect return values when chaining computations.
#'
#' @param prod.u function `(q1, q2) -> numeric` computing the inner product
#'   on \eqn{\mathcal{Y}}. Defaults to the standard Euclidean inner product.
#' @param norm.u function `q -> numeric` computing the induced norm. By
#'   default \eqn{\sqrt{\langle q, q\rangle}}; supplied as a separate
#'   argument so that callers may substitute a more efficient or more
#'   accurate computation if available.
#' @param precise function applied to scalar arithmetic in numerically
#'   sensitive steps; the identity by default.
#'
#' @return An object of class `"sphere"`: a list containing `prod.u`,
#'   `norm.u`, `precise`, and the four geometric primitives `exp.map`,
#'   `log.map`, `proj_tg`, `parall.transp` described above.
#'
#' @seealso [karcher_mean()] for computing the Fréchet/Karcher mean on a
#'   sphere object; [sphere_regression()] for fitting the kernel ridge
#'   regression model.
#' @export
#' @examples
#' # Euclidean unit sphere in R^3
#' S <- initialize_sphere()
#' p <- c(1, 0, 0); q <- c(0, 1, 0)
#' v <- S$log.map(p, q)        # tangent vector at p pointing toward q
#' S$exp.map(p, v)             # recovers q (up to numerical precision)
initialize_sphere <- function(prod.u = function (q1, q2) sum(q1*q2),
                              norm.u = function(q) sqrt(prod.u(q,q)),
                              precise = function(x) x ) {
  s <- list(prod.u = prod.u, norm.u = norm.u, precise = precise)

  s$exp.map = function (p,v) {
    # spherical Riemannian exponential map
    # projects a point v in TpS to S
    if (abs(prod.u(p,v)) >= 1e-7 ) {
      print('exp.map : this vector v is not in the tangent space at p')
      return()
    }
    nv = precise(norm.u(v))
    if (nv <= 1e-10) return(p)
    exp.pv = precise(precise(cos (nv) )* p + v* (precise(sin (nv) / nv)))
    return(exp.pv)
  }

  s$log.map = function (p,q) {
    # spherical Riemannian logarithmic map
    # projects a point q in S to TpS
    np = precise(norm.u(p))

    if (np >= (1+1e-9)| np <= (1-1e-9) ) {
      print(np)
      print('log.map : p is not on the sphere')
      return()
    }
    nq = precise(norm.u(q))
    if (nq >= (1+1e-8) | nq <= (1-1e-8) ) {
      print(np)
      print('log.map : q is not on the sphere')
      return()
    }
    npq = precise(norm.u(p-q) )
    if(npq <= 1e-9) return(0 * p)
    pq = precise(prod.u(p,q))
    if (pq <= (1+1e-9) & pq >= (1-1e-9) )return(0 * p)
    theta = precise(acos (pq))
    return ( precise(theta / (sin(theta)))* precise((q- p* pq)))
  }

  s$proj_tg = function (p, v) return(v - precise(prod.u(v, p)/prod.u(p,p))*p)  # orthogonal projection: projects a point v1 in the ambient space onto TpS


  s$parall.transp = function (p, q, v){
    # parallel transport p->q
    # v is a vector in TpS that has to be transported in TqS
    np = precise( norm.u(p))
    if (np >= (1+1e-10)| np <= (1-1e-10) ) {
      print('parall.transp: p is not on the sphere')
      return()
    }
    nq = precise( norm.u(q))
    if (nq >= (1+1e-10)| nq <= (1-1e-10) ) {
      print('parall.transp: q is not on the sphere')
      return()
    }
    if(prod.u(p,v) >= 1e-10) {print('parall.transp: v not in the tg space at p')
      return()}
    return(v- 2 * precise(prod.u(v, q)/(norm.u(p+q))^2) * (p+q))
  }

  class(s) <- c("sphere", class(s))
  return(s)
}


#' Low-rank isometric coordinatisation of a Hilbert-space sample
#'
#' For a sample \eqn{y_1, \ldots, y_n} of elements of a Hilbert space
#' \eqn{\mathcal{Y}}, compute coordinates in an optimally approximating
#' finite-dimensional subspace, isometrically identified with
#' \eqn{\mathbb{R}^r}. Internally this is the Karhunen-Loève / PCA
#' decomposition of the Gram matrix \eqn{G_{ij} = \langle y_i, y_j\rangle}:
#' its eigendecomposition supplies both the coordinates of the input
#' sample (`$coos`) and an isometry between the spanned subspace and
#' \eqn{\mathbb{R}^r} (`$trafo`).
#'
#' Used by [sphere_regression()] when the user supplies `subY`, in which
#' case the response sphere is replaced by its image in
#' \eqn{\mathbb{R}^r}. This greatly reduces the size of the coefficient
#' matrix when responses are high-dimensional (e.g.\ functional data
#' evaluated on a dense grid) at minimal accuracy cost --- typically only
#' a small number of components is needed to capture the bulk of the
#' response variability. See Section 5 ("Computational aspects") of the
#' accompanying paper.
#'
#' The inner product is supplied as a separate argument rather than via a
#' `sphere` object so that the function can be reused outside the
#' spherical setting (e.g.\ for tangent-space samples, where the relevant
#' geometry is linear).
#'
#' @param Y matrix with `n` columns, each holding one element of the
#'   Hilbert space.
#' @param prod function `(y1, y2) -> numeric` computing the inner product
#'   on columns of `Y`.
#' @param lowrank specification of the approximating subspace dimension:
#'   \itemize{
#'     \item `NULL` (the default): no rank reduction beyond dropping
#'       eigenvalues that are exactly zero.
#'     \item integer \eqn{\geq 1}: directly the maximum subspace
#'       dimension \eqn{r}.
#'     \item value in \eqn{[0, 1)}: a relative tolerance in the
#'       cumulative spectrum of the Gram matrix; the smallest dimension
#'       achieving \eqn{1 - \tau \leq |\text{lowrank}|} is selected.
#'       Note that passing `0` may already trigger a non-trivial
#'       reduction due to floating-point error in the eigendecomposition.
#'   }
#'
#' @return A list with elements
#'   \describe{
#'     \item{`coos`}{Matrix of size \eqn{r \times n}: coordinates of
#'       \eqn{y_1, \ldots, y_n} in the chosen subspace, with the property
#'       that \eqn{\langle\text{coos}_{\cdot i}, \text{coos}_{\cdot j}\rangle
#'       \approx \langle y_i, y_j\rangle}.}
#'     \item{`trafo(x, back)`}{Function mapping between
#'       \eqn{\mathcal{Y}} (or its truncation) and \eqn{\mathbb{R}^r}.
#'       With `back = FALSE`, projects a Hilbert-space element to its
#'       coordinates; with `back = TRUE`, lifts coordinates back to the
#'       spanned subspace of \eqn{\mathcal{Y}}. When `back` is missing
#'       its value is guessed from `nrow(x)`; this is only unambiguous
#'       when \eqn{r \neq m}, otherwise the function errors out.}
#'   }
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(50 * 5), 50, 5)              # 5 responses in R^50
#' iso <- sphereg2:::subspace_setup(Y, function(a, b) sum(a * b))
#' dim(iso$coos)                                  # at most 5 coordinates
#' max(abs(iso$trafo(iso$coos, back = TRUE) - Y)) # lossless round trip
#' @keywords internal
subspace_setup <- function(Y, prod, lowrank = NULL) {
  # obtain base representation
  n <- ncol(Y)
  m <- nrow(Y)
  G <- matrix(nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      G[i,j] <- prod(Y[,i], Y[,j])
    }
  }
  e <- eigen(G, symmetric = TRUE)

  r <- min(sum(e$values>0), nrow(Y))
  if(!is.null(lowrank)) {
    if(lowrank<1) {
      tau <- cumsum(e$values[seq_len(r)])
      tau <- tau/tail(tau, 1)
      lowrank <- which.max(1-tau <= abs(lowrank))
    }
    r <- min(r, lowrank)
  }

  e$vectors <- e$vectors[,seq_len(r), drop = FALSE]
  e$values <- e$values[seq_len(r)]
  Z <- sweep(e$vectors, 2, sqrt(e$values), `*`)

  # obtain maps
  Zinv <- sweep(e$vectors, 2, sqrt(e$values), `/`)
  YZinv <- Y %*% Zinv

  trafo <- function(x, back) {
    x <- as.matrix(x)
    if(missing(back)) {
      if(r != m)
        back <- (nrow(x) == r) else
          stop("Since r = n, I cannot guess whether to back-transform or not.")
    }

    if(back) {
      stopifnot(nrow(x) == r)
      return(YZinv %*% x)
    }

    stopifnot(nrow(x) == m)
    n_ <- ncol(x)
    G <- matrix(nrow = n, ncol = n_)
    for(i in 1:n) {
      for(j in 1:n_) {
        G[i,j] <- prod(Y[,i], x[,j])
      }
    }
    crossprod(Zinv, G)
  }

  list(coos = t(Z), trafo = trafo)
}


#' Average tangent vector for one Karcher-mean iteration
#'
#' Given a candidate sphere point `mu` and a sample of points `Y`,
#' compute the average of the logarithmic images of the sample at `mu`:
#' \deqn{\bar v(\mu) = \frac{1}{n}\sum_{i=1}^n \mathrm{Log}_\mu(Y_i).}
#' This is the gradient direction used by the Karcher-mean fixed-point
#' iteration: \eqn{\bar v = 0} characterises the Karcher mean.
#'
#' @param Y matrix with one sphere point per column.
#' @param mu sphere point at which to evaluate the average tangent.
#' @param S a sphere geometry as returned by [initialize_sphere()].
#'
#' @return Numeric vector of length `nrow(Y)`: the average tangent vector
#'   \eqn{\bar v(\mu) \in T_\mu\mathbb{S}}.
#'
#' @seealso [karcher_mean()].
#' @keywords internal
compute_barv = function(Y, mu, S = initialize_sphere()){
  n = ncol (Y)
  vs = matrix(NA, ncol = ncol(Y), nrow = nrow(Y))
  for (i in 1:n){
    vs[,i] = S$log.map(mu, Y[,i])
  }
  rowMeans(vs)
}
#' Karcher (Fréchet) mean of a sample on the sphere
#'
#' Compute the intrinsic mean of a sample of points on the sphere,
#' \deqn{\hat\mu = \arg\min_{m \in \mathbb{S}} \sum_{i=1}^n
#'   d_{\mathbb{S}}^2(m, Y_i),}
#' where \eqn{d_{\mathbb{S}}} is the geodesic (spherical) distance. The
#' minimiser is found by a damped fixed-point iteration: starting from
#' the (renormalised) extrinsic mean, repeatedly move along
#' \eqn{\bar v(\mu) = n^{-1} \sum_i \mathrm{Log}_\mu(Y_i)} via the
#' exponential map until \eqn{\|\bar v(\mu)\|} falls below `epsilon`,
#' as described by Srivastava and Klassen (2016).
#'
#' The iterates remain on \eqn{\mathbb{S}} by construction. Convergence
#' is guaranteed for samples concentrated within a geodesic ball of
#' radius \eqn{\pi/2}; for more dispersed samples the minimiser may not
#' be unique (cf.\ Afsari, 2011, and the bounded-support assumptions in
#' the accompanying paper).
#'
#' @param Y matrix containing one sphere point per column.
#' @param S sphere geometry as returned by [initialize_sphere()].
#' @param epsilon convergence threshold: the iteration stops once
#'   `S$norm.u(barv) < epsilon`.
#' @param step_size step size for the gradient update; the default of
#'   `0.5` provides damping for stability when the sample is spread out.
#' @param maxit maximum number of iterations before giving up. The
#'   function prints a notice and returns the current iterate if this is
#'   reached.
#' @param verbose logical; if `TRUE`, print the norm of the average
#'   tangent vector at each iteration.
#'
#' @return Numeric vector of length `nrow(Y)`: an estimate of the
#'   Karcher mean.
#'
#' @references
#' Srivastava, A. and Klassen, E. P. (2016).
#' *Functional and Shape Data Analysis*. Springer Series in Statistics.
#' Springer New York. \doi{10.1007/978-1-4939-4020-2}
#'
#' @seealso [initialize_sphere()].
#' @export
karcher_mean = function (Y, S = initialize_sphere(), epsilon = 1e-3, step_size = .5, maxit = 1e2, verbose = F){
  mu_ = rowMeans(Y)
  mu = mu_ / S$norm.u(mu_)

  for (it in 1:maxit){
    barv = compute_barv (Y, mu, S = S)
    if (verbose) cat('iteration', it, 'norm bar v', S$norm.u(barv), '\n')
    if (S$norm.u(barv) < epsilon ) return (mu)
    mu = S$exp.map(mu, step_size * barv)
  }
  cat('karcher mean: out for iterations')
  return(mu)
}

#' Trapezoidal quadrature on a regular grid
#'
#' Computes the trapezoidal-rule approximation to
#' \eqn{\int f(x)\, dx} given function values on an equispaced grid.
#' Useful as the `prod.u` inner product for sphere objects representing
#' the sphere of \eqn{L^2([a, b])}: for two functions discretised on a
#' common grid `x`, \eqn{\langle f, g\rangle_{L^2} \approx
#' \mathrm{trap}(x, f \cdot g)}.
#'
#' The implementation assumes the grid is regular --- it uses only
#' `x[2] - x[1]` as the step size --- so passing an unevenly-spaced `x`
#' produces a result that is approximate at best.
#'
#' @param x numeric vector of evaluation points, assumed regularly
#'   spaced and of length `>= 2`.
#' @param fx numeric vector of function values at `x`, same length as
#'   `x`.
#' @return Numeric scalar: the trapezoidal approximation to the integral.
#' @keywords internal
trap <- function(x, fx) {
  l = length(fx)
  (x[2] - x[1]) * (fx[1]/2 + sum(fx[2:(l-1)]) + fx[l]/2)
}


# compute gradients -------------------------------------------------------

#' Build the inner fitting routine for spherical kernel ridge regression
#'
#' Internal factory that, given a sphere geometry and a choice of
#' optimiser, returns a closure `sphere.regr(Xi, P, K, Y, Z, lambda, ...)`
#' which minimises the penalised empirical risk
#' \deqn{\mathcal{R}_n(f, \lambda) = \frac{1}{n}\sum_{i=1}^n
#'   d_{\mathbb{S}}^2\bigl(\mathrm{Exp}_o(f(X_i)), Y_i\bigr)
#'   + \lambda^2 \|f\|_{\mathcal{H}}^2}
#' over the VVRKHS representation
#' \eqn{f(\cdot) = \sum_i k(X_i, \cdot)\,\xi_i}. Here `P` plays the role
#' of the intercept \eqn{o \in \mathbb{S}}, `Xi` is the coefficient
#' matrix, `K` is the (possibly low-rank-decomposed) kernel matrix, `Y`
#' is the matrix of responses, `Z` carries any identifiability or
#' low-rank constraints, and `lambda` is the regularisation parameter.
#'
#' The closure shares the following helpers via lexical scope:
#' \describe{
#'   \item{`deriv.acos(mu, i, Y)`}{Scalar coefficient \eqn{\arccos
#'     (\langle \mu_i, Y_i\rangle) / \sqrt{1 - \langle\mu_i, Y_i\rangle^2}}
#'     arising in the gradient of \eqn{d_{\mathbb{S}}^2(\cdot, Y_i)};
#'     extended by continuity at \eqn{\langle\mu_i, Y_i\rangle = 1}.}
#'   \item{`deriv.exp.map(f, i, Y, q)`}{Derivative of
#'     \eqn{v \mapsto \mathrm{Exp}_q(v)} evaluated at \eqn{f(X_i)} in the
#'     ambient space.}
#'   \item{`grad.xi(...)`}{Riemannian gradient of \eqn{\mathcal{R}_n}
#'     with respect to the coefficient matrix \eqn{\Xi}, obtained by
#'     pre-projecting the ambient gradient onto \eqn{T_q\mathbb{S}}.}
#'   \item{`grad.q(...)`}{Gradient with respect to the intercept \eqn{q},
#'     computed via Jacobi-field decomposition (tangential / orthogonal
#'     components of the parallel-transported \eqn{\mathrm{Log}}).}
#'   \item{`obj.funct(...)`}{The penalised empirical risk, the value
#'     oracle for [stats::optim()].}
#' }
#'
#' The \eqn{\Xi}-update is carried out by [stats::optim()], with the
#' method taken from `control_optimizer$method` (default `"BFGS"`). When
#' `fix_p = FALSE` the intercept is updated by a separate Riemannian
#' gradient step before each `optim()` call, and the coefficients are
#' parallel-transported accordingly.
#'
#' # Elastic fits
#'
#' With `elastic = TRUE` the responses are treated as square-root
#' velocity functions of planar curves (see [fdasrvf::curve_to_q()]) and the fit
#' alternates between two steps until the fitted values stop moving:
#' \enumerate{
#'   \item *alignment*: every response is re-aligned to the current
#'     prediction \eqn{\mathrm{Exp}_p(f(X_i))} over rotations and
#'     reparametrisations with [align_vectorized_curves()], so that the
#'     residual measured afterwards is an elastic shape distance rather
#'     than an \eqn{L^2} distance between fixed parametrisations;
#'   \item *fitting*: the coefficients are updated against the re-aligned
#'     responses exactly as in the non-elastic case.
#' }
#' Because alignment is only meaningful in the curves' native
#' representation, `coordinates = TRUE` triggers a round trip through
#' `iso$trafo()` in each iteration: responses and intercept are lifted
#' back to the full representation, aligned, and projected down again.
#' This makes the elastic fit substantially more expensive than the
#' plain one --- hence the reduced inner `maxit` passed to
#' [stats::optim()] when `elastic = TRUE`. In the non-elastic case with
#' `fix_p = TRUE` there is nothing to alternate and the outer loop is
#' short-circuited to a single pass.
#'
#' @param S_intern sphere geometry the optimiser actually works in. When
#'   [sphere_regression()] fits in coordinates this is the Euclidean
#'   geometry of \eqn{\mathbb{R}^r}, otherwise it coincides with `S`.
#'   Its components are made available to the returned closure by
#'   lexical scope.
#' @param S sphere geometry of the responses in their native
#'   representation, as returned by [initialize_sphere()]. Only used for
#'   the elastic alignment step, which has to be carried out in that
#'   representation.
#' @param fix_p logical; if `TRUE` (default) the intercept `P` is held
#'   fixed throughout the fit (e.g.\ pinned to the Karcher mean of the
#'   responses, the usual choice in the paper).
#'   The option `fix_p=FALSE` is experimental and for future development.
#' @param control_optimizer list of control parameters passed verbatim to
#'   [stats::optim()] (so `method`, `control`, etc.).
#' @param maxit maximum number of outer iterations, i.e. of
#'   alignment/fitting alternations when `elastic = TRUE` and of
#'   intercept updates when `fix_p = FALSE`. `maxit = 0` returns the
#'   inputs unchanged.
#' @param eps_p convergence threshold on \eqn{\|\nabla_q \mathcal{R}\|},
#'   the gradient norm of the intercept update.
#' @param eps_diff convergence threshold on the largest change of any
#'   fitted value \eqn{f(X_i)} between two successive outer iterations.
#' @param eps_xi convergence threshold on the largest column norm of the
#'   gradient with respect to \eqn{\Xi}.
#' @param stp_p step size for the `P`-update; ignored when
#'   `fix_p = TRUE`.
#' @param elastic logical; if `TRUE`, alternate the fit with an elastic
#'   re-alignment of the responses, as described above.
#' @param coordinates logical; `TRUE` if the optimiser works in the
#'   low-dimensional coordinates supplied by `iso` rather than in the
#'   responses' native representation.
#' @param iso the output of [subspace_setup()] when `coordinates = TRUE`,
#'   supplying the isometry used to move between the two representations
#'   during elastic alignment.
#' @param align_eps,align_maxit convergence threshold and maximum number
#'   of iterations of the per-curve warping optimisation in the elastic
#'   alignment step, passed on to [align_vectorized_curves()] /
#'   [steep_desc()]. The defaults reproduce the previous fixed settings.
#' @param align_optimizer solver of the per-curve warping optimisation,
#'   `"lbfgs"` (default) or `"sd"`, see [opt.gamj()].
#' @param cl optional cluster as returned by [make_align_cluster()]; if
#'   supplied, the per-curve alignments of the elastic step are
#'   distributed over its workers.
#'
#' @return A function `sphere.regr(Xi, P, K, Y, Z, lambda, verbose,
#'   ...)` returning a list with components `Xi` (fitted coefficients)
#'   and `P` (fitted intercept, possibly equal to the input when
#'   `fix_p = TRUE`). Depending on the optimiser it additionally carries
#'   convergence diagnostics (`iterations`, `norm_pgradient`,
#'   `max_gradnorm_diffsuccit`) and, for elastic fits, the re-aligned
#'   responses `Y_elastic_align` together with the warpings
#'   `gamma_alignments` applied in each iteration.
#'
#' @keywords internal
initialize_sphere_fit = function(S_intern = initialize_sphere(), S = S_intern,
                                 fix_p = TRUE,
                                 control_optimizer = list(method = "BFGS"),
                                 maxit = 100, eps_p = 1e-6, eps_diff = 1e-6, eps_xi = 1e-6,
                                 stp_p = .1, elastic = FALSE, coordinates = FALSE, iso = c(),
                                 align_eps = 0.05, align_maxit = 50, cl = NULL,
                                 align_optimizer = c("lbfgs", "sd")) {
  align_optimizer <- match.arg(align_optimizer)

  prod.u <- S_intern$prod.u; norm.u <- S_intern$norm.u; precise <- S_intern$precise
  exp.map <- S_intern$exp.map; log.map <- S_intern$log.map
  proj_tg <- S_intern$proj_tg; parall.transp <- S_intern$parall.transp

      # get sphere geometry tools
       deriv.acos = function (mu, i, Y){ # derivaive of the acos^2 part of the loss function
         muy = prod.u(mu[,i], Y[,i])
         if (muy <= (1+1e-13) & muy >= (1-1e-13)) return (1) # lim_{x-> 1} (acos(x)/sqrt(1-x)) = 1
         return(precise(acos(muy)/sqrt(1-muy^2)))
       }

       deriv.exp.map = function (f, i, Y, q){ # derivative of the exp map
         nf = precise(norm.u(f[,i]))
         fy = precise(prod.u(f[,i], Y[,i]))
         return ( precise(sin (nf)/nf) * (-prod.u(q, Y[,i]) * f[,i] + Y[,i] -precise(fy/nf^2) * f[,i] )
                  + precise(cos(nf)/(nf^2)  )* fy * f[,i]) # we have to project because we are computing the gradient in an embedded manifold
       }

       grad.xi= function (Xi, K, Y, q, Z=NULL, lambda = 1e-2){ # compute the gradient of the loss function wrt Xi_l
         n = nrow(K)
         m = nrow(Xi)
         f = Xi %*% t(K) # list that contains f(X_i)
         mu = matrix(NA, nrow = m, ncol = n)
         der.acos = array(NA, n)
         der.exp.map = matrix(NA, nrow = m, ncol = n)
         for (i in 1:n){
           mu[,i] = exp.map(q, f[,i])
           der.acos [i] = deriv.acos(mu, i, Y)
           der.exp.map [,i] = proj_tg(q, as.numeric(der.acos [i]) * deriv.exp.map (f, i, Y, q))
         }
         # return gradient of empirical risk (i.e. dividing by n)
         if(is.null(Z)) return (( -  der.exp.map %*% K + lambda * f )/n)
         return (( -  der.exp.map %*% K + lambda * f %*% Z)/n) # regularization
       }

       grad_xi_srvf_fix_p = function (Xi, K, Y, q,  lambda = 1e-2){ # compute the gradient of the loss function wrt Xi_l
         grad.xi(Xi, K, Y, q, Z=NULL, lambda = lambda)
       }

       grad.q = function (Xi, K, Y, q){ # compute the gradient of the loss function wrt q using Jacobi fields
         n = nrow(K)
         f =  Xi  %*% t(K) # contains f(X_i)
         J = sapply(1:n, function(i) {
           g.dist2 = parall.transp(exp.map(q, f[,i]), q, log.map(exp.map(q, f[,i]), Y[, i]))
           if (prod.u(f[,i],f[,i]) == 0) g.tang.fi = f[,i] * 0
           else g.tang.fi = f[,i] * prod.u(f[,i], g.dist2)/prod.u(f[,i],f[,i])
           g.perp.fi = g.dist2 - g.tang.fi
           g.tang.fi + g.perp.fi * (  precise(cos (norm.u(f[,i]))) )})
         # return gradient of empirical risk (i.e. dividing by n)
         return (rowSums( - J )/n) # - sum_{i=1}^n cos(norm.u(f_i)) PT( perp log_{exp_qf_i}(y_i) ) +  PT( tang log_{exp_qf_i}(y_i) )
       }

       obj.funct = function (Y, K, Xi, q, Z, lambda = 0){
         # computes the objective function that we want to minimize
         n = nrow(K)
         m = ncol(Xi)

         fx = tcrossprod(Xi, K)
         nf2 = if(is.null(Z))
           sum ( sapply(1:n, function(i) crossprod(fx[,i], Xi[,i])) ) else {
            fxZ = fx %*% Z
             sum ( sapply(1:m, function(i) crossprod(fxZ[,i], Xi[,i])) )
           }
         return( sum ( sapply(1:n, function(i) norm.u(log.map(exp.map(q, fx[,i]), Y[,i]))^2 )) + lambda * nf2 )
       }


          sphere.regr = function (Xi, P, K, Y, Z = NULL, lambda,
                                  verbose = TRUE, eps_diff = 1e-8, eps_xi = 1e-8, ...){

            if(maxit==0) # immediately done if  maxit = 0
              return( list(Xi = Xi, P = P) )

            if(fix_p & !elastic) maxit <- 1


            for ( it in 1:maxit){
              pri = (it %% 3 == 0)


              gamma_alignments = list()
              Y_elastic_align =  NA
              # optimization step over p -------------------------------

              if(!fix_p) {
                # gradient step for intercept P
                gp = grad.q(Xi, K, Y, P)
                np = norm.u(gp)
                Pgood = (np < eps_p)
                p = exp.map(P, - stp_p * gp)
                # transport the Xi
                Xi = sapply(1:ncol(Xi), function(i) parall.transp(P, p, Xi[,i]))
              } else {
                Pgood <- TRUE
                np = 0
                p <- P
                }

              # elastic aligment if needed -------------------------------



              if (elastic) {
                # transform back to original dimension
                if (coordinates) {
                  Y <- iso$trafo(Y, back = TRUE)
                  Y <- sweep(Y, 2, apply(Y,2,S$norm.u), `/`)
                  p <- iso$trafo(p, back = TRUE)
                  }

                # align srvf
                fx_pred <- Xi %*% t(K)
                fx = fx_pred
                if (coordinates) fx_pred <- iso$trafo(fx_pred, back = TRUE)
                y_pred <- sapply(1:ncol(fx_pred), function(b) S$exp.map(p, fx_pred[, b]))

                aligned = align_sample(y_pred, Y, S = S, eps = align_eps, maxit = align_maxit, cl = cl, optimizer = align_optimizer)
                Y_aligned = aligned$Y_aligned
                gamma_alignments[[it]] = aligned$gammas
                Y_elastic_align = Y_aligned

                # back to low dimension
                if (coordinates) {
                  Y_aligned <- iso$trafo(Y_aligned, back = FALSE)
                  Y <- sweep(Y_aligned, 2, apply(Y_aligned,2,S_intern$norm.u), `/`)
                  p <- iso$trafo(p, back = FALSE)}else{
                  Y<- Y_aligned
                }
              }





              # optimization over Xi  ----------------------------------

              fn <- function(xi) { # xi <- c(Xi)
                Xi <- array(xi, dim = dim(Xi))
                Xi = sapply(1:ncol(Xi), function(i) proj_tg(p, Xi[,i]))

                obj.funct(Y = Y, K = K, Xi = Xi, q = p, Z = Z, lambda = lambda)
              }
              gr <- function(xi) {
                Xi <- array(xi, dim = dim(Xi))
                Xi = sapply(1:ncol(Xi), function(i) proj_tg(p, Xi[,i]))
                g <- grad.xi(array(Xi, dim = dim(Xi)), Y = Y, q = p, K = K, Z = Z, lambda = lambda)
                g <- array(g, dim = dim(Xi))
                c(sapply(1:ncol(Xi), function(i) proj_tg(p, g[,i])))
              }

              control_optimizer[c("par", "fn", "gr")] <- list(
                par = c(Xi),
                fn  = fn,
                gr  = gr
              )

              control_optimizer$method  <- "L-BFGS-B"
              control_optimizer$control <- list(
                maxit = if (elastic) 10 else 100,
                trace = if (pri) 1 else 0,
                factr = 1e7,
                pgtol = 1e-6
              )


              o <- do.call(optim, control_optimizer)

              Xi <- array(o$par, dim = dim(Xi))
              Xi = sapply(1:ncol(Xi), function(i) proj_tg(p, Xi[,i]))

              gxi = grad.xi(Xi = Xi, Y = Y, q = p, K = K, Z = Z, lambda = lambda)
              maxnxi = max(apply(gxi,2,norm.u))
              Xigood = maxnxi  < eps_xi

              # return result -----------------------------------------
              if (elastic) {
                fx0 <- Xi %*% t(K)
                max_norm_diff = max (sapply(1:ncol(fx0), function (i) norm.u( fx0[,i]- fx[,i] )))
                if (pri & verbose)  print(paste0('iteration ', it, ' max norm successive iterations f(x): ' , max_norm_diff, ' max grad xi norm ', maxnxi))# ,  '; loss function: ', obj.funct(Y, K, Xi, p)))

                if(Pgood & ((max_norm_diff< eps_diff) | Xigood)) return( list(Xi = Xi, P = p, norm_pgradient = np, iterations = it, Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments))

              } else {
                if(Pgood) return( list(Xi = Xi, P = p, norm_pgradient = np, iterations = it, Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments))
              }

              P <- p

            }
            if(verbose) cat('Out for iterations ', ' \n Iteration: ', it, '\n Max norm successive iterations f(x): ' ,  max_norm_diff, '\n Max grad xi norm: ', maxnxi, '\n')
            return( list(Xi = Xi, P = p, norm_pgradient = np, iterations = it,Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments))
          }

  return(sphere.regr)
} # function initialize_sphere_fit



#' Spherical kernel ridge regression
#'
#' Fit a kernel ridge regression model for responses on a sphere of
#' arbitrary (possibly infinite) dimension, given Polish-space
#' covariates. The conditional Fréchet mean of \eqn{Y \in \mathbb{S}}
#' given \eqn{X = x} is modelled intrinsically as
#' \deqn{\hat\mu(x) = \mathrm{Exp}_o\!\bigl(\hat f(x)\bigr),
#'   \qquad \hat f(\cdot) = \sum_{i=1}^n k(X_i, \cdot)\,\hat\xi_i \in
#'   \mathcal{H},}
#' where \eqn{\mathcal{H}} is a vector-valued RKHS of functions
#' \eqn{\mathcal{X} \to T_o\mathbb{S}} with separable kernel
#' \eqn{K(x, x') = k(x, x') \cdot \mathrm{id}}, and the coefficients
#' \eqn{\hat\xi_i} minimise the penalised empirical risk
#' \deqn{\mathcal{R}_n(f, \lambda) = \frac{1}{n}\sum_{i=1}^n
#'   d_{\mathbb{S}}^2\bigl(\mathrm{Exp}_o(f(X_i)), Y_i\bigr)
#'   + \lambda^2 \|f\|_{\mathcal{H}}^2.}
#' By default the intercept \eqn{o} is set to the Karcher mean of the
#' responses (see `start.p` to override).
#'
#' Computationally, the fit may be performed in two ways:
#' \itemize{
#'   \item In the response's native representation (`subY = NULL`).
#'   \item In low-dimensional coordinates obtained from the
#'     Karhunen-Loève decomposition of the responses' Gram matrix
#'     (`subY` non-`NULL`; the default whenever \eqn{m > n}, i.e.\
#'     the response dimension exceeds the sample size). This is
#'     usually a near-lossless approximation and avoids carrying the
#'     full response dimension into the optimiser, which is the main
#'     practical speed-up in high-dimensional or functional settings.
#' }
#' A complementary low-rank approximation of the kernel matrix `K` is
#' available via `subK`. Both reductions follow Section 5
#' ("Computational aspects") of the accompanying paper.
#'
#' The optimiser is selected via `optimizer`; the default `"optim"`
#' wraps base R's BFGS, which we found to be the most reliable choice.
#'
#' [cross_val_sphere_regression()] performs k-fold cross-validation on
#' a grid of kernel/penalty pairs to select tuning parameters.
#'
#' # Elastic shape responses
#'
#' Setting `elastic = TRUE` switches from an \eqn{L^2} fit of fixed
#' parametrisations to an elastic fit of shapes: responses are taken to
#' be square-root velocity functions of planar curves (see [fdasrvf::curve_to_q()])
#' and, in every iteration, re-aligned to the current prediction over
#' rotations and reparametrisations before the coefficients are updated.
#' The fitted model therefore approximates the conditional Fréchet mean
#' with respect to the elastic shape distance rather than the ambient
#' spherical one. See [initialize_sphere_fit()] for the details of the
#' alternation and [align_vectorized_curves()] for the alignment step;
#' [tg_space_regression_model()] provides the tangent-space counterpart
#' used for comparison in the simulations.
#'
#' @param X covariate matrix of size \eqn{n \times d} (or a vector,
#'   coerced to a one-column matrix).
#' @param Y response matrix of size \eqn{m \times n} with one
#'   observation on the sphere per column.
#' @param k a kernel function with signature `(x1, x2)`; the default is
#'   the Gaussian radial basis function, whose bandwidth argument
#'   `sigma2` is only ever taken at its default value, so a non-default
#'   bandwidth has to be baked into the closure (see
#'   [create_kernel_list()]). Alternatively, an already-evaluated
#'   \eqn{n \times n} kernel matrix; in that case the fit uses it
#'   directly and `model$k` is set to `NA`.
#' @param prod.u inner product defining the response geometry, passed
#'   to [initialize_sphere()] when `S` is `NULL`.
#' @param lambda penalty parameter \eqn{\lambda} in
#'   \eqn{\mathcal{R}_n(f, \lambda)}.
#' @param fix_p logical; if `TRUE` (default) the intercept is held
#'   fixed at `start.p`. Setting `FALSE` enables a joint optimisation
#'   over the intercept and the coefficients with the corresponding
#'   parallel transport (see [initialize_sphere_fit()]).
#' @param subK,subY low-rank approximation thresholds for the kernel
#'   matrix and for the response Gram matrix, respectively; see
#'   [subspace_setup()] for the encoding (`NULL` = no reduction,
#'   integer = subspace dimension, value in \eqn{[0, 1)} = relative
#'   spectral tolerance). The default `subY = if(nrow(Y) > ncol(Y)) 0`
#'   enables an automatic, near-lossless reduction in the
#'   high-dimensional regime.
#' @param S optional sphere geometry as returned by
#'   [initialize_sphere()]; if `NULL`, one is constructed from
#'   `prod.u`.
#' @param verbose logical; if `TRUE`, progress information is printed
#'   during fitting.
#' @param start.xi optional starting value for the coefficient matrix
#'   \eqn{\Xi}. If `NULL`, the closed-form tangent-space ridge
#'   solution at `start.p` is used.
#' @param start.p optional starting value for the intercept; if
#'   `NULL`, the Karcher mean of `Y` is used.
#' @param control_optimizer list of control parameters passed to
#'   [stats::optim()]; the default `list(method = "BFGS")` selects BFGS.
#' @param maxit maximum number of outer iterations, i.e. of
#'   alignment/fitting alternations when `elastic = TRUE` and of
#'   intercept updates when `fix_p = FALSE`.
#' @param eps_p,stp_p convergence threshold and step size for the
#'   intercept update; ignored when `fix_p = TRUE`.
#' @param eps_diff convergence threshold on the largest change of any
#'   fitted value \eqn{f(X_i)} between two successive outer iterations.
#' @param eps_xi convergence threshold on the largest column norm of the
#'   gradient with respect to \eqn{\Xi}.
#' @param elastic logical; if `TRUE`, fit shapes elastically by
#'   re-aligning the responses to the current prediction in every
#'   iteration, as described above.
#' @param align_eps,align_maxit convergence threshold and maximum number
#'   of iterations of the per-curve warping optimisation in the elastic
#'   alignment step (see [opt.gamj()]). `align_eps` only concerns the
#'   `"sd"` solver; `align_maxit` caps either solver. The defaults
#'   reproduce the previous fixed settings. With `"lbfgs"`, a small cap
#'   such as `align_maxit = 2` gives a shallow alignment comparable to
#'   the steepest descent at a fraction of its cost; a fully converged
#'   alignment is known to destabilise the alternating tangent-space fit
#'   of [tg_space_regression_model()].
#' @param align_optimizer solver of the per-curve warping optimisation:
#'   `"lbfgs"` (default, L-BFGS-B) or `"sd"` (the Riemannian steepest
#'   descent used before), see [opt.gamj()]. The two can settle in
#'   different local optima, so alternating fits may behave differently
#'   under the two solvers.
#' @param cl optional cluster as returned by [make_align_cluster()]; if
#'   supplied, the per-curve alignments of the elastic step are
#'   distributed over its workers. Remember to call
#'   `parallel::stopCluster()` when done.
#' @references
#' Matteo, B., Stöcker, A., and Tavakoli, S. (2026).
#' *Infinite-Dimensional Spherical Kernel Ridge Regression.*
#' arXiv preprint [arXiv:2606.00181](https://arxiv.org/abs/2606.00181).
#'
#' @return An object of class `"sphere_model"`: a list with components
#'   \describe{
#'     \item{`Xi`}{Coefficient matrix in the original response
#'       representation.}
#'     \item{`P`}{Fitted intercept on the sphere.}
#'     \item{`Xi_intern`, `P_intern`}{Coefficient matrix and intercept
#'       in the internal (possibly truncated) representation actually
#'       used by the optimiser.}
#'     \item{`S`, `S_intern`}{Original and internal sphere geometries.}
#'     \item{`K`, `Z_constrain`}{Kernel matrix and constraint matrix in
#'       the internal representation.}
#'     \item{`iso`}{If a coordinate reduction was applied, the output
#'       of [subspace_setup()] for use by [predict.sphere_model()].}
#'     \item{`k`, `lambda`, `X`}{The kernel (or `NA` if `k` was supplied
#'       as a matrix), penalty, and design matrix used in the fit.}
#'     \item{`Y_original`}{The responses as passed in, before any
#'       coordinate reduction or elastic re-alignment; used by
#'       [predict.sphere_model()] with `elastic = TRUE`.}
#'     \item{`subKdim`, `subYdim`}{Effective subspace dimensions of
#'       the kernel and response low-rank approximations. `subKdim`
#'       additionally carries the retained spectral share of `K` in its
#'       `"share"` attribute.}
#'     \item{`Y_elastic_align`, `gamma_alignments`}{For elastic fits,
#'       the re-aligned responses of the last iteration and the warpings
#'       applied in each iteration.}
#'   }
#' @export
#'
#' @examples
#' # square-root densities as responses: the rows of a smoothed surface
#' S <- initialize_sphere()
#' X <- seq(1, ncol(volcano), by = 6)
#' Y <- sqrt(volcano[, X]); Y <- sweep(Y, 2, apply(Y, 2, S$norm.u), "/")
#' fit <- sphere_regression(X, Y, k = function(a, b) exp(-sum((a - b)^2) / 50),
#'                          lambda = 0.1, verbose = FALSE)
#' pred <- predict(fit, newdata = matrix(X[-1] - 3))
#' apply(pred, 2, S$norm.u)                       # predictions lie on the sphere
sphere_regression = function (X, Y,
                              k = function(x1, x2, sigma2 = 1) exp(-sum((x1 - x2)^2) / (2 * sigma2)),
                              prod.u =  function (q1, q2) sum(q1*q2),
                              lambda = .1, fix_p= TRUE,
                              subK = 1e-6, subY = if(nrow(Y)>ncol(Y)) 0,
                              S = NULL,
                              verbose = TRUE, start.xi = NULL, start.p = NULL,
                              control_optimizer = list(method = "BFGS"),
                              maxit = 100, eps_p = 1e-6, eps_diff = 1e-8, eps_xi = 1e-8, stp_p = .1,
                              elastic = FALSE, align_eps = 0.05, align_maxit = 50, cl = NULL,
                              align_optimizer = c("lbfgs", "sd")){
  Y_original = Y
  align_optimizer <- match.arg(align_optimizer)

  coordinates <- !is.null(subY)

  if(!is.matrix(X)) X <- as.matrix(X)
  n <- ncol(Y)
  if (is.matrix(k)) K_ = k else K_ = sapply(1:n, function(i) sapply(1:n, function(j) k(X[i,], X[j,])))

  # initialize sphere geometry
  if(is.null(S))  S = initialize_sphere(prod.u = prod.u)

  iso = c()
  # if working in coordinates, transform the input
  if(coordinates) {
    iso <- subspace_setup(Y, S$prod.u, subY)
    Y <- iso$coos
    if(!is.null(start.xi))
      start.xi <- iso$trafo(start.xi, back = FALSE)
    S_intern <- initialize_sphere()
    if(!is.null(start.p)) {
      start.p <- iso$trafo(start.p, back = FALSE)
      start.p <- start.p/S_intern$norm.u(start.p)
    }
    # renormalize (for potential low-rank approximation)
    Y <- sweep(Y, 2, apply(Y,2,S_intern$norm.u), `/`)
  } else {
    S_intern <- S
  }

  # initialize fitting function
  sphere.regr = initialize_sphere_fit(S_intern, S, fix_p = fix_p,
                                      control_optimizer = control_optimizer,
                                      maxit = maxit,
                                      eps_p = eps_p, eps_diff = eps_diff, eps_xi = eps_xi,
                                      stp_p = stp_p, elastic = elastic, coordinates = coordinates, iso = iso,
                                      align_eps = align_eps, align_maxit = align_maxit, cl = cl,
                                      align_optimizer = align_optimizer)

    # get start.p and start.xi
  if ( is.null(start.p)) start.p <- karcher_mean(Y, S = S_intern)
  logY = sapply(1:n, function(i) S_intern$log.map(start.p, Y[, i]))



  # low rank approximation of K
  Z <- K <- NULL
  lowrank <- subK
  if(!is.null(lowrank)) {
    e <- eigen(K_, symmetric = TRUE)
    total <- sum(abs(e$values))
    if(lowrank<1) {
      r <- sum(e$values>0)
      e$values <- e$values[seq_len(r)]
      tau <- cumsum(e$values/sum(e$values))
      tau[length(tau)] <- 1 # avoid numerical issues
      lowrank <- which.max(1-tau <= abs(lowrank))
      if (lowrank == 1)
        warning('Rank of K is reduced to one dimension only! Please double-check!')
    }
    e$vectors <- e$vectors[, seq_len(lowrank), drop = FALSE]
      e$values <- e$values[seq_len(lowrank), drop = FALSE]
      attr(lowrank, "share") <- sum(abs(e$values)) / total
      Z <- sweep(e$vectors, 2, sqrt(e$values), `/`)
      K <- sweep(e$vectors, 2, sqrt(e$values), `*`)
  }


  # implement sum-to-zero identifiability constraint when p is not fixed
  if (fix_p) {
    if(is.null(K)) K <- K_
  } else {
    if(is.null(Z)) { # we want to project the columns of K_ into the orthogonal of the constant 1
      C = colSums(K_) # corresponds to the inner product with the constant: C = 1^T %*% K in R^(1xp)
      if (all(abs(C) < 1e-15)) {
        Z <- NULL
        K <- K_
      } else {
        Z = MASS::Null(C)  # make sure it's orthogonal to the constant: Z contains a basis of the null space of C: CZ =0, Z in R^(px(p-1))
        lf <- lsfit(K_, rep(1, nrow(K_)), intercept = FALSE)
        if (all(abs(lf$residuals) > 1e-15)) {
          Z <- cbind(Z, t(C))
        }
        K = K_ %*% Z # constrain to sum to 0: 1^TK = 0
      }

    } else {
      C = colSums(K)
      if(all(abs(C) > 1e-15)) {
        Z_ <- MASS::Null(C)
        lf <- lsfit(K, rep(1, nrow(K)), intercept = FALSE)
        if(all(abs(lf$residuals) > 1e-15)) {
          Z_ <- cbind(Z_, t(C))
        }
        Z <- Z %*% Z_
        K <- K %*% Z_
      }
    }
  }

  if(is.null(start.xi)) {
    if(is.null(Z)) {
      K_ <- chol(K + lambda * diag(n))
      start.xi = t(backsolve(K_, forwardsolve(t(K_), t(logY))))
      # faster version of: start.xi = logY %*% solve(K + lambda * diag(n))
    }  else {
        if(fix_p) { # i.e., if the only constraint matrix is due to low-rank approximation
          start.xi = logY %*% scale(K, scale = e$values[seq_len(lowrank)] + lambda, center = FALSE)
        } else {
          start.xi = logY %*% K %*% solve(crossprod(K) + lambda * crossprod(Z, K_) %*% Z)
        }
      }
    # make sure it is projected into the tangent space
    # start.xi <- apply(start.xi, 2, S$proj_tg, p = start.p)
  }


  model = sphere.regr(Xi = start.xi, P = start.p, K = K, Y = Y, Z = Z,
                      lambda = lambda, verbose = verbose, maxit = maxit, elastic = elastic)

  model$Z_constrain = Z
  model$K = K
  model$S = S
  model$S_intern <- S_intern
  model$P_intern <- model$P
  model$Xi_intern <- model$Xi
  model$subYdim <- nrow(model$Xi)
  if(coordinates) {
    model$iso <- iso
    model$P <- iso$trafo(model$P, back = TRUE)
    model$Xi <- iso$trafo(model$Xi, back = TRUE)
  }
  if (is.matrix(k)) model$k = NA else model$k = k
  model$lambda = lambda
  model$X = X
  model$Y_original = Y_original
  model$subKdim <- lowrank
  class(model) = c('sphere_model', class(model))
  return(model)
}



#' Predict from a spherical kernel ridge regression fit
#'
#' Evaluate a fitted [sphere_regression()] model, either at the training
#' covariates or at new ones. Two scales are available:
#' \describe{
#'   \item{`type = "link"`}{the tangent-space prediction
#'     \eqn{\hat f(x) = \sum_i k(X_i, x)\, \hat\xi_i \in T_o\mathbb{S}},
#'     i.e. the linear predictor before mapping to the sphere;}
#'   \item{`type = "response"`}{the prediction on the sphere itself,
#'     \eqn{\hat\mu(x) = \mathrm{Exp}_o(\hat f(x))}, which is the
#'     estimated conditional Fréchet mean.}
#' }
#'
#' Predictions are computed in the model's internal representation and,
#' unless `coordinates = TRUE`, lifted back to the responses' native
#' representation via the isometry stored in `object$iso`. Requesting
#' `coordinates = TRUE` is mainly useful inside cross-validation, where
#' the error is evaluated in coordinates anyway and the round trip would
#' be wasted work.
#'
#' For shape responses, `elastic = TRUE` additionally aligns each
#' prediction to the corresponding observed curve over rotations and
#' reparametrisations, so that the resulting residual is an elastic shape
#' distance. Since this requires the observed counterpart of every
#' prediction, it only applies to in-sample predictions and is silently
#' ignored when `newdata` is given.
#'
#' @param object a fitted model of class `sphere_model`, as returned by
#'   [sphere_regression()].
#' @param newdata optional matrix of new covariates with the same number
#'   of columns as the design matrix used in the fit. If `NULL` (the
#'   default), predictions are returned at the training covariates,
#'   reusing the stored kernel matrix.
#' @param type one of `"response"` (default) for predictions on the
#'   sphere, or `"link"` for the tangent-space linear predictor.
#' @param coordinates logical; if `TRUE`, return the prediction in the
#'   internal low-dimensional coordinates instead of the responses'
#'   native representation. Has no effect for models fitted without a
#'   coordinate reduction.
#' @param elastic logical; if `TRUE` and `newdata` is `NULL`, align each
#'   prediction to the corresponding observed curve with
#'   [align_vectorized_curves()].
#' @param ... further arguments, currently ignored.
#'
#' @return A matrix with one prediction per column: sphere points for
#'   `type = "response"`, tangent vectors at `object$P` for
#'   `type = "link"`.
#'
#' @seealso [sphere_regression()], [cross_val_sphere_regression()],
#'   [predict.tg_model()] for the tangent-space comparison model.
#' @export
predict.sphere_model <- function (object, newdata = NULL, type = c('response', 'link'), coordinates = FALSE, elastic = FALSE, ...){

  type = match.arg(type)

  if (is.null(newdata)){
    fx = object$Xi_intern %*% t(object$K)
  }else{
    if(is.null(ncol(newdata)))
      newdata <- as.matrix(newdata)
      stopifnot(ncol(newdata) == ncol(object$X))

    k_test = sapply(1:nrow(object$X), function(i) sapply(1:nrow(newdata), function(j) object$k(object$X[i,], newdata[j,])))

    if (!is.null(object$Z)) k_test = k_test %*% object$Z

    fx = object$Xi_intern %*% t(k_test)
  }

  if (type == 'link') {
    if(!coordinates & !is.null(object$iso))
      fx <- object$iso$trafo(fx, back = TRUE)
    return(fx)
  }

  y <- sapply(1:ncol(fx), function (j) object$S_intern$exp.map( object$P_intern,  fx [,j]))
  if(!is.matrix(y)) y <- t(y)
  if(!coordinates & !is.null(object$iso))
    y <- object$iso$trafo(y, back = TRUE)

  if (elastic) {
    if (is.null(newdata)){
      y = sapply ( 1: ncol(y), function (j) align_vectorized_curves (object$Y_original[,j] ,y[,j] , object$S)$aligned_q2)
    }
  }
  y
}



# Cross Validation --------------------------------------------------------

#' Evaluate an expression under a temporary random seed
#'
#' Sets the seed, evaluates `expr` and then restores the RNG state that
#' was in place before the call, so that a reproducible internal draw
#' (such as the choice of the folds to fit in the cross-validation
#' functions) does not restart the random stream of the caller. Without
#' this, simulation scripts that call the cross-validation repeatedly
#' would generate identical data in every repetition after the first.
#'
#' @param seed integer passed to [set.seed()].
#' @param expr expression to evaluate.
#' @return The value of `expr`.
#' @keywords internal
#' @noRd
with_local_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

#' Cross-validation of spherical kernel ridge regression
#'
#' Perform k-fold cross-validation of [sphere_regression()] over a grid
#' of kernel and penalty pairs `(kernels, lambdas)`, using the mean
#' out-of-fold geodesic prediction error
#' \deqn{\widehat{\mathrm{CV}}(k, \lambda) =
#'   \frac{1}{n_{\text{test}}}\sum_{i \in \text{test}}
#'   d_{\mathbb{S}}\bigl(\hat\mu(X_i), Y_i\bigr)}
#' averaged over folds as the selection criterion. With
#' `elastic = TRUE` the criterion becomes the elastic shape distance:
#' predictions and held-out responses are lifted back to their native
#' representation and aligned with [align_vectorized_curves()] before
#' the distance is taken.
#'
#' Candidate kernels are supplied as a *list of functions* of two
#' arguments rather than as a bandwidth grid, so that families other
#' than the plain Gaussian one can be screened; [create_kernel_list()]
#' builds the usual radial-basis (plus linear) grid.
#'
#' The reported `best_kernel` and `best_lambda` are the grid points
#' minimising the averaged criterion. Individual fits are wrapped in
#' `try(...)`, so a failure in one grid cell does not abort the sweep:
#' the corresponding entry of `mean_error` is averaged only over the
#' folds that succeeded, counted in `success`. Because a full sweep on
#' shape data can run for days, `partial_folds_fit` allows fitting only
#' a random subset of the folds, and `folds` allows reusing an
#' externally generated fold assignment across methods so that different
#' estimators are compared on identical splits.
#'
#' Note that the low-rank reductions `subK` and `subY` are recomputed
#' per fold from the training data, but the response coordinate system
#' `iso` is set up once from the *full* sample, so `subY` should be
#' understood as a computational device rather than as part of the
#' cross-validated pipeline.
#'
#' @param kernels list of candidate kernel functions, each with
#'   signature `(x1, x2)`; see [create_kernel_list()].
#' @param lambdas numeric vector of candidate penalties to evaluate.
#' @param folds_number number of folds (default 5). Ignored if `folds`
#'   is supplied.
#' @param random_folds logical; if `TRUE` (default), observations are
#'   assigned to folds at random, otherwise in blocks of consecutive
#'   indices.
#' @param partial_folds_fit optional number in \eqn{(0, 1]}: the share of
#'   folds actually fitted, drawn at random. Useful to keep the runtime
#'   of expensive elastic fits manageable. `NULL` (the default) fits all
#'   folds.
#' @param folds optional integer vector of length `nrow(X)` assigning
#'   each observation to a fold, overriding `folds_number` and
#'   `random_folds`.
#' @param align_eps,align_maxit convergence threshold and maximum number
#'   of iterations of the per-curve warping optimisation in the elastic
#'   *fitting* step (see [opt.gamj()]). With `"sd"`, relaxing them (e.g.
#'   `align_eps = 0.3`, `align_maxit = 10`) speeds the sweep up several
#'   fold at a small cost in alignment precision; with `"lbfgs"`,
#'   `align_maxit = 2` is the recommended shallow setting for the
#'   alternating fits.
#' @param align_eps_eval,align_maxit_eval alignment settings used when
#'   aligning held-out predictions for the CV *criterion*. Kept at the
#'   strict defaults independently of `align_eps`/`align_maxit`, so that
#'   error grids remain comparable across fitting settings; the
#'   evaluation step is a negligible share of the total cost.
#' @param cl optional cluster as returned by [make_align_cluster()] over
#'   which the per-curve alignments are distributed.
#' @param align_optimizer solver of the warping optimisation, `"lbfgs"`
#'   (default) or `"sd"`, used for both the fitting and the evaluation
#'   alignments; see [sphere_regression()].
#' @param ... further arguments passed on to the inner fitting routine
#'   built by [initialize_sphere_fit()].
#' @inheritParams sphere_regression
#'
#' @return A list with components
#'   \describe{
#'     \item{`mean_error`}{Matrix of average CV errors of size
#'       `length(kernels) x length(lambdas)`.}
#'     \item{`error_folds`}{List with one error matrix per fitted fold.}
#'     \item{`folds_to_fit`}{The folds that were actually fitted.}
#'     \item{`success`}{Matrix of the same shape as `mean_error`
#'       counting, for each grid cell, the number of folds in which the
#'       fit succeeded.}
#'     \item{`kernels`, `lambdas`, `folds_number`}{The input grids and
#'       fold count.}
#'     \item{`best_kernel`, `best_lambda`}{The kernel function and
#'       penalty minimising `mean_error`.}
#'     \item{`best_idx`}{Indices of the minimising grid cell(s).}
#'     \item{`time`}{Wall-clock time of the cross-validation loop.}
#'   }
#'
#' @seealso [sphere_regression()], [create_kernel_list()],
#'   [cross_val_tg_space_regression()] and [cross_val_vv_krr()] for the
#'   corresponding routines of the comparison methods.
#' @export
cross_val_sphere_regression <- function(X, Y,
                                        prod.u =  function (q1, q2) sum(q1*q2),
                                        kernels,
                                        lambdas = 10^seq(-5, -3, length.out = 10),
                                        folds_number = 5, fix_p= TRUE,
                                        subK = 1e-6, subY = if(nrow(Y)>ncol(Y)) 0,
                                        S = NULL, verbose = TRUE,
                                        control_optimizer = list(method = "BFGS"),
                                        start.xi = NULL, start.p = NULL,
                                        maxit = 100, eps_p = 1e-6, eps_xi = 1e-8, eps_diff = 1e-8,
                                        stp_p = .1, elastic = FALSE, random_folds = TRUE, partial_folds_fit = NULL, folds = NULL,
                                        align_eps = 0.05, align_maxit = 50,
                                        align_eps_eval = 0.05, align_maxit_eval = 50, cl = NULL,
                                        align_optimizer = c("lbfgs", "sd"), ...){
  align_optimizer <- match.arg(align_optimizer)

  coordinates <- !is.null(subY)

  # compute full design matrix
  if(!is.matrix(X)) X <- as.matrix(X)

  # sample folds
  n_total <- nrow(X)

  if (is.null(folds)){
    if (random_folds ) {
      folds <- sample(rep(1:folds_number, length.out = n_total))
    } else {
      folds = rep(1:folds_number, each = ceiling(n_total / folds_number))[1:n_total]
    }
  } else {
    folds_number = length(unique(folds))
  }

  # initialize the geometry
  if(is.null(S)) S = initialize_sphere(prod.u)

  # if working in coordinates, transform the input
  if(coordinates) {
    iso <- subspace_setup(Y, S$prod.u, subY)
    Y <- iso$coos
    if(!is.null(start.xi))
      start.xi <- iso$trafo(start.xi, back = FALSE)
    S_intern <- initialize_sphere()
    if(!is.null(start.p)) {
      start.p <- iso$trafo(start.p, back = FALSE)
      start.p <- start.p/S_intern$norm.u(start.p)
    }
    # renormalize (for potential low-rank approximation)
    Y <- sweep(Y, 2, apply(Y,2,S_intern$norm.u), `/`)
  } else {
    S_intern <- S
  }

  # initialize fitting function
  sphere.regr = initialize_sphere_fit(S_intern, S, fix_p = fix_p,
                                      control_optimizer = control_optimizer,
                                      maxit = maxit, eps_p = eps_p, eps_diff = eps_diff, eps_xi = eps_xi,
                                      stp_p = stp_p, elastic = elastic, coordinates = coordinates, iso = iso,
                                      align_eps = align_eps, align_maxit = align_maxit, cl = cl,
                                      align_optimizer = align_optimizer)

  prod.u = S_intern$prod.u; norm.u = S_intern$norm.u; log.map = S_intern$log.map; exp.map = S_intern$exp.map

  errors_grid <- success <- matrix(0, nrow = length(kernels), ncol = length(lambdas),
                        dimnames = list(paste0("kernel=", 1:length(kernels)), paste0("lambda=", lambdas)))

  start.xi_ <- start.xi
  start.p_ <- start.p

  error_folds = list()
  if (! is.null(partial_folds_fit) ) {
    # reproducible choice of the folds to fit, without resetting the caller's RNG stream
    folds_to_fit = with_local_seed(folds_number,
      sample(1:folds_number, replace = FALSE, size = (partial_folds_fit*folds_number)))
  } else {folds_to_fit = 1:folds_number}

  t1 = Sys.time()
  for (fold in folds_to_fit) {
    error_fold <- matrix(0, nrow = length(kernels), ncol = length(lambdas))
    if (verbose) cat('------------ fold :', fold, '------------', '\n')
    test_idx <- which(folds == fold)
    train_idx <- setdiff(1:n_total, test_idx)

    XY_train <- X[train_idx, , drop = FALSE]
    XY_test <- X[test_idx, , drop = FALSE]
    D_train <- Y[, train_idx, drop = FALSE]
    D_test <- Y[, test_idx, drop = FALSE]

    n_train <- ncol(D_train)
    if(is.null(start.p_)) start.p <- karcher_mean(D_train, S = S_intern)
    logY <- sapply(1:n_train, function(i) log.map(start.p, D_train[, i]))

    for (i in seq_along(kernels)) {
      k_function <- kernels[[i]]

      # Compute kernel matrices
      K_ <- sapply(1:n_train, function(a) sapply(1:n_train, function(b) k_function(XY_train[a, ], XY_train[b, ])))

      Z <- K <- NULL
      lowrank <- subK
      if(!is.null(lowrank)) {
        e <- eigen(K_, symmetric = TRUE)
        total <- sum(abs(e$values))
        if(lowrank<1) {
          r <- sum(e$values>0)
          e$values <- e$values[seq_len(r)]
          tau <- cumsum(e$values/sum(e$values))
          lowrank <- which.max(1-tau <= abs(lowrank))
        }
        e$vectors <- e$vectors[, seq_len(lowrank), drop = FALSE]
        e$values <- e$values[seq_len(lowrank), drop = FALSE]
        attr(lowrank, "share") <- sum(abs(e$values)) / total
        Z <- sweep(e$vectors, 2, sqrt(e$values), `/`)
        K <- sweep(e$vectors, 2, sqrt(e$values), `*`)
      }

      # implement sum-to-zero identifiability constraint when p is not fixed
      if (fix_p) {
        if(is.null(K)) K <- K_
      } else {
        if(is.null(Z)) { # we want to project the columns of K_ into the orthogonal of the constant 1
          C = colSums(K_) # corresponds to the inner product with the constant: C = 1^T %*% K in R^(1xp)
          if (all(abs(C) < 1e-15)) {
            Z <- NULL
            K <- K_
          } else {
            Z = MASS::Null(C)  # make sure it's orthogonal to the constant: Z contains a basis of the null space of C: CZ =0, Z in R^(px(p-1))
            lf <- lsfit(K_, rep(1, nrow(K_)), intercept = FALSE)
            if (all(abs(lf$residuals) > 1e-15)) {
              Z <- cbind(Z, t(C))
            }
            K = K_ %*% Z # constrain to sum to 0: 1^TK = 0
          }

        } else {
          C = colSums(K)
          if(all(abs(C) > 1e-15)) {
            Z_ <- MASS::Null(C)
            lf <- lsfit(K, rep(1, nrow(K)), intercept = FALSE)
            if(all(abs(lf$residuals) > 1e-15)) {
              Z_ <- cbind(Z_, t(C))
            }
            Z <- Z %*% Z_
            K <- K %*% Z_
          }
        }
      }


      for (j in seq_along(lambdas)) {
        run <- try({ #make sure the whole loop doesn't crash if it does for one combination
          lambda_st <- lambdas[j]
          if (verbose) cat('kernel:' , i, 'lambda: ', lambda_st, '\n')


          {
               start.xi = if(is.null(Z)) {
      logY %*% solve(K + lambda_st * diag(n_train))
    }  else {
        if(fix_p) { # i.e., if the only constraint matrix is due to low-rank approximation
          logY %*% scale(K, scale = e$values[seq_len(lowrank)] + lambda_st, center = FALSE)
        } else {
          logY %*% K %*% solve(crossprod(K) + lambda_st * crossprod(Z, K_) %*% Z)
        }
      }
          }

          if (verbose) cat('----------- model fit -----------', '\n')
          model_ <- sphere.regr(Xi = start.xi, P = start.p, K = K, Y = D_train, Z = Z, lambda = lambda_st, verbose = verbose,
                                maxit = maxit , eps_p = eps_p, aps_diff = eps_diff, eps_xi = eps_xi, stp_p = stp_p,elastic = elastic , ...)
          if (verbose) cat('-------------- end --------------', '\n')
          # only for predict function

          model_$X <- XY_train
          model_$k <- k_function
          model_$Z <- Z
          model_$S_intern <- S_intern
          model_$P_intern <- model_$P
          model_$Xi_intern <- model_$Xi
          model_$subYdim <- ncol(model_$Xi)
          if(coordinates) {
            model_$iso <- iso
          }

          # Prediction
          y_pred <- predict.sphere_model(model_, newdata = XY_test, coordinates = coordinates)

          if (elastic){
            D_test_original = iso$trafo(D_test, back = TRUE)
            y_pred_original = iso$trafo(y_pred, back = TRUE)
            y_pred_aligned = align_sample(D_test_original, y_pred_original, S = S,
                                          eps = align_eps_eval, maxit = align_maxit_eval, cl = cl,
                                          optimizer = align_optimizer)$Y_aligned
            err_test <- sapply(1:ncol(y_pred), function(b) S$norm.u(S$log.map(y_pred_aligned[, b], D_test_original[, b])))
          }else{
            err_test <- sapply(1:ncol(y_pred), function(b) norm.u(log.map(y_pred[, b], D_test[, b])))
          }

          mean_err <- mean(err_test)
        })
        if(!inherits(run, "try-error")) {
          error_fold[i,j] = mean_err
          errors_grid[i, j] <- errors_grid[i, j] + mean_err
          success[i, j] <- success[i, j] + 1
          cat(' mean error : ', mean_err, '\n')
        }
      }
    }
    error_folds[[fold]] = error_fold
  }

  # Average over the folds
  errors_grid <- errors_grid / success
  t2 = Sys.time()
  cat('time:', t2 - t1, '\n')


  # Get best parameters
  min_err_idx <- which(errors_grid == min(errors_grid, na.rm = TRUE), arr.ind = TRUE)
  best_idx = min_err_idx
  if (is.matrix(best_idx)) best_idx = best_idx[1,]
  best_kernel <- kernels[[best_idx[1]]]
  best_lambda <- lambdas[best_idx[2]]
  best_indices <- min_err_idx

  dimnames(errors_grid) <- dimnames(success) <- list(kernels = seq_along(kernels), lambdas = seq_along(lambdas))

  cat("Best kernel:", min_err_idx[1], "\n")
  cat("Best lambda:", best_lambda, "\n")
  cat("Minimum average CV error:", errors_grid[min_err_idx], "\n")

  # return computations
  list(
    mean_error = errors_grid,
    error_folds = error_folds,
    folds_to_fit = folds_to_fit,
    success = success,
    kernels = kernels, lambdas = lambdas, folds_number = folds_number,
    best_kernel = best_kernel, best_lambda = best_lambda, time = t2-t1, best_idx = best_indices
  )
}



#' Build a grid of candidate kernels for cross-validation
#'
#' Construct the list of kernel functions expected by the `kernels`
#' argument of [cross_val_sphere_regression()],
#' [cross_val_tg_space_regression()] and [cross_val_vv_krr()]. Each
#' entry is a Gaussian radial basis function plus an optional linear
#' term,
#' \deqn{k(x, x') = \exp\!\left(-\frac{\|x - x'\|^2}{2\sigma^2}\right)
#'   + \beta\, \langle x, x'\rangle,}
#' evaluated over all combinations of `sigmas` and `betas`. The linear
#' component lets the fitted function pick up a global trend that a
#' purely local kernel would have to approximate with a large bandwidth;
#' the default `betas = 0` switches it off and gives a plain
#' radial-basis grid.
#'
#' Bandwidths are frozen into each closure with [base::force()], so that
#' the returned functions take only two arguments and can be passed
#' around --- and stored in a fitted model --- without depending on the
#' loop variables they were created from.
#'
#' @param sigmas numeric vector of squared bandwidths \eqn{\sigma^2}.
#' @param betas numeric vector of weights for the linear kernel
#'   component.
#'
#' @return A list of `length(sigmas) * length(betas)` functions with
#'   signature `(x1, x2)`, varying `betas` fastest.
#'
#' @seealso [cross_val_sphere_regression()].
#' @export
#' @examples
#' kernels <- create_kernel_list(sigmas = c(1, 10), betas = 0)
#' kernels[[1]](c(0, 0), c(1, 1))
create_kernel_list <- function(sigmas = seq(1, 5, length.out = 4),
                               betas   = c(0)) {
  kernels <- list()

  # helper: create a single kernel and freeze its parameters
  make_rbf_plus_linear <- function(sigma2,  beta) {
    force(sigma2); force(beta)  # ensure eager binding
    function(x1, x2) {
      diff <- x1 - x2
      k_rbf    <- exp(-sum(diff^2) / (2 * sigma2))
      k_linear <- sum(x1 * x2)
      k_rbf + beta * k_linear
    }
  }


  rbf_plus_linear_kernels <- list()
  for (sigma2 in sigmas) {
    for (b in betas) {
      rbf_plus_linear_kernels[[length(rbf_plus_linear_kernels) + 1]] <-
        make_rbf_plus_linear(sigma2, b)
    }
  }

  kernels <- c(kernels, rbf_plus_linear_kernels)


  kernels
}



#######################################################
############# ELASTIC METHODS: COMPARISON #############
#######################################################




#######################################
########## Vector-valued KRR ##########
#######################################


#' Vector-valued kernel ridge regression (unconstrained baseline)
#'
#' Fit ordinary kernel ridge regression of the responses on the
#' covariates, ignoring the spherical geometry entirely: the responses
#' are treated as elements of \eqn{\mathbb{R}^m} and the coefficients
#' solve the linear system
#' \deqn{(K + \lambda I_n)\, A = Y^\top,}
#' so that \eqn{\hat f(x) = \sum_i k(X_i, x)\, A_{i\cdot}}.
#'
#' This is the naive comparison method in the simulation study: its
#' predictions do not in general lie on the sphere --- for shape data
#' they are not unit-norm square-root velocity functions and hence do not
#' correspond to curves of the right scale --- and its error is measured
#' in the ambient \eqn{L^2} metric rather than the geodesic one. It
#' quantifies how much is gained by respecting the geometry in
#' [sphere_regression()], and how much of that gain survives when only
#' the parametrisation is additionally quotiented out, as in
#' [tg_space_regression_model()].
#'
#' The system is solved with [base::solve()] on the matrix right-hand
#' side rather than by inverting \eqn{K + \lambda I} explicitly, which
#' is both faster and numerically better behaved.
#'
#' @param X_train covariate matrix of size \eqn{n \times d}, one
#'   observation per row.
#' @param y_train response matrix of size \eqn{m \times n}, one
#'   observation per column (the same orientation as elsewhere in the
#'   package).
#' @param lambda ridge penalty.
#' @param k a kernel function with signature `(x1, x2)`, or an
#'   already-evaluated \eqn{n \times n} kernel matrix.
#'
#' @return An object of class `"vv_krr_model"`: a list with the
#'   coefficient matrix `A` (of size \eqn{n \times m}), the training
#'   covariates `X_train`, the penalty `lambda`, and the kernel `k`.
#'
#' @seealso [predict.vv_krr_model()], [cross_val_vv_krr()],
#'   [sphere_regression()].
#' @export
vv_krr_train <- function(X_train, y_train, lambda = 1e-3, k = function(x1, x2, sigma2 = 1) exp(-sum((x1 - x2)^2) / (2 * sigma2))  ) {
  n <- nrow(X_train)

  # scalar kernel Gram matrix
  if (is.matrix(k)) K = k else K = sapply(1:n, function(i) sapply(1:n, function(j) k(X_train[i,], X_train[j,])))

  # solve (K + n*lambda*I_n) A = Y  for A (n x q)
  # using solve with RHS matrix is more stable than inverting explicitly
  A <- solve(K +  lambda * diag(n), t(y_train))

  model = list(A = A,
       X_train = X_train,
       lambda = lambda,
        k = k)
  class(model) = c('vv_krr_model', class(model))
  return(model)
}

#' Predict from a vector-valued kernel ridge regression fit
#'
#' Evaluate a [vv_krr_train()] model at new covariates,
#' \eqn{\hat f(x) = \sum_i k(X_i, x)\, A_{i\cdot}}. The result is
#' returned with one prediction per column, matching the response
#' orientation used throughout the package, whereas the coefficient
#' matrix `A` is stored the other way round.
#'
#' No projection onto the sphere is applied: predictions of this
#' baseline model generally have norm different from one, which is
#' exactly the deficiency it is included to illustrate.
#'
#' @param object a fitted model of class `vv_krr_model`.
#' @param newdata matrix of new covariates with the same number of columns
#'   as `model$X_train`.
#'
#' @param ... further arguments, currently ignored.
#' @return A matrix of size \eqn{m \times \mathrm{nrow}(X_{new})} with
#'   one prediction per column.
#'
#' @seealso [vv_krr_train()], [cross_val_vv_krr()].
#' @export
predict.vv_krr_model <- function(object, newdata, ...) {
  model <- object; X_new <- newdata

  K_new = sapply(1:nrow(model$X_train), function(i) sapply(1:nrow(X_new), function(j) model$k(model$X_train[i,], X_new[j,])))

  # f(X_new) = K_new %*% A  (m x n) (n x q) = (m x q)
  Y_hat <- K_new %*% model$A
  t(Y_hat)
}






#' Cross-validation of vector-valued kernel ridge regression
#'
#' k-fold cross-validation of [vv_krr_train()] over a grid of kernels
#' and penalties, structured exactly like
#' [cross_val_sphere_regression()] so that the two can be run on the
#' same fold assignment (pass the same `folds` to both) and their
#' selected models compared on equal footing.
#'
#' The selection criterion is the root mean squared error in the ambient
#' space, \eqn{\sqrt{m^{-1}\|\hat y_i - y_i\|^2}} averaged over the
#' held-out observations --- not a geodesic distance, since predictions
#' of this model need not lie on the sphere. Errors are therefore *not*
#' directly comparable in absolute terms to those reported by the
#' spherical and tangent-space routines.
#'
#' As there, each fit is wrapped in `try(...)` so a failure in one grid
#' cell does not abort the sweep.
#'
#' @param X covariate matrix of size \eqn{n \times d}, one observation
#'   per row.
#' @param Y response matrix of size \eqn{m \times n}, one observation per
#'   column.
#' @param kernels list of candidate kernel functions with signature
#'   `(x1, x2)`; see [create_kernel_list()].
#' @param lambdas numeric vector of candidate ridge penalties.
#' @param partial_folds_fit optional number in \eqn{(0, 1]}: the share of
#'   folds actually fitted, drawn at random. `NULL` (the default) fits
#'   all folds.
#' @param folds_number number of folds (default 5). Ignored if `folds`
#'   is supplied.
#' @param verbose logical; if `TRUE`, print progress information.
#' @param random_folds logical; if `TRUE` (default), observations are
#'   assigned to folds at random, otherwise in blocks of consecutive
#'   indices.
#' @param folds optional integer vector of length `nrow(X)` assigning
#'   each observation to a fold, overriding `folds_number` and
#'   `random_folds`.
#' @param ... further arguments, currently ignored.
#'
#' @return A list with components `error` (matrix of average CV errors),
#'   `success`, the inputs `kernels`, `lambdas` and `folds_number`, the
#'   selected `best_kernel` and `best_lambda`, the minimising `best_idx`,
#'   and the elapsed `time`.
#'
#' @seealso [vv_krr_train()], [cross_val_sphere_regression()].
#' @export
cross_val_vv_krr<- function(X, Y, kernels,
                                  lambdas = 10^seq(-5, -3, length.out = 10), partial_folds_fit = NULL,
                                  folds_number = 5, verbose = TRUE, random_folds = TRUE, folds = NULL, ...){


  # compute full design matrix
  if(!is.matrix(X)) X <- as.matrix(X)
  n_total <- nrow(X)

  # sample folds

  if (is.null(folds)){
    if (random_folds ) {
      folds <- sample(rep(1:folds_number, length.out = n_total))
    } else {
      folds = rep(1:folds_number, each = ceiling(n_total / folds_number))[1:n_total]
    }
  } else {
    folds_number = length(unique(folds))
  }

  # initialize the geometry
  errors_grid <- success <- matrix(0, nrow = length(kernels), ncol = length(lambdas),
                                   dimnames = list(paste0("kernels=", 1:length(kernels)), paste0("lambda=", lambdas)))

  if (! is.null(partial_folds_fit) ) {
    folds_to_fit = sample(1:folds_number, replace = FALSE, size = (partial_folds_fit*folds_number))
  } else {folds_to_fit = 1:folds_number}


  t1 = Sys.time()

  for (fold in folds_to_fit) {
    if (verbose) cat('------------ fold :', fold, '------------', '\n')
    test_idx <- which(folds == fold)
    train_idx <- setdiff(1:n_total, test_idx)

    XY_train <- X[train_idx, , drop = FALSE]
    XY_test <- X[test_idx, , drop = FALSE]
    D_train <- Y[, train_idx, drop = FALSE]
    D_test <- Y[, test_idx, drop = FALSE]

    n_train <- ncol(D_train)

    for (i in seq_along(kernels)) {
      k_function <- kernels[[i]]


      for (j in seq_along(lambdas)) {
        run <- try({ #make sure the whole loop doesn't crash if it does for one combination
          lambda_st <- lambdas[j]
          if (verbose) cat('kernel:' , i, 'lambda: ', lambda_st, '\n')


          if (verbose) cat('----------- model fit -----------', '\n')


          model_ <- vv_krr_train(X_train = XY_train, y_train = D_train,
                                 lambda = lambda_st,  k = k_function  )


          if (verbose) cat('-------------- end --------------', '\n')
          # only for predict function

          model_$X <- XY_train
          model_$k <- k_function


          # Prediction
          y_pred <- predict.vv_krr_model(model_, XY_test)

          err_test <- sqrt(colMeans((y_pred - D_test)^2))


          mean_err <- mean(err_test)
        })
        if(!inherits(run, "try-error")) {
          errors_grid[i, j] <- errors_grid[i, j] + mean_err
          success[i, j] <- success[i, j] + 1
          cat(' mean error : ', mean_err, '\n')
        }
      }
    }
  }


  # Average over the folds
  errors_grid <- errors_grid / success
  t2 = Sys.time()
  cat('time:', t2 - t1, '\n')


  # Get best parameters
  min_err_idx <- which(errors_grid == min(errors_grid, na.rm = TRUE), arr.ind = TRUE)
  best_idx = min_err_idx
  if (is.matrix(best_idx)) best_idx = best_idx[1,]
      best_kernel <- kernels[[best_idx[1]]]
      best_lambda <- lambdas[best_idx[2]]
  best_kernel <- kernels[[min_err_idx[1]]]
  best_lambda <- lambdas[min_err_idx[2]]
  best_indices <- min_err_idx

  dimnames(errors_grid) <- dimnames(success) <- list(kernels = seq_along(kernels), lambdas = seq_along(lambdas))

  cat("Best kernel:", min_err_idx[1], "\n")
  cat("Best lambda:", best_lambda, "\n")
  cat("Minimum average CV error:", errors_grid[min_err_idx], "\n")

  # return computations
  list(
    error = errors_grid,
    success = success,
    kernels = kernels, lambdas = lambdas, folds_number = folds_number,
    best_kernel = best_kernel, best_lambda = best_lambda, time = t2-t1, best_idx = best_indices
  )
}










######################################################
########## (elastic) tangent space regression ##########
######################################################



#' Tangent-space kernel ridge regression (linearised comparison)
#'
#' Fit a kernel ridge regression in a *single* tangent space instead of
#' intrinsically on the sphere. The responses are mapped once to
#' \eqn{T_o\mathbb{S}} via \eqn{v_i = \mathrm{Log}_o(Y_i)} at the Karcher
#' mean \eqn{o}, an ordinary kernel ridge regression is fitted there in
#' closed form,
#' \deqn{\hat\Xi = V (K + \lambda I_n)^{-1},}
#' and predictions are pushed back with \eqn{\mathrm{Exp}_o}. This is the
#' standard linearisation used in shape statistics, and the natural
#' comparison for [sphere_regression()]: the two coincide to first order
#' near \eqn{o} and differ increasingly as the responses spread out over
#' the sphere.
#'
#' With `elastic = TRUE` the fit is alternated with an elastic
#' re-alignment of the responses, exactly as in the intrinsic case: at
#' each iteration the responses are aligned to the current predictions
#' with [align_vectorized_curves()], re-projected to the tangent space,
#' and the closed-form solution recomputed. Iteration stops once the
#' largest change of any fitted value falls below `eps_`, or after
#' `maxit` rounds. Note that the linearisation point \eqn{o} stays fixed
#' at the Karcher mean of the *original* responses and is not updated
#' along with the alignment.
#'
#' With `maxit = 0` or `elastic = FALSE` the closed-form fit is returned
#' directly, and the model carries `iterations = 0`. The intercept
#' cannot be optimised over: `fix_p = FALSE` prints a message and returns
#' the unrefined fit.
#'
#' Unlike [sphere_regression()], this routine offers no low-rank
#' reduction of either the responses or the kernel matrix, so it operates
#' throughout in the responses' native representation.
#'
#' @param X covariate matrix of size \eqn{n \times d}, one observation
#'   per row.
#' @param Y response matrix of size \eqn{m \times n} with one sphere
#'   point per column.
#' @param k a kernel function with signature `(x1, x2)`, or an
#'   already-evaluated \eqn{n \times n} kernel matrix.
#' @param prod.u inner product defining the response geometry, passed to
#'   [initialize_sphere()] when `S` is `NULL`.
#' @param lambda ridge penalty.
#' @param fix_p logical; must be `TRUE`, as optimisation over the
#'   linearisation point is not implemented here.
#' @param S optional sphere geometry as returned by
#'   [initialize_sphere()].
#' @param verbose logical; if `TRUE`, print progress information.
#' @param start.xi optional starting value for the coefficient matrix;
#'   if `NULL`, the closed-form ridge solution is used.
#' @param start.p optional linearisation point; if `NULL`, the Karcher
#'   mean of `Y` is used.
#' @param maxit maximum number of alignment/fitting alternations for
#'   elastic fits; `0` returns the closed-form fit.
#' @param eps_ convergence threshold on the largest change of any fitted
#'   value between two successive iterations.
#' @param elastic logical; if `TRUE`, alternate the fit with an elastic
#'   re-alignment of the responses.
#' @param align_eps,align_maxit convergence threshold and maximum number
#'   of iterations of the per-curve warping optimisation in the elastic
#'   alignment step (see [opt.gamj()]); `align_eps` only concerns the
#'   `"sd"` solver, `align_maxit` caps either. The alternation of this
#'   model is only stable with a shallow alignment: with `"lbfgs"` use a
#'   small cap such as `align_maxit = 2`, since a fully converged warping
#'   makes the fit drift away from the data over the iterations.
#' @param align_optimizer solver of the per-curve warping optimisation:
#'   `"lbfgs"` (default, L-BFGS-B) or `"sd"` (Riemannian steepest
#'   descent), see [opt.gamj()].
#' @param cl optional cluster as returned by [make_align_cluster()] over
#'   which the per-curve alignments are distributed.
#'
#' @return An object of class `"tg_model"`: a list with the coefficient
#'   matrix `Xi`, the linearisation point `P`, the number of
#'   `iterations` performed, the kernel matrix `K` and kernel `k`, the
#'   geometry `S`, the penalty `lambda`, the design matrix `X`, the
#'   untouched responses `Y_original`, and --- for elastic fits --- the
#'   re-aligned responses `Y_elastic_align` together with the warpings
#'   `gamma_alignments` of each iteration.
#'
#' @seealso [predict.tg_model()], [cross_val_tg_space_regression()],
#'   [sphere_regression()] for the intrinsic counterpart.
#' @export
tg_space_regression_model = function (X, Y,
                                k = function(x1, x2, sigma2 = 1) exp(-sum((x1 - x2)^2) / (2 * sigma2)),
                                prod.u =  function (q1, q2) sum(q1*q2),
                                lambda = 1e-3, fix_p= TRUE,
                                S = NULL,
                                verbose = TRUE, start.xi = NULL, start.p = NULL,
                                maxit = 0, eps_ = 1e-6,
                                elastic = FALSE, align_eps = 0.05, align_maxit = 50, cl = NULL,
                                align_optimizer = c("lbfgs", "sd")){
  align_optimizer <- match.arg(align_optimizer)



  n <- nrow(X)
  if (is.matrix(k)) K = k else K = sapply(1:n, function(i) sapply(1:n, function(j) k(X[i,], X[j,])))
  if(is.null(S))  S = initialize_sphere(prod.u = prod.u)
  prod.u = S$prod.u; norm.u = S$norm.u; log.map = S$log.map; exp.map = S$exp.map

  # get start.p and start.xi
  if ( is.null(start.p)) start.p <- karcher_mean(Y, S = S)
  logY = sapply(1:n, function(i) S$log.map(start.p, Y[, i]))

  if(is.null(start.xi)) { # COMMENT: fastest version of the below, using positive definiteness
    R <- chol(K + lambda * diag(n))           # A = t(R) %*% R, R upper triangular
    start.xi <- t(backsolve(R, forwardsolve(t(R), t(logY))))
  }
    # { # COMMENT: this is the simplified version of the below, working only for squared design matrices
    #   start.xi = logY %*% solve(K + lambda * diag(n)) # a faster version would be: t(solve(K + lambda * diag(n), t(logY)))
    # }
    # { # COMMENT: this would be the kernel ridge regression version of the below, with penalty matrix K implementing the quadratic norm penalty
    #   start.xi = logY %*% K %*% solve(crossprod(K) + lambda * K)
    # }
  # { # COMMENT: this implements the solution of PLS with design matrix K and ridge penalty directly on the coefficients
  #   start.xi =  t(solve(t(K) %*% K + lambda * diag(n)) %*% t(K) %*% t(logY))
  # }

  gamma_alignments = list()
  Y_elastic_align =  NA

  model = list(Xi = start.xi, P = start.p, iterations = 0 ,
               K = K, k = k, S = S,  lambda = lambda, X = X, Y_original = Y,
               Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments)
  class(model) = c('tg_model', class(model))

  if (!fix_p) {
    cat('This method is not implemented for optimization over the origin')
    return(model)
  }


  if (maxit == 0 | elastic == FALSE){
    return(model)
  }

  Xi = start.xi


  if (elastic){
    for (it in 1:maxit){

        # align srvf
      fx_pred <- Xi %*% t(K)
      y_pred <- sapply(1:ncol(fx_pred), function(b) exp.map(start.p, fx_pred[, b]))

      aligned = align_sample(y_pred, Y, S = S, eps = align_eps, maxit = align_maxit, cl = cl, optimizer = align_optimizer)
      Y_aligned = aligned$Y_aligned
      gamma_alignments[[it]] = aligned$gammas
      Y_elastic_align = Y_aligned



      logY = sapply(1:n, function(i) log.map(start.p, Y_aligned[, i]))
      # Xi0 = Xi
      Xi =  t(solve(t(K) %*% K + lambda * diag(n)) %*% t(K) %*% t(logY))
      fx0 <- Xi %*% t(K)

      max_norm_diff = max (sapply(1:ncol(fx0), function (i) norm.u( fx0[,i]- fx_pred[,i] )))
      model = list(Xi = Xi, P = start.p, iterations = it ,
                   K = K, k = k, S = S,  lambda = lambda,  X = X, Y_original = Y,
                   Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments)

      class(model) = c('tg_model', class(model))
      if (max_norm_diff <  eps_) return(model)

    }
    model = list(Xi = Xi, P = start.p, iterations = it ,
                 K = K, k = k, S = S,  lambda = lambda, X = X, Y_original = Y,
                 Y_elastic_align = Y_elastic_align, gamma_alignments = gamma_alignments)
    class(model) = c('tg_model', class(model))
    return(model)
  }
}



#' Predict from a tangent-space regression fit
#'
#' Evaluate a [tg_space_regression_model()] fit at the training
#' covariates or at new ones, mapping the tangent-space prediction back
#' to the sphere with \eqn{\mathrm{Exp}_o}. Unlike
#' [predict.sphere_model()] there is no `type` argument: only the
#' response scale is returned, since the model carries no coordinate
#' reduction and the linear predictor is simply `object$Xi %*% t(K)`.
#'
#' As for the intrinsic model, `elastic = TRUE` aligns each prediction to
#' the corresponding observed curve and therefore only applies to
#' in-sample predictions; it is silently ignored when `newdata` is
#' given.
#'
#' @param object a fitted model of class `tg_model`.
#' @param newdata optional matrix of new covariates with the same number
#'   of columns as `object$X`. If `NULL` (the default), the stored
#'   kernel matrix is reused.
#' @param elastic logical; if `TRUE` and `newdata` is `NULL`, align each
#'   prediction to the corresponding observed curve with
#'   [align_vectorized_curves()].
#' @param ... further arguments, currently ignored.
#'
#' @return A matrix with one predicted sphere point per column.
#'
#' @seealso [tg_space_regression_model()], [predict.sphere_model()].
#' @export
predict.tg_model <- function (object, newdata = NULL, elastic = FALSE, ...){

  if (is.null(newdata)){
    fx = object$Xi %*% t(object$K)
  }else{
    if(is.null(ncol(newdata)))
      newdata <- as.matrix(newdata)
    stopifnot(ncol(newdata) == ncol(object$X))

    k_test = sapply(1:nrow(object$X), function(i) sapply(1:nrow(newdata), function(j) object$k(object$X[i,], newdata[j,])))

    fx = object$Xi%*% t(k_test)
  }
  y <- sapply(1:ncol(fx), function (j) object$S$exp.map( object$P,  fx [,j]))

  if (elastic) {
    if (is.null(newdata)){
      y = sapply ( 1: ncol(y), function (j) align_vectorized_curves (object$Y_original[,j] ,y[,j] , object$S)$aligned_q2)
    }
  }
  y
}




#' Cross-validation of tangent-space regression
#'
#' k-fold cross-validation of [tg_space_regression_model()] over a grid
#' of kernels and penalties, mirroring [cross_val_sphere_regression()]
#' so that the intrinsic and the linearised estimator can be compared on
#' the same folds (pass the same `folds` to both).
#'
#' The selection criterion is the mean out-of-fold geodesic distance
#' between prediction and held-out response, and with `elastic = TRUE`
#' the elastic shape distance obtained after aligning the two with
#' [align_vectorized_curves()]. Since this model works in the responses'
#' native representation throughout, no coordinate round trip is needed
#' before aligning --- in contrast to [cross_val_sphere_regression()].
#'
#' The linearisation point is recomputed per fold as the Karcher mean of
#' the training responses, unless `start.p` is supplied.
#'
#' @param X covariate matrix of size \eqn{n \times d}, one observation
#'   per row.
#' @param Y response matrix of size \eqn{m \times n} with one sphere
#'   point per column.
#' @param kernels list of candidate kernel functions with signature
#'   `(x1, x2)`; see [create_kernel_list()].
#' @param lambdas numeric vector of candidate ridge penalties.
#' @param folds_number number of folds (default 5). Ignored if `folds`
#'   is supplied.
#' @param random_folds logical; if `TRUE` (default), observations are
#'   assigned to folds at random, otherwise in blocks of consecutive
#'   indices.
#' @param partial_folds_fit optional number in \eqn{(0, 1]}: the share of
#'   folds actually fitted, drawn at random. `NULL` (the default) fits
#'   all folds.
#' @param folds optional integer vector of length `nrow(X)` assigning
#'   each observation to a fold, overriding `folds_number` and
#'   `random_folds`.
#' @param align_eps_eval,align_maxit_eval alignment settings used when
#'   aligning held-out predictions for the CV criterion; kept at the
#'   strict defaults independently of the fitting settings
#'   `align_eps`/`align_maxit`.
#' @param align_optimizer solver of the warping optimisation, `"lbfgs"`
#'   (default) or `"sd"`, used for both the fitting and the evaluation
#'   alignments; see [tg_space_regression_model()].
#' @param ... further arguments, currently ignored.
#' @inheritParams tg_space_regression_model
#'
#' @return A list with components `error` (matrix of average CV errors),
#'   `error_folds`, `folds_to_fit`, `success`, the inputs `kernels`,
#'   `lambdas` and `folds_number`, the selected `best_kernel` and
#'   `best_lambda`, the minimising `best_idx`, and the elapsed `time`.
#'
#' @seealso [tg_space_regression_model()],
#'   [cross_val_sphere_regression()], [cross_val_vv_krr()].
#' @export
cross_val_tg_space_regression <- function(X, Y,
                                        prod.u =  function (q1, q2) sum(q1*q2),
                                        kernels,
                                        lambdas = 10^seq(-5, -3, length.out = 10),
                                        folds_number = 5, fix_p= TRUE,
                                        S = NULL, verbose = TRUE,
                                        start.xi = NULL, start.p = NULL,
                                        maxit = 20, eps_ = 1e-5, elastic = FALSE, random_folds = TRUE,
                                        partial_folds_fit = NULL, folds = NULL,
                                        align_eps = 0.05, align_maxit = 50,
                                        align_eps_eval = 0.05, align_maxit_eval = 50, cl = NULL,
                                        align_optimizer = c("lbfgs", "sd"), ...){
  align_optimizer <- match.arg(align_optimizer)

  # compute full design matrix
  if(!is.matrix(X)) X <- as.matrix(X)

  # sample folds
  n_total <- nrow(X)

  if (is.null(folds)){
    if (random_folds ) {
      folds <- sample(rep(1:folds_number, length.out = n_total))
    } else {
      folds = rep(1:folds_number, each = ceiling(n_total / folds_number))[1:n_total]
    }
  } else {
    folds_number = length(unique(folds))
  }


  # initialize the geometry
  if(is.null(S)) S = initialize_sphere(prod.u)

  prod.u = S$prod.u; norm.u = S$norm.u; log.map = S$log.map; exp.map = S$exp.map

  errors_grid <- success <- matrix(0, nrow = length(kernels), ncol = length(lambdas),
                                   dimnames = list(paste0("kernels=", 1:length(kernels)), paste0("lambda=", lambdas)))

  start.xi_ <- start.xi
  start.p_ <- start.p

  error_folds = list()
  if (! is.null(partial_folds_fit) ) {
    # reproducible choice of the folds to fit, without resetting the caller's RNG stream
    folds_to_fit = with_local_seed(folds_number,
      sample(1:folds_number, replace = FALSE, size = (partial_folds_fit*folds_number)))
  } else {folds_to_fit = 1:folds_number}


  t1 = Sys.time()

  for (fold in folds_to_fit) {
    error_fold <-  matrix(0, nrow = length(kernels), ncol = length(lambdas))
    if (verbose) cat('------------ fold :', fold, '------------', '\n')
    test_idx <- which(folds == fold)
    train_idx <- setdiff(1:n_total, test_idx)

    XY_train <- X[train_idx, , drop = FALSE]
    XY_test <- X[test_idx, , drop = FALSE]
    D_train <- Y[, train_idx, drop = FALSE]
    D_test <- Y[, test_idx, drop = FALSE]

    n_train <- ncol(D_train)
    if(is.null(start.p_)) start.p <- karcher_mean(D_train, S = S)
    logY <- sapply(1:n_train, function(i) log.map(start.p, D_train[, i]))

    for (i in seq_along(kernels)) {
      k_function <- kernels[[i]]


      for (j in seq_along(lambdas)) {
        run <- try({ #make sure the whole loop doesn't crash if it does for one combination
          lambda_st <- lambdas[j]
          if (verbose) cat('kernel number: ' , i, 'lambda: ', lambda_st, '\n')


          if (verbose) cat('----------- model fit -----------', '\n')

          model_ <- tg_space_regression_model (X = XY_train, Y = D_train,
                                               k = k_function,
                                               prod.u =  prod.u , lambda = lambda_st, fix_p= fix_p,
                                               S = S,
                                               verbose = TRUE, start.xi = NULL, start.p = start.p ,
                                               maxit =maxit, eps_ = eps_,
                                               elastic = elastic,
                                               align_eps = align_eps, align_maxit = align_maxit, cl = cl,
                                      align_optimizer = align_optimizer)

          if (verbose) cat('-------------- end --------------', '\n')
          # only for predict function

          model_$X <- XY_train
          model_$k <- k_function


          # Prediction
          y_pred <- predict.tg_model(object = model_, newdata = XY_test, elastic = elastic)

          if (elastic){
            y_pred_aligned = align_sample(D_test, y_pred, S = S,
                                          eps = align_eps_eval, maxit = align_maxit_eval, cl = cl,
                                          optimizer = align_optimizer)$Y_aligned
            err_test <- sapply(1:ncol(y_pred), function(b) norm.u(log.map(y_pred_aligned[, b], D_test[, b])))
          }else{
            err_test <- sapply(1:ncol(y_pred), function(b) norm.u(log.map(y_pred[, b], D_test[, b])))
          }

          mean_err <- mean(err_test)
        })
        if(!inherits(run, "try-error")) {

          error_fold [i,j] = mean_err
          errors_grid[i, j] <- errors_grid[i, j] + mean_err
          success[i, j] <- success[i, j] + 1

          cat(' mean error : ', mean_err, '\n')
        }
      }
    }
    error_folds [[fold]] = error_fold
  }

  # Average over the folds
  errors_grid <- errors_grid / success
  t2 = Sys.time()
  cat('time:', t2 - t1, '\n')


  # Get best parameters
  min_err_idx <- which(errors_grid == min(errors_grid, na.rm = TRUE), arr.ind = TRUE)
  best_idx = min_err_idx
  if (is.matrix(best_idx)) best_idx = best_idx[1,]
      best_kernel <- kernels[[best_idx[1]]]
      best_lambda <- lambdas[best_idx[2]]
  best_kernel <- kernels[[min_err_idx[1]]]
  best_lambda <- lambdas[min_err_idx[2]]
  best_indices <- min_err_idx

  dimnames(errors_grid) <- dimnames(success) <- list(kernels = seq_along(kernels), lambdas = seq_along(lambdas))

  cat("Best kernel:", min_err_idx[1], "\n")
  cat("Best lambda:", best_lambda, "\n")
  cat("Minimum average CV error:", errors_grid[min_err_idx], "\n")

  # return computations
  list(
    error = errors_grid,
    error_folds = error_folds,
    folds_to_fit = folds_to_fit,
    success = success,
    kernels = kernels, lambdas = lambdas, folds_number = folds_number,
    best_kernel = best_kernel, best_lambda = best_lambda, time = t2-t1, best_idx = best_indices
  )
}







