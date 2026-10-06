# sphereg2

Kernel ridge regression for responses on the unit sphere of a Hilbert space,
with an elastic variant for planar shapes represented by their square-root
velocity functions (SRVFs). The package builds on `sphereg`, the
implementation of the spherical kernel ridge regression of Matteo, Stöcker and
Tavakoli (2026, arXiv:2606.00181), and adds the elastic shape regression of

> Matteo, B., Stöcker, A., Pigoli, D. and Tavakoli, S. *Elastic kernel ridge regression, with applications in phonetics* (submitted).

The regression model, its estimation and the alignment algorithms are
described in the paper and in its Web Appendix; the package documentation
(`?sphere_regression`, `?opt.gamj`, ...) describes the implementation. The
scripts and the derived data reproducing the simulations, the applications and
the figures of the paper are distributed with it as supplementary material.

## Installation

```r
devtools::install_github("bmatteo98/sphereg2")
```

or, from the package tarball distributed as supplementary material,

```r
install.packages("sphereg2_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Contents

- `R/source_sphere.R`: the geometry of the sphere (`initialize_sphere()`
  returns the exponential and logarithmic maps and the parallel transport for
  a user-supplied inner product), the spherical kernel ridge regression
  `sphere_regression()` with its `predict` method and the cross-validation
  `cross_val_sphere_regression()`, and two comparison models, the
  tangent-space regression `tg_space_regression_model()` and the
  vector-valued kernel ridge regression `vv_krr_train()`, each with its
  cross-validation. Responses are stored one per column; a planar SRVF sampled
  at `m` points is a vector of length `2m`, the `x` coordinates followed by
  the `y` coordinates.
- `R/source_srvf.R`: the elastic shape machinery: optimal rotation `opt.O()`,
  optimal reparametrization `opt.gamj()` with the three solvers selected by
  its `optimizer` argument (`"lbfgs"`, the default, L-BFGS-B on the
  scale-invariant extension of the cost; `"sd"`, Riemannian steepest descent;
  `"rlbfgs"`, Riemannian limited-memory BFGS on the sphere of warpings),
  the alignment of one curve to another (`align_curves()`,
  `align_vectorized_curves()`) or of a sample to its Karcher mean
  (`align_multiple_curves()`), and `make_align_cluster()` to spread the
  alignments of an elastic fit over several workers.
- `R/contour_ordering.R`: the pre-processing of segmented images used in the
  applications of the paper: `order.contour()` orders the pixels labelled as
  one object along its contour (diameter of the Euclidean minimum spanning
  tree), `orient_segments()` and `bind_with_connector()` join several
  contours into one open curve. Not exported.

## Notes

- Every elastic fit alternates between fitting the regression and re-aligning
  the responses to the current fit; `align_maxit` caps the warping solver
  inside the alternation, and the cross-validation functions distinguish these
  settings from those used to score the held-out curves (`align_maxit_eval`).
- Cross-validating an elastic fit costs several times a non-elastic one, and
  the cost of a fit varies strongly over the grid of hyper-parameters.
