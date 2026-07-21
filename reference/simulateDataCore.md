# Generate Random Panel Data for FETWFE Simulations (core)

Generates a random panel data set for simulation studies of the fused
extended two-way fixed effects (FETWFE) estimator. The function creates
a balanced panel with \\N\\ units over \\T\\ time periods, assigns
treatment status across \\G\\ treated cohorts (with equal marginal
probabilities for treatment and non-treatment), and constructs a design
matrix along with the corresponding outcome. When `gen_ints = TRUE` the
full design matrix is returned (including interactions between
covariates and fixed effects and treatment indicators). When
`gen_ints = FALSE` the design matrix is generated in a simpler format
(with no interactions) as expected by
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).
Moreover, the covariates are generated according to the specified
`distribution`: by default, covariates are drawn from a normal
distribution; if `distribution = "uniform"`, they are drawn uniformly
from \\\[-\sqrt{3}, \sqrt{3}\]\\.

When \\d = 0\\ (i.e. no covariates), no covariate-related columns or
interactions are generated.

See the simulation studies section of Faletto (2025) for details.

## Usage

``` r
simulateDataCore(
  N,
  T,
  G = NULL,
  d,
  sig_eps_sq,
  sig_eps_c_sq,
  beta,
  seed = NULL,
  gen_ints = FALSE,
  distribution = "gaussian",
  guarantee_rank_condition = FALSE,
  assignment_type = "marginal",
  assignment_coefs = NULL,
  R = NULL
)
```

## Arguments

- N:

  Integer. Number of units in the panel.

- T:

  Integer. Number of time periods.

- G:

  Integer. Number of treated cohorts (with treatment starting in periods
  2 to T).

- d:

  Integer. Number of time-invariant covariates.

- sig_eps_sq:

  Numeric. Variance of the idiosyncratic (observation-level) noise.

- sig_eps_c_sq:

  Numeric. Variance of the unit-level random effects. Must be
  non-negative; `0` is allowed (yields a panel with no unit-level random
  effects).

- beta:

  Numeric vector. Coefficient vector for data generation. Its required
  length depends on the value of `gen_ints`:

  - If `gen_ints = TRUE` and `d > 0`, the expected length is \\p = G +
    (T-1) + d + dG + d(T-1) + num\\treats + num\\treats \times d\\,
    where \\num\\treats = T \times G - \frac{G(G+1)}{2}\\.

  - If `gen_ints = TRUE` and `d = 0`, the expected length is \\p = G +
    (T-1) + num\\treats\\.

  - If `gen_ints = FALSE`, the expected length is \\p = G + (T-1) + d +
    num\\treats\\.

- seed:

  (Optional) Integer. Seed for reproducibility.

- gen_ints:

  Logical. If `TRUE`, generate the full design matrix with interactions;
  if `FALSE` (the default), generate a design matrix without any
  interaction terms.

- distribution:

  Character. Distribution to generate covariates. Defaults to
  `"gaussian"`. If set to `"uniform"`, covariates are drawn uniformly
  from \\\[-\sqrt{3}, \sqrt{3}\]\\. To obtain a matching ground truth
  from
  [`getTes`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md),
  pass the same `distribution` value there.

- guarantee_rank_condition:

  (Optional). Logical. If TRUE, the returned data set is guaranteed to
  have at least `d + 1` units per cohort, which is necessary for the
  final design matrix to have full column rank. Default is FALSE, in
  which case only `>= 1` unit per cohort is required – this permits
  small cohorts and therefore high-dimensional (`p > NT`) panels, which
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  fits in its regularized regime (and on which the high-dimensional
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  path can be exercised, using the known data-generating coefficient
  vector; recovering the truth requires an adequate number of units).

- assignment_type:

  Character. One of `"marginal"` (default), `"multinomial"`, or
  `"ordered"`. Selects the cohort-assignment DGP. `"marginal"` preserves
  the pre-1.14.0 behavior. The non-marginal types require a non-NULL
  `assignment_coefs` argument (typically pulled from a `FETWFE_coefs`
  object built with the matching `assignment_type`).

- assignment_coefs:

  Optional list returned by `.gen_assignment_coefs()` (an internal
  helper). Required when `assignment_type != "marginal"`.

- R:

  Deprecated. The former name for `G`; still accepted with a deprecation
  warning, and will be removed in a future release. Use `G`.

## Value

An object of class `"FETWFE_simulated"`, which is a list containing:

- pdata:

  A dataframe containing generated data that can be passed to
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).

- X:

  The design matrix. When `gen_ints = TRUE`, \\X\\ has \\p\\ columns
  with interactions; when `gen_ints = FALSE`, \\X\\ has no interactions.

- y:

  A numeric vector of length \\N \times T\\ containing the generated
  responses.

- covs:

  A character vector containing the names of the generated features (if
  \\d \> 0\\), or simply an empty vector (if \\d = 0\\)

- time_var:

  The name of the time variable in pdata

- unit_var:

  The name of the unit variable in pdata

- treatment:

  The name of the treatment variable in pdata

- response:

  The name of the response variable in pdata

- coefs:

  The coefficient vector \\\beta\\ used for data generation.

- first_inds:

  A vector of indices indicating the first treatment effect for each
  treated cohort.

- N_UNTREATED:

  The number of never-treated units.

- assignments:

  A vector of counts (of length \\G+1\\) indicating how many units fall
  into the never-treated group and each of the \\G\\ treated cohorts.

- indep_counts:

  Independent cohort assignments (for auxiliary purposes).

- p:

  The number of columns in the design matrix \\X\\.

- N:

  Number of units.

- T:

  Number of time periods.

- G:

  Number of treated cohorts.

- R:

  Deprecated alias for `G`, retained for backward compatibility;
  populated with the same value. Use `G`. Will be removed in a future
  release.

- d:

  Number of covariates.

- sig_eps_sq:

  The idiosyncratic noise variance.

- sig_eps_c_sq:

  The unit-level noise variance.

## Details

When `gen_ints = TRUE`, the function constructs the design matrix by
first generating base fixed effects and a long-format covariate matrix
(via `generateBaseEffects()`), then appending interactions between the
covariates and cohort/time fixed effects (via `generateFEInts()`) and
finally treatment indicator columns and treatment-covariate interactions
(via `genTreatVarsSim()` and `genTreatInts()`). When `gen_ints = FALSE`,
the design matrix consists only of the base fixed effects, covariates,
and treatment indicators.

The argument `distribution` controls the generation of covariates. For
`"gaussian"`, covariates are drawn from `rnorm`; for `"uniform"`, they
are drawn from `runif` on the interval \\\[-\sqrt{3}, \sqrt{3}\]\\.

When \\d = 0\\ (i.e. no covariates), the function omits any
covariate-related columns and their interactions.

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

## Examples
