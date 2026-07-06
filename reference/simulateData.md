# Generate Random Panel Data for FETWFE Simulations

Generates a random panel data set for simulation studies of the fused
extended two-way fixed effects (FETWFE) estimator by taking an object of
class `"FETWFE_coefs"` (produced by
[`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md))
and using it to simulate data. The function creates a balanced panel
with \\N\\ units over \\T\\ time periods, assigns treatment status
across \\G\\ treated cohorts (with equal marginal probabilities for
treatment and non-treatment), and constructs a design matrix along with
the corresponding outcome. The covariates are generated according to the
specified `distribution`: by default, covariates are drawn from a normal
distribution; if `distribution = "uniform"`, they are drawn uniformly
from \\\[-\sqrt{3}, \sqrt{3}\]\\. When \\d = 0\\ (i.e. no covariates),
no covariate-related columns or interactions are generated. See the
simulation studies section of Faletto (2025) for details.

## Usage

``` r
simulateData(
  coefs_obj,
  N,
  sig_eps_sq,
  sig_eps_c_sq,
  distribution = "gaussian",
  guarantee_rank_condition = FALSE,
  seed = NULL
)
```

## Arguments

- coefs_obj:

  An object of class `"FETWFE_coefs"` containing the coefficient vector
  and simulation parameters.

- N:

  Integer. Number of units in the panel.

- sig_eps_sq:

  Numeric. Variance of the idiosyncratic (observation-level) noise.

- sig_eps_c_sq:

  Numeric. Variance of the unit-level random effects. Must be
  non-negative; `0` is allowed (yields a panel with no unit-level random
  effects).

- distribution:

  Character. Distribution to generate covariates. Defaults to
  `"gaussian"`. If set to `"uniform"`, covariates are drawn uniformly
  from \\\[-\sqrt{3}, \sqrt{3}\]\\.

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

- seed:

  (Optional) Controls the random-number generator for the simulated
  panel. As of fetwfe 1.24.0 the default is `NULL`, which draws from the
  ambient random-number generator (respecting any preceding
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) and emits a
  warning; pass an integer for a reproducible panel, or `NA` to draw
  from the ambient generator silently. `simulateData()` no longer reuses
  `coefs_obj$seed` (the seed
  [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
  used to build the coefficients). Default `NULL`.

## Value

An object of class `"FETWFE_simulated"`, which is a list containing:

- pdata:

  A dataframe containing generated data that can be passed to
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).

- X:

  The design matrix \\X\\, with \\p\\ columns with interactions.

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

This function extracts simulation parameters from the `FETWFE_coefs`
object and passes them, along with additional simulation parameters, to
the internal function
[`simulateDataCore()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateDataCore.md).
It validates that all necessary components are returned and assigns the
S3 class `"FETWFE_simulated"` to the output.

The random draw is controlled by the `seed` argument, not by
`coefs_obj$seed`. By default (`seed = NULL`) `simulateData()` draws from
the ambient random-number generator (so a preceding
[`set.seed()`](https://rdrr.io/r/base/Random.html) is respected and
repeated calls return different panels) and emits a warning noting that
this default changed in fetwfe 1.24.0. Pass an integer `seed` for a
reproducible panel (the same integer always yields the same panel), or
`seed = NA` to use the ambient generator without the warning. To vary
the panel across Monte Carlo replications, pass a different `seed` each
replication.

Passing an explicit numeric `seed` calls `set.seed(seed)` internally and
**leaves the global RNG advanced** after the call returns. This is
deliberate — it makes a simulation both reproducible (the same `seed`
always yields the same panel) and varying (subsequent draws consume the
advanced stream). Use `seed = NA` (or the default `seed = NULL`) to draw
from / preserve the ambient stream without calling
[`set.seed()`](https://rdrr.io/r/base/Random.html).

The argument `distribution` controls the generation of covariates. For
`"gaussian"`, covariates are drawn from `rnorm`; for `"uniform"`, they
are drawn from `runif` on the interval \\\[-\sqrt{3}, \sqrt{3}\]\\
(which ensures that the covariates have unit variance regardless of
which distribution is chosen).

When \\d = 0\\ (i.e. no covariates), the function omits any
covariate-related columns and their interactions.

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate coefficients
  coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)

  # Simulate data using the coefficients
  sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
} # }
```
