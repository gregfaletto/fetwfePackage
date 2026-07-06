# Compute pooled event-time treatment-effect estimates

For a fitted object from
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
or
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
computes the pooled event-time treatment-effect estimates `tau_E(e)`,
defined as cohort-weighted averages of the cell-level treatment-effect
estimates at each post-treatment event time `e = t - g` (where `t` is
calendar time and `g` is the cohort's first-treated calendar time).
Weights are sample-cohort-size weights (matching
`did::aggte(type = "dynamic")` convention).

Standard errors combine two terms, mirroring the package's existing
overall-ATT SE machinery: `var_1(e)` from regression-coefficient noise
(computed via the same `gram_inv` machinery the package uses for cohort
SEs, or the cluster-robust sandwich under `se_type = "cluster"`), and
`var_2(e)` from cohort-probability noise (analog of the existing
`getSecondVarTermOLS` / `getSecondVarTermDataApp` machinery, with the
multinomial Jacobian restricted to cohorts valid at event time `e`).
Combined as `sqrt(var_1 + var_2)` by default (asymptotically exact under
paper Theorem `te.asym.norm.thm`(\\c'\\) / Assumption (\\\Psi\\-IF),
which the package's default cohort-sample-proportions estimator
satisfies); the conservative Cauchy-Schwarz bound
`sqrt(var_1 + var_2 + 2 sqrt(var_1 * var_2))` is available via
`se_type = "conservative"` (for users with non-(\\\Psi\\-IF)
propensity-score estimators). When `indep_counts` was supplied at fit
time, the tight formula applies regardless of `se_type` (two-sample
regime, Theorem (b)).

## Usage

``` r
eventStudy(x, alpha = NULL, ci_type = NULL)
```

## Arguments

- x:

  A fitted object of class `"fetwfe"`, `"etwfe"`, or `"betwfe"`.

- alpha:

  (Optional) Significance level for confidence intervals. Defaults to
  `x$alpha` (the alpha used at fit time).

- ci_type:

  (Optional) Character; one of `"simultaneous"` or `"pointwise"`, or
  `NULL` (default). Controls whether the returned `ci_low` / `ci_high`
  are the event-study-family simultaneous (family-wise, uniform) band or
  the per-event-time pointwise Wald intervals. `NULL` inherits the fit's
  stored `ci_type` (so a fit made with the default
  `ci_type = "simultaneous"` yields simultaneous event-study bands,
  which `print` / `summary` / `plot` then display); for objects fitted
  before version 1.16.0 (no `ci_type` slot) `NULL` falls back to
  `"pointwise"`. The `se` column is identical under both settings; the
  interval bounds, their critical-value multiplier, and the `p_value`
  differ (under `"simultaneous"` the `p_value` is the single-step max-T
  multiplicity- adjusted dual of the band, under `"pointwise"` the
  per-effect Wald p-value; \#200). Degenerate event times (no valid
  contributing cohorts) carry `NA` bounds and `NA` `p_value` under both
  settings.

## Value

A data frame with class `c("eventStudy", "data.frame")` and columns:

- event_time:

  Integer; event time `e = t - g`, ranging from 0 to `T - 2`.

- n_cohorts:

  Integer; number of cohorts contributing to the pooled estimate at
  event time `e`.

- estimate:

  Numeric; the pooled event-time ATT estimate.

- se:

  Numeric; combined standard error.

- ci_low:

  Numeric; lower bound of the (1 - alpha) Wald CI.

- ci_high:

  Numeric; upper bound of the (1 - alpha) Wald CI.

- p_value:

  Numeric; follows the fit's `ci_type`. Under `"pointwise"`, the
  two-sided Wald p-value (`2 * pnorm(-|estimate / se|)`); under
  `"simultaneous"` (the default), the single-step max-T
  multiplicity-adjusted (family-wise) p-value matching the simultaneous
  band (#200). `NA` when `se` is `0` or `NA`.

Only post-treatment event times (`e >= 0`) are included; pre-treatment
placebo periods would require an extended regression specification and
are out of scope for this initial release.

## See also

[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
for the parallel per-cohort accessor;
[`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md)
for the fully disaggregated per-(cohort, time) accessor.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- fetwfeWithSimulatedData(dat)
  eventStudy(res)
} # }
```
