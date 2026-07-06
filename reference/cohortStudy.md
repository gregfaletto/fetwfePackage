# Per-cohort average treatment effects

Extracts per-cohort ATT estimates from a fitted FETWFE / ETWFE / BETWFE
/ twfeCovs object as a tidy data frame. Parallel to
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
for event-time aggregation: the per-cohort information is already
available via `result$catt_df`, and `cohortStudy()` surfaces it through
a discoverable function with its own help page (so users can reach it
via `?cohortStudy` without having to know the slot name).

The function is a pass-through on `result$catt_df` modulo class: the
columns, their values, and their order are unchanged. The returned
object carries class `c("cohortStudy", "catt_df", "data.frame")`. The
`cohortStudy` class dispatches the broom tidier
[`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md);
the `catt_df` class preserves the helpful-error layer (introduced in
fetwfe 1.11.0) that intercepts pre-1.11.0 Title-Case column names
(`Cohort`, `Estimated TE`, etc.) with a migration message; the
`data.frame` base preserves the standard data-frame methods (`print`,
`head`, `nrow`,
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html),
etc.).

## Usage

``` r
cohortStudy(result)
```

## Arguments

- result:

  A fitted object from
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
  or
  [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
  (or their `*WithSimulatedData()` wrapper analogs, which return the
  same classes).

## Value

A data frame with class `c("cohortStudy", "catt_df", "data.frame")`
containing one row per treated cohort and columns:

- cohort:

  Character; the cohort label (the calendar time at which the cohort
  first received treatment).

- estimate:

  Numeric; the per-cohort ATT estimate.

- se:

  Numeric; standard error for the per-cohort ATT (`NA` when the Gram
  matrix is singular or, for
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  /
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
  the bridge penalty zeroed out the cohort).

- ci_low, ci_high:

  Numeric; the stored lower and upper confidence-interval bounds,
  reflecting the fit's `ci_type` – simultaneous (family-wise) by
  default, or pointwise `1 - alpha` Wald bounds when the fit used
  `ci_type = "pointwise"` (`alpha` is the value passed at fit time).

- p_value:

  Numeric; follows the fit's `ci_type`. Under `"pointwise"`, the
  two-sided Wald p-value (`2 * pnorm(-|estimate / se|)`); under
  `"simultaneous"` (the default), the single-step max-T
  multiplicity-adjusted (family-wise) p-value matching the simultaneous
  band (#200). `NA` when `se` is `0` or `NA`.

- selected:

  ([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  /
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
  only.) Logical; `TRUE` when the bridge penalty left the cohort's ATT
  nonzero. Absent for
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  and
  [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md),
  which do not perform selection.

Use `tidy(cohortStudy(result))` (with the `broom` package loaded) to
reshape to broom convention (`term`, `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low`, `conf.high`, optionally `selected`);
see
[`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md).

## See also

[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
for the parallel event-time accessor;
[`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md)
for the fully disaggregated per-(cohort, time) accessor;
[`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md)
for broom-shape translation.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- fetwfeWithSimulatedData(dat)
  cs <- cohortStudy(res)
  cs
  # Broom-shape translation:
  if (requireNamespace("broom", quietly = TRUE)) {
    broom::tidy(cs)
  }
} # }
```
