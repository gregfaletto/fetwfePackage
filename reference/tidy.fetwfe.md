# Tidy an `fetwfe` fitted object

Returns a `broom`-style tidy data frame for an object of class
`"fetwfe"`. Row 1 is the overall ATT (`term = "ATT"`); subsequent rows
are the cohort-specific ATTs (`term = "Cohort <adoption-time>"`), one
per treated cohort, sorted by ascending cohort label. Standard error,
z-statistic, and p-value reflect the value of `se_type` used at fit time
(model-based by default, cluster-robust under `se_type = "cluster"`).
Cohorts that the bridge penalty zeroed out (`selected = FALSE`) carry
`NA` for `std.error` / `statistic` / `p.value`.

## Usage

``` r
# S3 method for class 'fetwfe'
tidy(x, conf.int = TRUE, conf.level = 1 - x$alpha, ...)
```

## Arguments

- x:

  An object of class `"fetwfe"` returned by
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).

- conf.int:

  Logical. If `TRUE` (default), `conf.low` and `conf.high` columns are
  included.

- conf.level:

  Numeric in (0, 1). Applies only to the overall-ATT row (row 1), whose
  CI is recomputed at this level; defaults to `1 - x$alpha` (faithful to
  the alpha used at fit time; deviates from
  [`broom::tidy.lm`](https://broom.tidymodels.org/reference/tidy.lm.html)'s
  `0.95` default by design). The cohort rows pass through the fit-time
  `catt_df` bounds (reflecting the fit's `ci_type` and the fit-time
  alpha) and are NOT recomputed at `conf.level` (#197) — mirroring
  [`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md)'s
  stored-bounds behavior.

- ...:

  Unused; present for S3 compatibility.

## Value

A data frame with `G + 1` rows and columns `term`, `estimate`,
`std.error`, `statistic`, `p.value`, optionally `conf.low` /
`conf.high`, and `selected` (logical).

## Details

The cohort-row `conf.low` / `conf.high` columns pass through the fit's
stored `catt_df` bounds, so they reflect the fit's `ci_type` (#197):
simultaneous (family-wise) by default, or pointwise when the fit used
`ci_type = "pointwise"`. They are NOT recomputed from `conf.level` (see
the `conf.level` note). The overall-ATT row (row 1) is a scalar, so its
CI is the pointwise Wald interval at `conf.level` (pointwise ==
simultaneous for a single effect).

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- fetwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(res)
} # }
```
