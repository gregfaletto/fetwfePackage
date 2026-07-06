# Tidy an `eventStudy` object

Returns a `broom`-style tidy data frame for the output of
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md).
Renames existing columns to broom conventions (`se` \\\to\\ `std.error`,
`p_value` \\\to\\ `p.value`) and adds a `term` column
(`"e<event_time>"`) plus a `statistic` column (`estimate / std.error`)
so the schema matches `tidy.<estimator>()` for downstream `bind_rows()`
consumers.

## Usage

``` r
# S3 method for class 'eventStudy'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  An object of class `"eventStudy"` returned by
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md).

- conf.int:

  Logical; include `conf.low` / `conf.high` columns.

- conf.level:

  Numeric in (0, 1). Retained for `broom`-convention parity (default
  `0.95`) but no longer recomputes the event-study CIs: as of version
  1.16.0 (#197) the `conf.low` / `conf.high` columns pass through the
  `eventStudy` object's stored `ci_low` / `ci_high` (reflecting the
  fit's `ci_type`). To change the confidence level, recompute
  `eventStudy(fit, alpha = ...)` at the desired alpha first.

- ...:

  Unused.

## Value

A data frame with one row per event-time and columns `term`,
`event_time`, `n_cohorts`, `estimate`, `std.error`, `statistic`,
`p.value`, and (when `conf.int = TRUE`) `conf.low` / `conf.high`.

## Details

The
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
output stores its confidence-interval bounds (`ci_low` / `ci_high`),
which reflect the fit's `ci_type` (#197): simultaneous (family-wise,
uniform) by default, or pointwise when the fit used
`ci_type = "pointwise"`. When `conf.int = TRUE` (the default),
`conf.low` / `conf.high` PASS THROUGH those stored bounds rather than
recomputing from `estimate +/- z * se` — so the tidied event-study CIs
agree with `print` / `summary` / `plot` and with
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
under the default. When `conf.int = FALSE`, the CI columns are omitted.
(Degenerate event times carry `NA` bounds under both `ci_type`
settings.)

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- fetwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(eventStudy(res))
} # }
```
