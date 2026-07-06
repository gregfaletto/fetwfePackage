# Tidy a `FETWFE_tes` simulation truth object

Returns a `broom`-style tidy data frame for the population-truth object
returned by
[`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md).
Row 1 is the overall true ATT (`term = "ATT_true"`); subsequent rows are
the true cohort ATTs (`term = "Cohort <adoption-time>"`, using the
simulator's convention that cohort `g` adopts at calendar time `g + 1`,
so the labels match what `tidy.<estimator>` uses on a fitted panel
generated from the same `FETWFE_coefs`). Standard error / statistic /
p-value columns are always `NA_real_` — there is no sampling
distribution for a population truth. When `conf.int = TRUE` (default,
matching the sibling tidy methods), `conf.low` / `conf.high` columns are
included and also set to `NA_real_`. When `conf.int = FALSE`, those
columns are omitted.

## Usage

``` r
# S3 method for class 'FETWFE_tes'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  An object of class `"FETWFE_tes"` returned by
  [`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md).

- conf.int:

  Logical; include `conf.low` / `conf.high` columns. Defaults to `TRUE`
  to match the sibling tidy methods and preserve pre-#84 backward
  compatibility. Population-truth objects have no sampling distribution,
  so the CI columns are always filled with `NA_real_` when included.

- conf.level:

  Numeric in (0, 1). Accepted for broom-convention parity but unused (no
  CIs to compute for a population truth); validated regardless. Defaults
  to `0.95` (`FETWFE_tes` objects do not carry an alpha slot, so there
  is no fitted-object value to default to).

- ...:

  Unused.

## Value

A data frame with `G + 1` rows and columns `term`, `estimate`,
`std.error`, `statistic`, `p.value`, and (when `conf.int = TRUE`)
`conf.low` / `conf.high`.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  broom::tidy(getTes(coefs))
} # }
```
