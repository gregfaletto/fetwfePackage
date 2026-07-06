# Tidy an `etwfe` fitted object

Like
[`tidy.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.fetwfe.md)
but for an ETWFE fit. Has no `selected` column (ETWFE does no
regularized selection).

## Usage

``` r
# S3 method for class 'etwfe'
tidy(x, conf.int = TRUE, conf.level = 1 - x$alpha, ...)
```

## Arguments

- x:

  An object of class `"etwfe"` returned by
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md).

- conf.int:

  Logical; include CI columns.

- conf.level:

  Numeric in (0, 1); defaults to `1 - x$alpha`. Applies only to the
  overall-ATT row; the cohort rows pass through the fit-time `catt_df`
  bounds (reflecting the fit's `ci_type`) and are not recomputed at
  `conf.level` (#197). See
  [`tidy.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.fetwfe.md).

- ...:

  Unused.

## Value

A data frame with `G + 1` rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- etwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(res)
} # }
```
