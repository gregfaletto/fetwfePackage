# Tidy a `betwfe` fitted object

Like
[`tidy.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.fetwfe.md)
but for a BETWFE fit. Includes the `selected` column reflecting BETWFE's
bridge-penalized selection.

## Usage

``` r
# S3 method for class 'betwfe'
tidy(x, conf.int = TRUE, conf.level = 1 - x$alpha, ...)
```

## Arguments

- x:

  An object of class `"betwfe"` returned by
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md).

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
  res <- betwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(res)
} # }
```
