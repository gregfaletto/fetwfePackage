# Augment user-supplied data with fitted values and residuals from an etwfe fit

Same shape as
[`augment.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/augment.fetwfe.md),
dispatched on class `"etwfe"`. `data` is auto-sorted by `(unit, time)`
and any first-period-treated units are auto-trimmed; pass the same raw
`pdata` you handed to
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md).

## Usage

``` r
# S3 method for class 'etwfe'
augment(x, data, ...)
```

## Arguments

- x:

  An object of class `"etwfe"`.

- data:

  A panel `data.frame` with the response column under
  `x$response_col_name`. Any sort order; first-period-treated units are
  auto-trimmed.

- ...:

  Unused.

## Value

`data` with `.fitted` and `.resid` columns appended.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5,
                               eff_size = 2),
                      N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- etwfeWithSimulatedData(sim)
  broom::augment(res, data = sim$pdata)
} # }
```
