# Augment user-supplied data with fitted values and residuals from a betwfe fit

Same shape as
[`augment.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/augment.fetwfe.md),
dispatched on class `"betwfe"`. `data` is auto-sorted by `(unit, time)`
and any first-period-treated units are auto-trimmed; pass the same raw
`pdata` you handed to
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md).

## Usage

``` r
# S3 method for class 'betwfe'
augment(x, data, ...)
```

## Arguments

- x:

  An object of class `"betwfe"`.

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
  res <- betwfeWithSimulatedData(sim)
  broom::augment(res, data = sim$pdata)
} # }
```
