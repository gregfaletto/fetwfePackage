# Glance an `etwfe` fitted object

Like
[`glance.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/glance.fetwfe.md)
but omits the `lambda_star` / `lambda_star_model_size` columns — ETWFE
has no regularization.

## Usage

``` r
# S3 method for class 'etwfe'
glance(x, ...)
```

## Arguments

- x:

  An object of class `"etwfe"`.

- ...:

  Unused.

## Value

A one-row data frame with 11 columns.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- etwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::glance(res)
} # }
```
