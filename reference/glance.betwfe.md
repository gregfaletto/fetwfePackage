# Glance a `betwfe` fitted object

Same schema as
[`glance.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/glance.fetwfe.md)
(BETWFE also has regularization).

## Usage

``` r
# S3 method for class 'betwfe'
glance(x, ...)
```

## Arguments

- x:

  An object of class `"betwfe"`.

- ...:

  Unused.

## Value

A one-row data frame with 16 columns.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- betwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::glance(res)
} # }
```
