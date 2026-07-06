# Glance a `twfeCovs` fitted object

Like
[`glance.etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/glance.etwfe.md)
(and with the same schema): omits the `lambda_star` /
`lambda_star_model_size` columns, since `twfeCovs` performs no
regularization.

## Usage

``` r
# S3 method for class 'twfeCovs'
glance(x, ...)
```

## Arguments

- x:

  An object of class `"twfeCovs"`.

- ...:

  Unused.

## Value

A one-row data frame with 11 columns.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- twfeCovsWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::glance(res)
} # }
```
