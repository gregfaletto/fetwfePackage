# Glance an `fetwfe` fitted object

Returns a one-row `broom`-style summary data frame with model-level
scalars: panel-shape counts (`nobs`, `n_units`, `n_periods`,
`n_cohorts`, `n_covs`, `n_features`), bridge-regression tuning
(`lambda_star`, `lambda_star_model_size`), variance components
(`sig_eps_sq`, `sig_eps_c_sq`), and inference settings (`alpha`,
`se_type`, `indep_counts_used`).

## Usage

``` r
# S3 method for class 'fetwfe'
glance(x, ...)
```

## Arguments

- x:

  An object of class `"fetwfe"`.

- ...:

  Unused.

## Value

A one-row data frame with 16 columns.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- fetwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::glance(res)
} # }
```
