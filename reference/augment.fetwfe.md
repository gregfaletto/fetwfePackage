# Augment user-supplied data with fitted values and residuals from a fetwfe fit

Computes `.fitted = X %*% beta_hat + x$y_mean` and
`.resid = data[[x$response_col_name]] - .fitted`, then column-binds
those two columns onto `data`. The response mean and column name are
stored on the fitted object during fitting (the estimator internally
centers `y` before solving), so fitted values come back on the
original-response scale without the caller having to remember either.

## Usage

``` r
# S3 method for class 'fetwfe'
augment(x, data, ...)
```

## Arguments

- x:

  An object of class `"fetwfe"`.

- data:

  A panel `data.frame` with one row per unit-period (any sort order —
  augment auto-sorts), containing the response column under the same
  name used at fit time (see `x$response_col_name`).
  First-period-treated units, if present, are auto-trimmed.

- ...:

  Unused.

## Value

A copy of `data` with two extra numeric columns: `.fitted` and `.resid`.

## Details

`data` is auto-handled to match the fitted design: rows are auto-sorted
by `(unit, time)`, and any first-period-treated units (whose treatment
effect cannot be identified by the estimator) are auto-trimmed via
`idCohorts()`. So you can pass the same raw `pdata` you handed to
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
— the method takes care of alignment. The only hard requirement is that
`data` contains the response column under its original name.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5,
                               eff_size = 2),
                      N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- fetwfeWithSimulatedData(sim)
  broom::augment(res, data = sim$pdata)
} # }
```
