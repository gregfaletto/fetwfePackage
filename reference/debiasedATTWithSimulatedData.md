# Run debiasedATT() on simulated data

Convenience wrapper mirroring the `*WithSimulatedData()` pattern: fits
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
on a `"FETWFE_simulated"` object (via
[`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md))
and returns
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
on the fit.

## Usage

``` r
debiasedATTWithSimulatedData(simulated_obj, q = 0.5, alpha = NULL, ...)
```

## Arguments

- simulated_obj:

  An object of class `"FETWFE_simulated"` (from
  [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)).

- q:

  Numeric; the `L_q` bridge penalty exponent. Must be `< 1` (the
  debiased SE requires bridge selection). Defaults to `0.5`.

- alpha:

  Numeric in `(0, 1)` or `NULL`; confidence level `1 - alpha`. Default
  `NULL`, inheriting the `alpha` stored on the fit.

- ...:

  Further arguments forwarded to
  [`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md)
  (e.g. `fusion_structure`, `lambda_selection`, `cv_seed`).

## Value

An object of class `"debiased_att"` (see
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)).

## See also

[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md),
[`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 1)
  dat <- simulateData(coefs, N = 200, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
  debiasedATTWithSimulatedData(dat)
} # }
```
