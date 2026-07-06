# Plot CATT or event-study estimates from a fitted BETWFE

Parallel to
[`plot.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/plot.fetwfe.md).
BETWFE uses the same bridge-penalty selection mechanism, so `selected`
is encoded the same way (TRUE = nonzero, FALSE = zeroed by the bridge
penalty).

## Usage

``` r
# S3 method for class 'betwfe'
plot(x, type = c("event_study", "catt"), conf_int = TRUE, alpha = NULL, ...)
```

## Arguments

- x:

  A fitted object from
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
  or
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md).
  (`twfeCovs` is not currently supported – see GitHub issue \#58 for the
  broader treatment of `twfeCovs` class methods.)

- type:

  Character; either `"event_study"` (default; event-time pooled
  coefficients from
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md))
  or `"catt"` (per-cohort ATTs from `result$catt_df`).

- conf_int:

  Logical; if `TRUE` (default), include confidence interval error bars.

- alpha:

  Numeric; overrides the fit's alpha for CI computation. `NULL`
  (default) uses the fit's alpha (so 95% CIs when the fit's alpha is
  0.05). With the default `NULL`, both views show the fit's `ci_type`
  band (simultaneous by default): the `"catt"` view reads the bounds
  stored in `catt_df`, and the `"event_study"` view re-derives them via
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
  (which inherits the fit's `ci_type`). **Asymmetry under an explicit
  `alpha`:** for the `"event_study"` view, the explicit alpha is
  forwarded to
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md),
  which re-runs the joint machinery and so still produces a
  *simultaneous* band at that alpha (when the fit's `ci_type` is
  `"simultaneous"`); for the `"catt"` view, an explicit alpha recomputes
  `ci_low` / `ci_high` from `estimate +/- qnorm(1 - alpha/2) * se`, i.e.
  a *pointwise* band at that alpha (the stored simultaneous bounds are
  at the fit's alpha and are not re-derived at a new alpha for the catt
  view). To plot simultaneous catt bands at a different alpha, refit at
  that alpha or call
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
  directly.

- ...:

  Currently unused; reserved for future arguments.

## Value

A `ggplot` object.

## See also

[`plot.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/plot.fetwfe.md)
for the full documentation;
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md);
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- betwfeWithSimulatedData(dat)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot(res)
  }
} # }
```
