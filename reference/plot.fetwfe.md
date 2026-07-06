# Plot CATT or event-study estimates from a fitted FETWFE / ETWFE / BETWFE

Returns a ggplot object showing either event-study coefficients
(default; `type = "event_study"`) or per-cohort average treatment
effects (`type = "catt"`) from a fitted estimator. Mirrors the
visualization style of
[`did::ggdid()`](https://bcallaway11.github.io/did/reference/ggdid.html)
from the Callaway-Sant'Anna `did` package, providing a single-call route
from a fitted object to a publication-ready visualization.

For `fetwfe` / `betwfe` (the bridge-penalty estimators) in the CATT
view, points are shape- and color-coded by whether the bridge penalty
left that cohort's ATT nonzero (`selected = TRUE`) or zeroed it out
(`selected = FALSE`). For `etwfe` (no selection), all points are
uniformly styled.

## Usage

``` r
# S3 method for class 'fetwfe'
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

A `ggplot` object. Users can customize further via standard ggplot
layer-addition syntax (e.g., `plot(res) + ggplot2::theme_classic()`).

## See also

[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
for the per-cohort accessor;
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
for the event-time accessor;
[`plot.etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/plot.etwfe.md),
[`plot.betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/plot.betwfe.md)
for the parallel methods.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- fetwfeWithSimulatedData(dat)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot(res)                          # default: event-study coefficients
    plot(res, type = "catt")           # per-cohort ATTs
    plot(res, conf_int = FALSE)        # point estimates only
    plot(res, alpha = 0.1)             # 90% CIs
  }
} # }
```
