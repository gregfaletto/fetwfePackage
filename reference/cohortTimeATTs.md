# Per-(cohort, time) average treatment effects

Extracts the fully disaggregated treatment-effect estimates from a
fitted FETWFE / ETWFE / BETWFE object: one row for every
`(cohort, time)` cell, with no averaging over cohorts or over event
time. This is the finest-grained view of the estimated effects — the
`num_treats` underlying parameters themselves — complementing
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
(which averages each cohort's cells over time) and
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
(which averages over cohorts at each event time).

Like
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md),
this accessor is not available for
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
objects: that estimator has a single treatment-effect parameter per
cohort (no per-time disaggregation), so its finest granularity is
already
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md).

Standard errors are the per-cell regression standard errors
\\\sqrt{\sigma\_\varepsilon^2 \\ \psi' G^{-1} \psi / (NT)}\\, recomputed
at call time from the fit's stored design (the same Gram-matrix
machinery
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
uses; nothing is added to the fitted object). Because each cell is a
single cohort-time parameter, the cohort-probability sampling variance
that contributes to the aggregated
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
/
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
standard errors is identically zero here, so a cell's SE is a single
coefficient's regression SE.

Confidence intervals and p-values are *pointwise* \\1 - \alpha\\ Wald
quantities (`estimate +/- z * se`). For simultaneous (family-wise) bands
over the cell family, use
`simultaneousCIs(result, family = "all_post_treatment")`.

## Usage

``` r
cohortTimeATTs(result, alpha = NULL)
```

## Arguments

- result:

  A fitted object from
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
  or
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
  (or their `*WithSimulatedData()` wrapper analogs, which return the
  same classes).
  [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
  objects are not supported (see Description).

- alpha:

  Numeric in `(0, 1)`; the pointwise confidence level is `1 - alpha`.
  Defaults to the `alpha` stored on the fit.

## Value

A data frame with class `c("cohortTimeATTs", "data.frame")` containing
one row per `(cohort, time)` treatment-effect cell, sorted by cohort
then time, with columns:

- cohort:

  Character; the cohort label (the calendar time at which the cohort
  first received treatment), matching
  [`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md).

- time:

  Numeric; the calendar time of the cell, equal to the cohort's adoption
  time plus the event time (`0, 1, ...`). Real panels carry their actual
  calendar times. For synthetic
  [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
  /
  [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
  fixtures (whose panel runs `1, ..., T`, so the stored first year is
  `1`) this coincides with the 1-based panel-time index. (Only a
  hand-built or legacy fit with no stored first year falls back to that
  panel-time index directly.)

- estimate:

  Numeric; the cell's ATT estimate.

- se:

  Numeric; the pointwise standard error. `0` for a cell zeroed out by
  the fusion/bridge penalty while other cells survive
  ([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  /
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md));
  `NA` when standard errors are unavailable — the fit was computed with
  `q >= 1`, the Gram matrix on the selected support is singular, or the
  penalty zeroed the *entire* treatment block (no cells selected, so
  there is no support to recompute the Gram from; this matches
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)).

- ci_low, ci_high:

  Numeric; the pointwise `1 - alpha` Wald bounds
  `estimate -/+ qnorm(1 - alpha/2) * se`. `(0, 0)` for a fused-away
  cell; `NA` when `se` is `NA`.

- p_value:

  Numeric; the two-sided pointwise Wald p-value
  `2 * pnorm(-|estimate / se|)`. `NA` when `se` is `0` or `NA`.

- selected:

  ([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  /
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
  only.) Logical; `TRUE` when the bridge penalty left the cell's
  estimate nonzero. Absent for
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
  which does not perform selection.

Use `tidy(cohortTimeATTs(result))` (with the `broom` package loaded) to
reshape to broom convention; see
[`tidy.cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortTimeATTs.md).

## Details

The cell standard error is computed from `psi`, the cell's row of the
(selected) treatment-effect design — for
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
the relevant row of the inverse fusion transform \\D^{-1}\\ in the
transformed (theta) coordinate space, for
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
/
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
a unit selector in the original (beta) coordinate space restricted to
the selected support. It is never gated on the point estimate: a cell
whose estimate is exactly zero because the penalty fused it away has an
all-zero `psi` and therefore `se = 0` (the correct degenerate value),
while a cell whose estimate happens to be near zero for other reasons
still receives its proper nonzero SE.

## See also

[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
for the per-cohort (time-averaged) accessor;
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
for the per-event-time (cohort-averaged) accessor;
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
for simultaneous (family-wise) bands over the cell family
(`family = "all_post_treatment"`);
[`tidy.cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortTimeATTs.md)
for broom-shape translation.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  res <- fetwfeWithSimulatedData(dat)
  cta <- cohortTimeATTs(res)
  cta
  # Broom-shape translation:
  if (requireNamespace("broom", quietly = TRUE)) {
    broom::tidy(cta)
  }
} # }
```
