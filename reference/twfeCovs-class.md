# TWFE-With-Covariates Output Class

S3 class for the output of
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md).
Carries the same styled `print` / `summary` / `coef` surface as the
three sibling estimators, plus `tidy` / `glance` (broom) and
`simultaneousCIs`. `plot` and `augment` are intentionally not provided:
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
estimates one pooled effect per cohort, so it has no per-(cohort, time)
/ event-study basis to plot, and its coefficient vector is in a reduced
basis that `augment()`'s fitted-value path does not match (#58). Both
raise an informative error.
