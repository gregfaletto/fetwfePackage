# cran-comments.md

**Package:** fetwfe
**Version:** 1.10.0
**Date:** 2026-05-21

This is a feature-and-bugfix update. The previous CRAN release is 1.5.0; the
1.5.1-1.9.33 versions listed in `NEWS.md` are unreleased development.

## What's new since 1.5.0

### New functionality
- `event_study()`, with `plot()` methods for the estimator classes, for
  pooled event-time treatment-effect estimates and plots.
- broom `tidy()`, `glance()`, and `augment()` methods for the estimator
  outputs (and for `event_study()` and `getTes()`).
- An experimental cluster-robust standard-error option (`se_type =
  "cluster"`) on all four estimators; default behavior is unchanged.
- An `allow_no_never_treated` argument that auto-truncates panels with no
  never-treated units instead of erroring.
- Per-cohort `P_value` columns, and `selected` columns for the bridge
  estimators, with a vignette section on interpreting them.
- S3 classes for the `betwfe()` and `getTes()` outputs; two new vignettes.

### Bug fixes affecting numerical output
- Most significant for existing users: the internal variance-component
  estimator was corrected. It previously collapsed the unit-level variance
  component to nearly zero, degenerating the GLS weighting step to the
  identity; it is now estimated by REML. Both standard errors and point
  estimates shift after this update on panels with non-negligible
  unit-level variance. This is a bug fix, not an API change.
- An off-by-one error in a variance-component Jacobian was corrected,
  shifting `fetwfe()`'s reported overall standard error.
- Smaller correctness fixes (e.g. `betwfe()` standard errors under partial
  selection; cohort ordering in `catt_df`; column preservation in
  `augment()`).

### Internal
- Extensive internal hardening since 1.5.0 — object and method-entry
  validators, a source-tree reorganization, refactors removing duplicated
  logic, and expanded test coverage — with no user-facing behavior change.

## Test environments
- Local: macOS 15.5 (Apple Silicon), R 4.5.0.
- win-builder: R-devel, via `devtools::check_win_devel()`.
- macOS builder: R-release, via `devtools::check_mac_release()`.

## R CMD check results
0 errors | 0 warnings | 0 notes.

A local `R CMD check --as-cran` additionally reports the environmental NOTE
"unable to verify current time"; this arises only because the local check
machine cannot reach a time server, and does not occur on CRAN's build
infrastructure.

## Reverse dependencies
`fetwfe` has no reverse dependencies on CRAN.

## Dependency changes since 1.5.0
- `generics` added to `Imports` (used by the new broom methods).
- `broom`, `ggplot2`, and `lme4` added to `Suggests`.
- `expm` moved from `Imports` to `Suggests`.

## Other notes
- No compiled code; the package is pure R. No new system requirements.
- Requires R (>= 4.1.0).

–––
Thank you for your time reviewing this submission.
