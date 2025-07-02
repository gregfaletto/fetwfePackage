# cran-comments.md

**Package:** fetwfe  
**Version:** 1.4.1  
**Date:** 2025-07-01

## What’s new since version 1.0.0 (2025-05-14)

### Core estimation enhancements
- Introduced an S3 class for `fetwfe()` to enable method dispatch and future extensions.

### New estimators and wrappers
- **Extended TWFE**  
  - `etwfe()` for extended two-way fixed effects estimation  
  - `etwfeWithSimulatedData()`, a convenience wrapper for running `etwfe()` on simulated data
- **Bridge-penalized TWFE**  
  - `btwfe()` implementing a bridge-penalized extended two-way fixed effects estimator  
  - `btwfeWithSimulatedData()`, a wrapper for `btwfe()`
- **Two-way FE with covariates**  
  - `twfeCovs()` for standard two-way FE models including covariates  
  - `twfeCovsWithSimulatedData()`, a wrapper for `twfeCovs()`

### Simulation and coefficient utilities
- `genCoefs()` / `genCoefsCore()` to generate coefficient vectors for simulation studies
- `simulateData()` / `simulateDataCore()` to create panel data from a coefficient object  
  - New `guarantee_rank_condition` argument in `simulateData()` ensures the design matrix has full column rank
- `fetwfeWithSimulatedData()` to streamline estimation on simulated datasets
- `getTes()` to compute “true” treatment effects from a coefficient object

### Data-format conversion helpers
- `attgtToFetwfeDf()` converts `did::att_gt()` outputs into a format suitable for `fetwfe()`
- `ewtfeToFetwfeDf()` converts `etwfe::etwfe()` outputs into a `fetwfe()`‐compatible data frame

### Regularization & numerical stability
- Added `add_ridge` argument to `fetwfe()` to include a small ridge penalty on coefficients
- Automatic centering and scaling of covariates before applying the ridge penalty

### Bug fixes & documentation updates
- Corrected mathematical definitions in the estimator  
- Fixed edge-case bugs (e.g., covariate centering/scaling)  
- Thoroughly updated Rd files, examples, and vignettes to illustrate all new features

## CRAN checks

- **R CMD check --as-cran**: no errors, warnings, or notes on macOS (Apple Silicon; M4 Pro)  
- **Windows devel** (via `devtools::check_win_devel()`): no errors, warnings, or notes  
- **Examples**: all run cleanly within recommended time limits

## Dependency changes

- No changes to `Imports` or `Suggests` since version 1.0.0.

## Other notes

- No compiled code changes or new system requirements.  
- All vignettes and example datasets have been rebuilt.  
- Maintains compatibility with R ≥ 4.1.0.

–––  
Please let me know if you need any additional information.  
Thank you for your time!  