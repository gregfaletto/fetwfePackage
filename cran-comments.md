# cran-comments.md

**Package:** fetwfe  
**Version:** 1.0.0 
**Date:** 2025-05-14

## What’s new

- **New exported functions**  
  - `simulateData()` for simulating panel data from a `FETWFE_coefs` object  
  - `genCoefs()` to generate coefficient vectors for simulation studies  
  - `fetwfeWithSimulatedData()` wrapper to run `fetwfe()` on simulated data  
  - `getTes()` to compute true treatment effects from a `FETWFE_coefs` object
  - `simulateDataCore()`, called under the hood by `simulateData()` and exported for convenience
  - `genCoefsCore()`, called under the hood by `genCoefs()` and exported for convenience`
- **Enhancements to `fetwfe()`**  
  - Centering and scaling of covariates before ridge penalty for improved numerical stability  
  - Support for factor covariates (`covs` now accepts factors, automatically encoded)  
  - New `add_ridge` argument to include a small ridge regularization term  
- **Bug fixes & refactoring**  
  - Fixed centering/scaling edge-case (issue #1)  
  - Various minor bug fixes and code refactoring for clarity and maintainability  
- **Documentation**  
  - Updated Rd documentation and examples for all new functions  
  - Enhanced vignette with simulation workflows and examples using `simulateData()` and `getTes()`

## CRAN checks

- **R CMD check --as-cran**: no errors, warnings, or notes on macOS (Apple Silicon; M4 Pro)
- **Windows devel** (via `devtools::check_win_devel()`): no errors, warnings, or notes
- **Examples**: all examples run cleanly within time limits

## Dependency changes

- No changes to `Imports` or `Suggests` since version 0.4.4

## Other notes

- No compiled code changes or new system requirements  
- Vignettes and example datasets regenerated  
- Maintains compatibility with R ≥ 3.6

–––  
Please let me know if you need any additional information.  
Thank you for your time!