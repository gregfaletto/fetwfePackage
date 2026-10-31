# cran-comments.md

**Package:** fetwfe
**Version:** 1.10.1
**Date:** 2026-10-31

This is a patch release responding to CRAN's report of 2026-10-31 that
fetwfe 1.10.0 fails the additional MKL check. It changes the tests and the
version bookkeeping only; there is no change to the package's code.

- Changed the tests that compared two separately computed results exactly to
  compare them at the default tolerance instead.

For both failing comparisons (`test-twfeCovs.R:1030` and `:1031` in 1.10.0),
CRAN's MKL log reads "Objects equal but not identical". The tests run under
testthat edition 2, where `expect_identical()` prints that message exactly
when the comparison `expect_equal()` makes finds the two objects equal, so on
the run CRAN observed every rewritten comparison passes. Later runs rest on a
local margin: on R's reference BLAS, perturbing one of the two compared
results by a relative 1e-14 fails every rewritten comparison as 1.10.0 wrote
it and none as 1.10.1 writes it, while a relative 1e-6 perturbation still
fails every rewritten comparison. This margin was measured locally; CRAN's MKL
environment was not reproduced.

## Test environments
- Local: macOS 26.6.2 (Apple Silicon), R 4.5.0, R's reference BLAS.

## R CMD check results
0 errors | 0 warnings | 0 notes.

## Reverse dependencies
`fetwfe` has no reverse dependencies on CRAN.

–––
Thank you for your time reviewing this submission.
