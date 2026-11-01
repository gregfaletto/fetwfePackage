# cran-comments.md

**Package:** fetwfe
**Version:** 1.10.1
**Date:** 2026-11-01

This is a patch release responding to CRAN's report of 2026-10-31 that
fetwfe 1.10.0 fails the additional MKL check. It changes the tests and the
version bookkeeping only; there is no change to the package's code.

- Changed the tests that compared two separately computed results exactly to
  compare them at the default tolerance instead.

For both failing comparisons (`test-twfeCovs.R:1030` and `:1031` in 1.10.0),
CRAN's MKL log reads "Objects equal but not identical". The rewritten tests run
under testthat edition 2, where `expect_identical()` prints that message
exactly when the comparison `expect_equal()` makes finds the two objects
equal, so every rewritten comparison would pass on the values CRAN's run
produced. Later runs rest on a
local margin: on R's reference BLAS, perturbing one of the two compared
results by a relative 1e-14 fails every rewritten comparison as 1.10.0 wrote
it and none as 1.10.1 writes it, while a relative 1e-6 perturbation still
fails every rewritten comparison. Computing one side of each comparison with
R's own matrix products instead of the BLAS (`options(matprod = "internal")`)
fails 16 of the 24 rewritten comparisons as 1.10.0 wrote them and none as
1.10.1 writes them. This margin was measured locally; CRAN's MKL environment
was not reproduced.

## Test environments
- Local: macOS 26.6.2 (Apple Silicon), R 4.5.0, R's reference BLAS.
- R-hub, 2026-11-01, `R CMD check --as-cran --no-manual` with R-devel on each
  platform:
  - `mkl`: Fedora Linux 44, R-devel r90820, serial Intel MKL 2025.3
    (`/opt/intel/oneapi/mkl/2025.3/lib/libmkl_gf_lp64.so.2`);
  - `linux`: Ubuntu 24.04.5 LTS, R-devel r90805;
  - `windows`: Windows Server 2025 (build 26100), R-devel r90820 (ucrt);
  - `macos`: macOS Sequoia 15.7.9 (x86_64), R-devel r90820.

R-hub's MKL is 2025.3, not the 2026.2.0 of CRAN's check, and it did not
reproduce the failure: 1.10.0's own sources also pass on the `mkl` platform,
in a control run on the same image. So the R-hub runs are a smoke test on
R-devel, not evidence about MKL.

## R CMD check results
0 errors | 0 warnings | 0 notes, locally with `devtools::check()` and on each
R-hub platform. Neither builds the PDF manual or runs CRAN's incoming checks.

## Reverse dependencies
`fetwfe` has no reverse dependencies on CRAN.

–––
Thank you for your time reviewing this submission.
