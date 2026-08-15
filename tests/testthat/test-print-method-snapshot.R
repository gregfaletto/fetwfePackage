library(testthat)
library(fetwfe)

# Snapshot-based guardrail for `print.<class>(x)` and `print.summary.<class>(x)`
# across all four estimators (`fetwfe`, `etwfe`, `betwfe`, `twfeCovs`). Locks
# the current printed output of each method against ten markdown goldens under
# tests/testthat/_snaps/print-method-snapshot.md so any byte drift surfaces as
# a failing test. Issue #77 step 1 of 2 (the follow-up consolidates the three
# `R/*_class.R` print/summary bodies into a shared helper; this guardrail
# protects that refactor). Extended by #439, which consolidated five rendering
# pieces shared between the print and summary bodies.
#
# TWO fixtures, deliberately. Do not merge them.
#
# `pdata` (N=30, T=5, R=2, seed=123) is the only parameter set all four
# estimators can fit without error: `etwfe()` rejects larger (T, R) due to
# pure-OLS rank deficiency (the bridge regression in `fetwfe()`/`betwfe()`
# absorbs it). But it yields G == d == 2, and the Model Details block prints
# `Treated cohorts (G)` and `Covariates (d)` as `%d` on same-typed integers --
# so on this fixture alone a G/d transposition renders identically and is
# invisible to the entire suite. Measured during #439's plan review: swapping
# those two rows left the suite at FAIL 0 / PASS 526.
#
# `pdata_distinct` (N=40, T=6, R=3, seed=123, one covariate) exists solely to
# break that symmetry: N=40, T=6, G=3, d=1, p=41 are five distinct values, so
# any transposition among them changes the rendered output. It is fetwfe-only,
# which is fine -- `.cat_model_details()` is class-agnostic.
#
# `generate_panel_data()` itself lives in `tests/testthat/helper-panel-fixture.R`,
# which testthat sources before any test file (#91 consolidated what had been 14
# inline copies). This comment previously said the body was duplicated verbatim
# here; that stopped being true at #91.
#
# Snapshot tests require testthat edition 3 and skip under CRAN by default
# (so `R CMD check --as-cran` shows SKIPs for these tests; `devtools::test()`
# runs them).

testthat::local_edition(3)

# generate_panel_data() is defined in tests/testthat/helper-panel-fixture.R
# (sourced by testthat before this file runs; issue #91).

pdata <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

fit_fetwfe <- fetwfe(
	pdata = pdata,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = c("cov1", "cov2"),
	verbose = FALSE,
	# Pin to BIC so the snapshot is stable across the v1.13.0 default
	# change. The CV path's stochastic fold assignment would produce a
	# different `Lambda*` in the printed output (only deterministic given
	# the seed, but locking the seed here couples the snapshot to the
	# CV-implementation details rather than to the print-method format,
	# which is what this guardrail is for).
	lambda_selection = "bic"
)

fit_etwfe <- etwfe(
	pdata = pdata,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = c("cov1", "cov2"),
	verbose = FALSE
)

fit_betwfe <- betwfe(
	pdata = pdata,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = c("cov1", "cov2"),
	verbose = FALSE,
	# Pin to BIC for the same reason as fit_fetwfe above.
	lambda_selection = "bic"
)

# twfeCovs (#439). The fourth estimator class had no golden at all, and it is
# the ONLY caller that passes `include_event_study = FALSE` -- so it is the one
# class whose output must contain no event-study section, and the only golden
# anywhere covering that branch of the two refactored functions.
fit_twfecovs <- twfeCovs(
	pdata = pdata,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = c("cov1", "cov2"),
	verbose = FALSE
)

# The distinct-dimension fixture (#439). See the header comment for why this is
# a separate fixture rather than a tweak to `pdata`. `covs = "cov1"` (one
# covariate, not two) is what makes d = 1 rather than 2.
pdata_distinct <- generate_panel_data(N = 40, T = 6, R = 3, seed = 123)

fit_distinct <- fetwfe(
	pdata = pdata_distinct,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = "cov1",
	verbose = FALSE,
	# Pin to BIC for the same reason as fit_fetwfe above. This call must stay
	# byte-identical to `fit_distinct` in
	# tests/testthat/test-print-summary-single-source-439.R, which asserts the
	# same Model Details lines these goldens lock.
	lambda_selection = "bic"
)

test_that("print.fetwfe output is stable", {
	expect_snapshot(print(fit_fetwfe))
})

test_that("print.summary.fetwfe output is stable", {
	expect_snapshot(print(summary(fit_fetwfe)))
})

test_that("print.etwfe output is stable", {
	expect_snapshot(print(fit_etwfe))
})

test_that("print.summary.etwfe output is stable", {
	expect_snapshot(print(summary(fit_etwfe)))
})

test_that("print.betwfe output is stable", {
	expect_snapshot(print(fit_betwfe))
})

test_that("print.summary.betwfe output is stable", {
	expect_snapshot(print(summary(fit_betwfe)))
})

# ------------------------------------------------------------------------------
# #439: the fourth estimator class, and the distinct-dimension fixture.
# ------------------------------------------------------------------------------

test_that("print.twfeCovs output is stable", {
	expect_snapshot(print(fit_twfecovs))
})

test_that("print.summary.twfeCovs output is stable", {
	expect_snapshot(print(summary(fit_twfecovs)))
})

test_that("print.fetwfe output is stable on a distinct-dimension fit", {
	expect_snapshot(print(fit_distinct))
})

test_that("print.summary.fetwfe output is stable on a distinct-dimension fit", {
	expect_snapshot(print(summary(fit_distinct)))
})
