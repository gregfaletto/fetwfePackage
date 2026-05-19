library(testthat)
library(fetwfe)

# Snapshot-based guardrail for `print.<class>(x)` and `print.summary.<class>(x)`
# across the three estimators (`fetwfe`, `etwfe`, `betwfe`). Locks the current
# printed output of each method against six markdown goldens under
# tests/testthat/_snaps/print-method-snapshot.md so any byte drift surfaces as
# a failing test. Issue #77 step 1 of 2 (the follow-up consolidates the three
# `R/*_class.R` print/summary bodies into a shared helper; this guardrail
# protects that refactor).
#
# Fixture parameters (N=30, T=5, R=2, seed=123) are the only ones all three
# estimators can fit without error: `etwfe()` rejects larger (T, R) due to
# pure-OLS rank deficiency (the bridge regression in `fetwfe()`/`betwfe()`
# absorbs it). The `generate_panel_data()` body is verbatim from
# `tests/testthat/test-fetwfe.R`; the duplication is acceptable so the
# guardrail file is self-contained (see .plans/feat-print-snapshot-77/PLAN.md
# Decision Log D3 / D4).
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
	verbose = FALSE
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
	verbose = FALSE
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
