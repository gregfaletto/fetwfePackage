library(testthat)
library(fetwfe)

# Issue #174, anti-silent-swallow guard.
#
# Predecessor PR #156 (closing #138) added an "Event Study (preview):"
# section to `print.<class>()` and `summary.<class>()` by wrapping the
# `eventStudy(x)` call in `tryCatch(..., error = function(e) NULL)` and
# silently omitting the section on error. That swallow masked the
# scattered-cohort bug on `bacondecomp::divorce` for two releases (the
# eventStudy() D-inverse-matrix construction errored, the tryCatch ate
# it, the user saw nothing where the section should have been). #174
# removed the swallow under the strict policy
#
#     `eventStudy()`'s contract is "succeeds on any fit produced by
#     `fetwfe()` / `betwfe()` / `etwfe()` / `twfeCovs()`". If it fails
#     on a valid fit, that is always a bug.
#
# This test file guards that contract: a real `eventStudy()` failure on
# a (deliberately-mutated) fit must propagate through `print(res)` and
# `summary(res)` rather than being silently swallowed. If a future PR
# re-introduces a defense-in-depth `tryCatch(eventStudy(x), ...)` at
# the print/summary layer, these tests fail.

.no_swallow_setup <- function() {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	simulateData(coefs, N = 60, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

# Failure mode: blank-out the cohort_probs_overall slot. The
# `.event_study_etwfe_betwfe()` body explicitly errors with a stable
# message when `cohort_probs_overall` is NULL (see R/event_study.R
# `eventStudy(): cohort_probs_overall missing from object.`); the
# fetwfe path has a parallel guard on `internal$theta_hat`. Both are
# stable substrings to match against.

# Implementation note: under #174 `eventStudy()`'s preconditon
# (`.check_for_event_study()`) runs `.assert_estimator_object()`, which
# in turn runs `.validate_<class>()`. The validator catches the
# missing-slot mutation BEFORE the eventStudy body's own `stop()` on the
# same condition fires. Either origin proves the silent-swallow
# tryCatch is gone (the prior tryCatch swallowed both); the test
# matches on the validator's stable wording (`Missing slot(s):`).

test_that("print(res) propagates eventStudy() error for ETWFE (no silent swallow)", {
	sim <- .no_swallow_setup()
	res <- etwfeWithSimulatedData(sim)
	# The print path will hit `eventStudy()` ->
	# `.check_for_event_study()` -> `.validate_etwfe()` -> stop with
	# `Missing slot(s): cohort_probs_overall`.
	res$cohort_probs_overall <- NULL
	expect_error(
		print(res),
		"Missing slot.s.: cohort_probs_overall"
	)
})

test_that("summary(res) propagates eventStudy() error for ETWFE (no silent swallow)", {
	sim <- .no_swallow_setup()
	res <- etwfeWithSimulatedData(sim)
	res$cohort_probs_overall <- NULL
	expect_error(
		summary(res),
		"Missing slot.s.: cohort_probs_overall"
	)
})

test_that("print(res) propagates eventStudy() error for FETWFE (no silent swallow)", {
	sim <- .no_swallow_setup()
	res <- fetwfeWithSimulatedData(sim, q = 0.5)
	# FETWFE's `internal$theta_hat` is in
	# `.EXPECTED_INTERNAL_SLOTS_FETWFE`, so the validator flags it.
	res$internal$theta_hat <- NULL
	expect_error(
		print(res),
		"Missing slot.s.: theta_hat"
	)
})

test_that("summary(res) propagates eventStudy() error for FETWFE (no silent swallow)", {
	sim <- .no_swallow_setup()
	res <- fetwfeWithSimulatedData(sim, q = 0.5)
	res$internal$theta_hat <- NULL
	expect_error(
		summary(res),
		"Missing slot.s.: theta_hat"
	)
})

test_that("print(res) propagates eventStudy() error for BETWFE (no silent swallow)", {
	sim <- .no_swallow_setup()
	res <- betwfeWithSimulatedData(sim, q = 0.5)
	# Same trigger as ETWFE: validator stops on the missing slot.
	res$cohort_probs_overall <- NULL
	expect_error(
		print(res),
		"Missing slot.s.: cohort_probs_overall"
	)
})
