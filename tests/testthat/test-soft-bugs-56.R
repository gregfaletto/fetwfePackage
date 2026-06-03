# Tests for the small-consistency cleanups bundled in issue #56 (v1.9.5).
#
# Items covered in this file:
#   1. R/event_study.R p-value mask: replace inline `is.na(ses) | ses == 0`
#      mask with a call to the canonical .compute_p_values() (the canonical
#      mask also handles negative SE — unreachable from package code paths
#      but a defensive guard against future refactors).
#   3. R/variance_machinery.R threshold: change `10e-6` (= 1e-5)
#      to `1e-6` to match the typed-intent literal. The change *loosens*
#      the rejection band from width 1e-5 to width 1e-6 — both are well
#      above floating-point roundoff and the guard is unreachable in
#      normal operation either way. Also adds the symmetric guard to
#      R/variance_machinery.R::getSecondVarTermDataApp() (which previously
#      had no such guard).
#   4. R/utility.R::idCohorts() panel-balance error message paren balance.
#
# Item 2 (`getSecondVarTermDataApp()` in R/variance_machinery.R, `<` ->
# `<=`) is not tested here. The assertion is structural; under
# `first_inds = getFirstInds(R, T)` block sizes are strictly decreasing
# so both forms pass. test-fetwfe-var2-fix.R already exercises the
# surrounding function on a standard fixture, so any regression in the
# relaxation would surface there.

# ------------------------------------------------------------------------------
# Item 1: .assemble_event_study_df() canonical NA mask for zero / negative SE.
# ------------------------------------------------------------------------------

test_that(".assemble_event_study_df masks zero and negative SEs to NA p-value", {
	out <- fetwfe:::.assemble_event_study_df(
		event_times = 0:3,
		n_cohorts = c(2L, 2L, 1L, 1L),
		estimates = c(0.5, -0.3, 0.2, 0.1),
		ses = c(0.2, 0, -0.1, 0.15),
		z = stats::qnorm(0.975)
	)
	# ses > 0: finite p-value.
	expect_false(is.na(out$p_value[1]))
	expect_false(is.na(out$p_value[4]))
	# ses == 0 and ses < 0: NA p-value. The negative-SE row is unreachable
	# from the package's own SE machinery (`sqrt(...)` is non-negative); the
	# test exercises the defensive guard.
	expect_true(is.na(out$p_value[2]))
	expect_true(is.na(out$p_value[3]))
})

# ------------------------------------------------------------------------------
# Item 3 — OLS path: getSecondVarTermOLS() asserts sum(cohort_probs_overall)
# < 1 - 1e-6 (the post-fix typed-intent literal).
# ------------------------------------------------------------------------------

.minimal_ols_call <- function(cohort_probs_overall) {
	# Minimal valid call shape for getSecondVarTermOLS: R=2, num_treats=3.
	# Post-#180 (DC2) the function no longer accepts `first_inds` — the
	# pre-#180 helper passed `first_inds = c(1L, 3L)` only to feed a
	# tautological `stopifnot(all.equal(unlist(sel_inds), 1:num_treats))`
	# assertion inside the function, which was dropped.
	fetwfe:::getSecondVarTermOLS(
		psi_mat = matrix(c(0.5, 0.5, 0, 0, 0, 1), nrow = 3, ncol = 2),
		tes = c(1, 2, 3),
		cohort_probs_overall = cohort_probs_overall,
		num_treats = 3L,
		N = 100L,
		T = 3L,
		G = 2L
	)
}

test_that("getSecondVarTermOLS fires the post-fix guard for sums in the new rejection band", {
	# sum = 1 - 5e-7 is inside the new [1 - 1e-6, 1) rejection band.
	cpo <- c(0.5, 0.5 - 5e-7)
	stopifnot(abs(sum(cpo) - (1 - 5e-7)) < 1e-12)
	expect_error(
		.minimal_ols_call(cpo),
		"1 - 1e-06",
		fixed = TRUE
	)
})

test_that("getSecondVarTermOLS does NOT fire on sums inside the old (loose) band but outside the new (loose) one", {
	# sum = 1 - 5e-6 was inside the OLD rejection band (10e-6 = 1e-5);
	# under the new threshold it is NOT in the [1 - 1e-6, 1) band, so the
	# guard must NOT fire. This documents the loosening explicitly so a
	# future reverter sees the looser behavior is intended.
	cpo <- c(0.5, 0.5 - 5e-6)
	stopifnot(abs(sum(cpo) - (1 - 5e-6)) < 1e-12)
	expect_silent(.minimal_ols_call(cpo))
})

# ------------------------------------------------------------------------------
# Item 3 — FETWFE path: getSecondVarTermDataApp() gains the same guard.
# ------------------------------------------------------------------------------

.minimal_fetwfe_call <- function(cohort_probs_overall) {
	# Minimal valid call shape for getSecondVarTermDataApp: R=2, num_treats=3,
	# first_inds = c(1, 3) so sel_inds = list(1:2, 3:3) and unlist = 1:3.
	# d_inv_treat_sel is the inverse fusion-transform matrix (3x3 for these
	# dimensions); we build it via the package's helper so it has the right
	# structure.
	d_inv <- fetwfe:::genInvTwoWayFusionTransformMat(
		n_vars = 3L,
		first_inds = c(1L, 3L),
		G = 2L
	)
	fetwfe:::getSecondVarTermDataApp(
		sel_treat_inds_shifted = 1:3,
		cohort_probs_overall = cohort_probs_overall,
		first_inds = c(1L, 3L),
		theta_hat_treat_sel = c(0.1, 0.2, 0.3),
		num_treats = 3L,
		N = 100L,
		T = 3L,
		G = 2L,
		d_inv_treat_sel = d_inv
	)
}

test_that("getSecondVarTermDataApp fires the new symmetric guard for sums in the new rejection band", {
	cpo <- c(0.5, 0.5 - 5e-7)
	expect_error(
		.minimal_fetwfe_call(cpo),
		"1 - 1e-06",
		fixed = TRUE
	)
})

test_that("getSecondVarTermDataApp does NOT fire for sums outside the new rejection band", {
	cpo <- c(0.5, 0.5 - 5e-6)
	expect_silent(.minimal_fetwfe_call(cpo))
})

# ------------------------------------------------------------------------------
# Item 4: idCohorts() unbalanced-panel error message has balanced parens.
# ------------------------------------------------------------------------------

test_that("idCohorts unbalanced-panel error message has balanced parens", {
	df <- data.frame(
		unit = c("A", "A", "A", "B", "B"),
		time = c(1, 2, 3, 1, 2),
		treat = c(0, 0, 1, 0, 1)
	)
	err <- tryCatch(
		fetwfe:::idCohorts(
			df = df,
			time_var = "time",
			unit_var = "unit",
			treat_var = "treat"
		),
		error = function(e) conditionMessage(e)
	)
	# Sanity: the canonical message text is preserved.
	expect_match(err, "Panel does not appear to be balanced")
	# Parens balance. Counting via nchar(gsub(non-paren-chars-removed))
	# avoids gregexpr's "-1 with length 1 on no match" gotcha.
	expect_equal(
		nchar(gsub("[^(]", "", err)),
		nchar(gsub("[^)]", "", err))
	)
})
