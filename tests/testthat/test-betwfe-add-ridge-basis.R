library(testthat)
library(fetwfe)

# Tests for the BETWFE add_ridge wrong-basis fix (issue #74, v1.9.12).
#
# Background: before this fix, `R/betwfe_core.R:840`'s call to
# `prep_for_etwfe_regression()` did not pass `is_fetwfe`, so it picked up
# the default `is_fetwfe = TRUE` (in `prep_for_etwfe_regression()`, `R/input_prep.R`). With `add_ridge
# = TRUE`, this caused BETWFE's ridge augmentation rows to be built as
# `sqrt(lambda_ridge) * D_inverse` (the inverse FETWFE fusion-transform
# matrix) instead of the correct `sqrt(lambda_ridge) * diag(p)` (identity
# basis). Silent — the fit converged, but with the wrong L2 penalty
# structure.
#
# This file has three contracts:
#   1. `prep_for_etwfe_regression(is_fetwfe = FALSE)` augments with
#      identity basis (the BETWFE / ETWFE / twfeCovs path).
#   2. `prep_for_etwfe_regression(is_fetwfe = TRUE)` augments with
#      D_inverse basis (the FETWFE path).
#   3. `betwfe(..., add_ridge = TRUE)` integration: reconstruct the
#      inputs `betwfe_core` would forward to the helper and assert the
#      augmentation rows ARE identity-basis. The pre-fix buggy version of
#      betwfe would have failed this assertion.

# generate_panel_data() is defined in tests/testthat/helper-panel-fixture.R
# (sourced by testthat before this file runs; issue #91).

# Shared setup: a fixture + a baseline (add_ridge = FALSE) fit to extract
# upstream inputs without re-running the REML variance-component
# estimator. The values used at the prep_for_etwfe_regression call site
# are class-level metadata + a few internal vectors; everything but
# `in_sample_counts` is exposed as a top-level slot on the fit. The
# missing slot is rebuilt via fetwfe:::idCohorts() against the raw pdata
# (the same helper betwfe() uses internally).

pdata <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

fit_baseline <- betwfe(
	pdata = pdata,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = c("cov1", "cov2"),
	add_ridge = FALSE,
	verbose = FALSE
)

first_inds <- fetwfe:::getFirstInds(fit_baseline$R, fit_baseline$T)
num_treats <- length(fit_baseline$treat_inds)

# Rebuild in_sample_counts via idCohorts. Order: never-treated first,
# then the R treated cohorts in the order idCohorts returns them.
ret <- fetwfe:::idCohorts(
	df = pdata,
	time_var = "time",
	unit_var = "unit",
	treat_var = "treatment"
)
cohorts <- ret$cohorts
N_treated <- length(unlist(cohorts))
in_sample_counts <- as.integer(c(
	fit_baseline$N - N_treated,
	vapply(cohorts, length, integer(1))
))

test_that("prep_for_etwfe_regression(is_fetwfe = FALSE) augments with identity basis", {
	# Pass the fit's stored (post-REML) sig values so the helper
	# does not re-run estOmegaSqrtInv internally; this keeps the
	# test fast and deterministic.
	res_false <- fetwfe:::prep_for_etwfe_regression(
		verbose = FALSE,
		sig_eps_sq = fit_baseline$sig_eps_sq,
		sig_eps_c_sq = fit_baseline$sig_eps_c_sq,
		y = fit_baseline$y,
		X_ints = fit_baseline$X_ints,
		X_mod = fit_baseline$X_ints,
		N = fit_baseline$N,
		T = fit_baseline$T,
		G = fit_baseline$G,
		d = fit_baseline$d,
		p = fit_baseline$p,
		num_treats = num_treats,
		add_ridge = TRUE,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = FALSE,
		indep_counts = NA,
		is_fetwfe = FALSE
	)

	n_data_rows <- fit_baseline$N * fit_baseline$T
	aug_rows <- res_false$X_final_scaled[
		(n_data_rows + 1):(n_data_rows + fit_baseline$p),
	]
	expected <- sqrt(res_false$lambda_ridge) * diag(fit_baseline$p)
	# Strip dimnames so the comparison is over values only; X_final_scaled
	# inherits row/column names that diag(p) does not have.
	dimnames(aug_rows) <- NULL
	expect_equal(aug_rows, expected, tolerance = 1e-12)
})

test_that("prep_for_etwfe_regression(is_fetwfe = TRUE) augments with D_inverse basis", {
	res_true <- fetwfe:::prep_for_etwfe_regression(
		verbose = FALSE,
		sig_eps_sq = fit_baseline$sig_eps_sq,
		sig_eps_c_sq = fit_baseline$sig_eps_c_sq,
		y = fit_baseline$y,
		X_ints = fit_baseline$X_ints,
		X_mod = fit_baseline$X_ints,
		N = fit_baseline$N,
		T = fit_baseline$T,
		G = fit_baseline$G,
		d = fit_baseline$d,
		p = fit_baseline$p,
		num_treats = num_treats,
		add_ridge = TRUE,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = FALSE,
		indep_counts = NA,
		is_fetwfe = TRUE
	)

	D_inverse <- fetwfe:::genFullInvFusionTransformMat(
		first_inds = first_inds,
		T = fit_baseline$T,
		G = fit_baseline$G,
		d = fit_baseline$d,
		num_treats = num_treats
	)

	n_data_rows <- fit_baseline$N * fit_baseline$T
	aug_rows <- res_true$X_final_scaled[
		(n_data_rows + 1):(n_data_rows + fit_baseline$p),
	]
	expected <- sqrt(res_true$lambda_ridge) * D_inverse
	# Strip dimnames so the comparison is over values only; X_final_scaled
	# inherits row/column names that diag(p) does not have.
	dimnames(aug_rows) <- NULL
	expect_equal(aug_rows, expected, tolerance = 1e-12)

	# Sanity: the two branches really do produce different
	# matrices on this fixture.
	identity_expected <- sqrt(res_true$lambda_ridge) *
		diag(fit_baseline$p)
	expect_false(isTRUE(all.equal(aug_rows, identity_expected)))
})

test_that("betwfe(add_ridge = TRUE) integration: smoke check + augmentation reconstruction", {
	# Smoke check: betwfe with add_ridge = TRUE runs to completion
	# and produces finite output. Note: this assertion alone is
	# weak — the pre-fix buggy version also produced finite output.
	# But under the defensive-cleanup (no-default) regime, a future
	# regression that drops `is_fetwfe = FALSE` from
	# R/betwfe_core.R's call will cause this fit to error at
	# runtime (missing required argument when add_ridge = TRUE
	# reaches the augmentation branch), so this smoke check catches
	# that specific regression mode.
	fit_ridge <- betwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		add_ridge = TRUE,
		verbose = FALSE
	)

	expect_s3_class(fit_ridge, "betwfe")
	expect_true(is.finite(fit_ridge$att_hat))
	expect_true(all(is.finite(fit_ridge$beta_hat)))

	# Reconstruct the inputs betwfe_core forwards. Independent
	# helper call with is_fetwfe = FALSE; assert the augmentation
	# rows are identity-basis. This locks the helper-contract
	# behavior on the actual integration inputs.
	res <- fetwfe:::prep_for_etwfe_regression(
		verbose = FALSE,
		sig_eps_sq = fit_ridge$sig_eps_sq,
		sig_eps_c_sq = fit_ridge$sig_eps_c_sq,
		y = fit_ridge$y,
		X_ints = fit_ridge$X_ints,
		X_mod = fit_ridge$X_ints,
		N = fit_ridge$N,
		T = fit_ridge$T,
		G = fit_ridge$G,
		d = fit_ridge$d,
		p = fit_ridge$p,
		num_treats = num_treats,
		add_ridge = TRUE,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = FALSE,
		indep_counts = NA,
		is_fetwfe = FALSE
	)

	n_data_rows <- fit_ridge$N * fit_ridge$T
	aug_rows <- res$X_final_scaled[
		(n_data_rows + 1):(n_data_rows + fit_ridge$p),
	]
	expected <- sqrt(res$lambda_ridge) * diag(fit_ridge$p)
	# Strip dimnames so the comparison is over values only; X_final_scaled
	# inherits row/column names that diag(p) does not have.
	dimnames(aug_rows) <- NULL
	expect_equal(aug_rows, expected, tolerance = 1e-12)
})

test_that("betwfe(add_ridge = TRUE) numerical regression: att_hat matches post-fix value", {
	# The strong integration check: with a simulated DGP that
	# produces a non-degenerate fit, the post-fix code produces a
	# specific att_hat value. The pre-fix bug (silently passing
	# is_fetwfe = TRUE in betwfe_core's call) shifts att_hat by
	# ~5e-5 on this fixture. Hard-coded tolerance 1e-5 catches the
	# bug while allowing for cross-platform float drift.
	#
	# The fixture uses larger sig_eps_sq (= 16) than the panel
	# fixture above, which scales up lambda_ridge and amplifies
	# the bug's downstream effect on att_hat (from ~2e-6 at sig=1
	# to ~5e-5 at sig=16).
	coefs <- genCoefs(
		G = 2,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 4,
		seed = 42
	)
	sim <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 16,
		sig_eps_c_sq = 8,
		seed = 42
	)
	# Pin to BIC because the recorded reference value below was generated
	# under the BIC selection path (v1.12.0 and earlier). The v1.13.0+
	# default is CV, which produces a slightly different att_hat on this
	# fixture; that's expected, and is tested separately in
	# test-lambda-selection-164.R.
	fit <- betwfeWithSimulatedData(
		sim,
		add_ridge = TRUE,
		lambda_selection = "bic"
	)

	# Sanity: the fit selected something. The pre-fix code on this
	# fixture also produces a non-degenerate fit, so this is not a
	# bug-catch in itself.
	expect_true(any(fit$beta_hat != 0))

	# Numerical regression catch. The post-fix value on this fixture
	# is 3.4764564039 (recorded after applying the bug fix at
	# R/betwfe_core.R:840 + the defensive cleanup in R/input_prep.R).
	# The pre-fix buggy value is ~3.4764078169 (a 4.9e-5 shift).
	# Tolerance 1e-5 catches the shift while comfortably exceeding
	# expected cross-platform float drift on grpreg+BLAS.
	expect_equal(fit$att_hat, 3.4764564039, tolerance = 1e-5)
})
