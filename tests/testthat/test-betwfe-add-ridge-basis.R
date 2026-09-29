library(testthat)
library(fetwfe)

# Tests for the BETWFE add_ridge wrong-basis fix (issue #74, v1.9.12).
#
# Background: before this fix, `betwfe_core()`'s call to
# `prep_for_etwfe_regression()` did not pass `is_fetwfe`, so it picked up
# the default `is_fetwfe = TRUE` (in `prep_for_etwfe_regression()`, `R/input_prep.R`). With `add_ridge
# = TRUE`, this caused BETWFE's ridge augmentation rows to be built as
# `sqrt(lambda_ridge) * D_inverse` (the inverse FETWFE fusion-transform
# matrix) instead of the correct `sqrt(lambda_ridge) * diag(p)` (identity
# basis). Silent — the fit converged, but with the wrong L2 penalty
# structure.
#
# This file's contracts:
#   1. `prep_for_etwfe_regression(is_fetwfe = FALSE)` augments with
#      identity basis (the BETWFE / ETWFE / twfeCovs path).
#   2. `prep_for_etwfe_regression(is_fetwfe = TRUE)` augments with
#      D_inverse basis (the FETWFE path).
#   3. The helper contract of 1, checked on a `betwfe(add_ridge = TRUE)`
#      fit's own inputs.
#   4. The #74 guard: a recorder around `.append_ridge_rows()` asserts that
#      a `betwfe(add_ridge = TRUE)` fit passes it `is_fetwfe = FALSE`. So it
#      fails if `betwfe_core()` passes `is_fetwfe = TRUE`, and errors if
#      `betwfe_core()` omits the argument, which then has no value to record.

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
	# and produces finite output.
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

	# The helper contract on this fit's own inputs: a helper call with
	# is_fetwfe = FALSE builds identity-basis augmentation rows.
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

test_that("betwfe(add_ridge = TRUE) builds its ridge rows in the identity basis (recorder)", {
	real <- fetwfe:::.append_ridge_rows
	rec <- new.env(parent = emptyenv())
	rec$is_fetwfe <- list()
	testthat::with_mocked_bindings(
		betwfe(
			pdata = pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2"),
			add_ridge = TRUE,
			verbose = FALSE
		),
		.append_ridge_rows = function(...) {
			a <- list(...)
			rec$is_fetwfe[[length(rec$is_fetwfe) + 1L]] <- a$is_fetwfe
			real(...)
		},
		.package = "fetwfe"
	)
	expect_identical(rec$is_fetwfe, list(FALSE))
})

test_that("betwfe(add_ridge = TRUE) att_hat pin", {
	# A plain pin of one bridge fit with the ridge applied.
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
	# under the BIC selection path. The default is CV, which produces a
	# slightly different att_hat on this fixture; that's expected, and is
	# tested separately in test-lambda-selection-164.R.
	fit <- betwfeWithSimulatedData(
		sim,
		add_ridge = TRUE,
		lambda_selection = "bic"
	)

	# The fit selected something, so the pin is not on a null model.
	expect_true(any(fit$beta_hat != 0))

	# Tolerance 1e-5 comfortably exceeds cross-platform float drift on
	# grpreg + BLAS.
	expect_equal(fit$att_hat, 3.4803274345, tolerance = 1e-5)
})
