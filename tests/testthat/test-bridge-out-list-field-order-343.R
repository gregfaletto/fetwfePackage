library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #343: pin the EXACT top-level and $internal field ORDER of the assembled
# betwfe()/fetwfe() fit objects. test-ols-out-list-field-order-270.R pins the OLS
# pair (etwfe/twfeCovs); the bridge pair had no order pin -- the class validators
# (test-class-validators.R) and doc-slot-parity test check the field SET, not the
# order. identical()-based consumers, print snapshots, and the planned
# .assemble_bridge_estimator() consolidation (#344) all depend on this order, so
# lock it against silent drift. (Note the two estimators legitimately differ:
# betwfe carries X_ints/y/X_final/y_final at the top level while fetwfe nests
# them in $internal and instead carries fusion_structure/fusion_matrix at the
# top; fetwfe also orders alpha before calc_ses, and its $internal adds theta_hat
# and d_inv_treat.)
# ------------------------------------------------------------------------------
test_that("betwfe()/fetwfe() top-level and internal field order is stable (#343)", {
	cf <- genCoefs(G = 4, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 270)
	sim <- simulateData(
		cf,
		N = 160,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 270
	)
	base <- list(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		verbose = FALSE
	)
	b <- do.call(betwfe, base)
	f <- do.call(fetwfe, base)

	betwfe_top <- c(
		"att_hat",
		"att_se",
		"att_p_value",
		"att_selected",
		"catt_hats",
		"catt_ses",
		"cohort_probs",
		"cohort_probs_overall",
		"catt_df",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"lambda_selection",
		"cv_folds",
		"cv_seed",
		"X_ints",
		"y",
		"X_final",
		"y_final",
		"N",
		"T",
		"G",
		"R",
		"d",
		"p",
		"calc_ses",
		"alpha",
		"se_type",
		"indep_counts_used",
		"y_mean",
		"response_col_name",
		"time_var",
		"unit_var",
		"treatment",
		"covs",
		"ci_type",
		"internal"
	)
	fetwfe_top <- c(
		"att_hat",
		"att_se",
		"att_p_value",
		"att_selected",
		"catt_hats",
		"catt_ses",
		"cohort_probs",
		"cohort_probs_overall",
		"catt_df",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"lambda_selection",
		"cv_folds",
		"cv_seed",
		"fusion_structure",
		"fusion_matrix",
		"N",
		"T",
		"G",
		"R",
		"d",
		"p",
		"alpha",
		"calc_ses",
		"se_type",
		"indep_counts_used",
		"y_mean",
		"response_col_name",
		"time_var",
		"unit_var",
		"treatment",
		"covs",
		"ci_type",
		"internal"
	)

	expect_identical(names(b), betwfe_top)
	expect_identical(names(f), fetwfe_top)

	expect_identical(
		names(b$internal),
		c(
			"X_ints",
			"y",
			"X_final",
			"y_final",
			"calc_ses",
			"variance_components",
			"first_year"
		)
	)
	expect_identical(
		names(f$internal),
		c(
			"X_ints",
			"y",
			"X_final",
			"y_final",
			"theta_hat",
			"calc_ses",
			"variance_components",
			"first_year",
			"d_inv_treat"
		)
	)
})
