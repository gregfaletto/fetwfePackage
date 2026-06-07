library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #270: etwfe() and twfeCovs() now share .assemble_ols_estimator() for their entry
# bodies. The two estimators' top-level out-lists have the SAME 35 fields but a
# DIFFERENT historical order at positions 24-26 (etwfe: alpha, calc_ses, se_type;
# twfeCovs: calc_ses, se_type, alpha). identical()-based consumers and backward
# compatibility depend on these orders, and the shared assembler reorders by name
# for twfeCovs -- so pin both orders against a silent regression.
# ------------------------------------------------------------------------------
test_that("etwfe()/twfeCovs() top-level field order is preserved by the shared assembler (#270)", {
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
	e <- do.call(etwfe, base)
	t <- do.call(twfeCovs, base)

	etwfe_order <- c(
		"att_hat",
		"att_se",
		"att_p_value",
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
		"alpha",
		"calc_ses",
		"se_type", # positions 24-26 (etwfe order)
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
	# twfeCovs differs from etwfe ONLY at positions 24-26.
	twfecovs_order <- etwfe_order
	twfecovs_order[24:26] <- c("calc_ses", "se_type", "alpha")

	expect_identical(names(e), etwfe_order)
	expect_identical(names(t), twfecovs_order)
	# Same field SET, differing only in the alpha/calc_ses/se_type block.
	expect_setequal(names(e), names(t))
})
