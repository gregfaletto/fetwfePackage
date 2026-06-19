# Tests for `.build_selected_out_result()` in R/result_assembly.R (issue #80).
#
# The helper consolidates the four pre-1.9.14 inline "no features selected"
# early-exit blocks that lived in `R/betwfe_core.R` (betwfe_core(), two
# blocks) and `R/fetwfe_core.R` (fetwfe_core(), two blocks). Per the plan
# (`.plans/refactor-selected-out-helper-80/PLAN.md`), the 4 blocks share
# the same 33-/34-field return list and only differ in (a) the `beta_hat`
# value, (b) whether `theta_hat` is included, and (c) the verbose-message
# text. These tests lock the helper's contract directly with hand-crafted
# inputs rather than constructing fixtures that trigger each early-exit
# path from an end-to-end fit (much simpler; the existing test suite
# already exercises the integration).
#
# Snapshot guardrail covers Blocks 1 + 3 (BETWFE and FETWFE intercept-only)
# via `tests/testthat/_snaps/print-method-snapshot.md`; these tests cover
# Blocks 2 + 4 (no-treatment) and lock the helper's contract directly.

# Builder for the minimal helper-input keyword-arg list. Returns a list
# with the union of all arguments at deterministic values. Callers spread
# this via `do.call(fetwfe:::.build_selected_out_result, ...)` and override
# specific fields per test.
.make_helper_args <- function() {
	G <- 2L
	N <- 10L
	T <- 3L
	d <- 0L
	p <- 7L # contrived but consistent; helper doesn't validate p arithmetic
	list(
		message_text = "test message",
		verbose = FALSE,
		G = G,
		c_names = c("2", "3"),
		q = 0.5,
		beta_hat = rep(0, p),
		treat_inds = 3:5,
		treat_int_inds = integer(0),
		cohort_probs = c(0.5, 0.5),
		indep_cohort_probs = c(0.5, 0.5),
		cohort_probs_overall = 0.7,
		indep_cohort_probs_overall = 0.7,
		sig_eps_sq = 1.0,
		sig_eps_c_sq = 0.5,
		lambda.max = 0.10,
		lambda.max_model_size = 5L,
		lambda.min = 0.001,
		lambda.min_model_size = 20L,
		lambda_star = 0.05,
		lambda_star_model_size = 1L,
		X_ints = matrix(0, N * T, p),
		y = rep(0, N * T),
		X_final = matrix(0, N * T, p),
		y_final = rep(0, N * T),
		N = N,
		T = T,
		d = d,
		p = p
	)
}

# Field ordering the helper must produce for BETWFE-shape returns (Blocks
# 1 + 2). 38 fields after v1.13.0 (#164) added `cv_seed_used` for
# CV-path provenance (37 in v1.12.0 after #141/#146 added 4
# variance-component slots). No `theta_hat`.
.expected_betwfe_field_order <- c(
	"in_sample_att_hat",
	"in_sample_att_se",
	"in_sample_att_se_no_prob",
	"in_sample_att_var_1",
	"in_sample_att_var_2",
	"indep_att_hat",
	"indep_att_se",
	"indep_att_var_1",
	"indep_att_var_2",
	"catt_hats",
	"catt_ses",
	"catt_df",
	"beta_hat",
	"treat_inds",
	"treat_int_inds",
	"cohort_probs",
	"indep_cohort_probs",
	"cohort_probs_overall",
	"indep_cohort_probs_overall",
	"sig_eps_sq",
	"sig_eps_c_sq",
	"lambda.max",
	"lambda.max_model_size",
	"lambda.min",
	"lambda.min_model_size",
	"lambda_star",
	"lambda_star_model_size",
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"N",
	"T",
	"G",
	"d",
	"p",
	"calc_ses",
	"cv_seed_used"
)

# Field ordering the helper must produce for FETWFE-shape returns (Blocks
# 3 + 4). 39 fields after v1.13.0 added `cv_seed_used` (38 in v1.12.0).
# `theta_hat` slots between `catt_df` and `beta_hat`.
.expected_fetwfe_field_order <- c(
	"in_sample_att_hat",
	"in_sample_att_se",
	"in_sample_att_se_no_prob",
	"in_sample_att_var_1",
	"in_sample_att_var_2",
	"indep_att_hat",
	"indep_att_se",
	"indep_att_var_1",
	"indep_att_var_2",
	"catt_hats",
	"catt_ses",
	"catt_df",
	"theta_hat",
	"beta_hat",
	"treat_inds",
	"treat_int_inds",
	"cohort_probs",
	"indep_cohort_probs",
	"cohort_probs_overall",
	"indep_cohort_probs_overall",
	"sig_eps_sq",
	"sig_eps_c_sq",
	"lambda.max",
	"lambda.max_model_size",
	"lambda.min",
	"lambda.min_model_size",
	"lambda_star",
	"lambda_star_model_size",
	"X_ints",
	"y",
	"X_final",
	"y_final",
	"N",
	"T",
	"G",
	"d",
	"p",
	"calc_ses",
	"cv_seed_used"
)

test_that("BETWFE intercept-only block (Block 1) shape: 37 fields, expected ordering, no theta_hat", {
	args <- .make_helper_args()
	# Block 1 input shape: all-zero beta_hat (slopes), no theta_hat,
	# q < 1 -> ret_se = 0.
	args$beta_hat <- rep(0, args$p)
	args$message_text <- "No features selected (or only intercept); all treatment effects estimated to be 0."

	res <- do.call(fetwfe:::.build_selected_out_result, args)

	expect_type(res, "list")
	expect_length(res, 38L)
	expect_identical(names(res), .expected_betwfe_field_order)
	expect_false("theta_hat" %in% names(res))

	# Spot-check key field values.
	expect_identical(res$in_sample_att_hat, 0)
	expect_identical(res$indep_att_hat, 0)
	expect_identical(res$in_sample_att_se, 0) # q < 1 -> ret_se = 0
	expect_identical(res$indep_att_se, 0)
	expect_true(res$calc_ses)

	expect_named(res$catt_hats, c("2", "3"))
	expect_identical(unname(res$catt_hats), c(0, 0))
	expect_identical(unname(res$catt_ses), c(0, 0))

	expect_s3_class(res$catt_df, "data.frame")
	expect_s3_class(res$catt_df, "catt_df")
	expect_named(
		res$catt_df,
		c(
			"cohort",
			"estimate",
			"se",
			"ci_low",
			"ci_high",
			"p_value",
			"selected"
		)
	)
	expect_equal(nrow(res$catt_df), args$G)
	expect_identical(res$catt_df$cohort, args$c_names)
	expect_identical(res$catt_df$selected, rep(FALSE, args$G))
	expect_identical(res$catt_df$p_value, rep(NA_real_, args$G))

	expect_identical(res$beta_hat, rep(0, args$p))
	expect_identical(res$N, args$N)
	expect_identical(res$T, args$T)
	expect_identical(res$G, args$G)
	expect_identical(res$d, args$d)
	expect_identical(res$p, args$p)
})

test_that("BETWFE no-treatment block (Block 2) shape: caller-applied add_ridge scaling propagates", {
	args <- .make_helper_args()
	args$message_text <- "No treatment features selected; all treatment effects estimated to be 0."
	# Block 2 input shape: non-zero beta_hat (some non-treatment
	# coefs selected), no theta_hat. Simulate caller-applied
	# `add_ridge` scaling: beta_hat * (1 + lambda_ridge).
	base_beta <- c(0.5, -0.3, 0, 0, 0, 0.2, 0)
	lambda_ridge <- 0.1
	args$beta_hat <- base_beta * (1 + lambda_ridge)

	res <- do.call(fetwfe:::.build_selected_out_result, args)

	expect_length(res, 38L)
	expect_identical(names(res), .expected_betwfe_field_order)
	expect_false("theta_hat" %in% names(res))

	expect_identical(res$beta_hat, base_beta * 1.1)
	expect_identical(res$in_sample_att_hat, 0)
	expect_identical(unname(res$catt_hats), c(0, 0))
})

test_that("FETWFE intercept-only block (Block 3) shape: 39 fields, theta_hat between catt_df and beta_hat", {
	args <- .make_helper_args()
	args$message_text <- "No features selected (or only intercept); all treatment effects estimated to be 0."
	args$beta_hat <- rep(0, args$p)
	args$theta_hat <- c(1.5, rep(0, args$p)) # intercept + zero slopes
	args$include_theta <- TRUE

	res <- do.call(fetwfe:::.build_selected_out_result, args)

	expect_length(res, 39L)
	expect_identical(names(res), .expected_fetwfe_field_order)

	# theta_hat present, immediately between catt_df (pos 12) and beta_hat (pos 14)
	# (positions shifted by 4 after the 4 new variance-component slots).
	expect_true("theta_hat" %in% names(res))
	expect_identical(which(names(res) == "theta_hat"), 13L)
	expect_identical(which(names(res) == "beta_hat"), 14L)
	expect_identical(which(names(res) == "catt_df"), 12L)

	expect_identical(res$theta_hat, args$theta_hat)
	expect_identical(res$beta_hat, rep(0, args$p))
})

test_that("FETWFE no-treatment block (Block 4) shape: caller passes untransformed beta_hat plus theta_hat", {
	args <- .make_helper_args()
	args$message_text <- "No treatment features selected; all treatment effects estimated to be 0."
	# Block 4 caller computes beta_hat via `untransformCoefImproved()`
	# then optionally scales by (1 + lambda_ridge). Simulate the
	# final value the helper receives.
	args$beta_hat <- c(0.4, -0.2, 0, 0, 0, 0.1, 0) * 1.05 # post-untransform + add_ridge
	args$theta_hat <- c(2.0, 0.3, -0.1, 0, 0, 0, 0.15, 0) # full theta with intercept
	args$include_theta <- TRUE

	res <- do.call(fetwfe:::.build_selected_out_result, args)

	expect_length(res, 39L)
	expect_identical(names(res), .expected_fetwfe_field_order)
	expect_identical(res$theta_hat, args$theta_hat)
	expect_identical(res$beta_hat, args$beta_hat)
})

test_that("q gating: q >= 1 -> ret_se = NA, calc_ses = FALSE; q < 1 -> ret_se = 0, calc_ses = TRUE", {
	args_low <- .make_helper_args()
	args_low$q <- 0.5
	res_low <- do.call(fetwfe:::.build_selected_out_result, args_low)

	expect_identical(res_low$in_sample_att_se, 0)
	expect_identical(res_low$indep_att_se, 0)
	expect_identical(unname(res_low$catt_ses), c(0, 0))
	expect_identical(res_low$catt_df$se, c(0, 0))
	expect_true(res_low$calc_ses)

	args_high <- .make_helper_args()
	args_high$q <- 1
	res_high <- do.call(fetwfe:::.build_selected_out_result, args_high)

	expect_identical(res_high$in_sample_att_se, NA)
	expect_identical(res_high$indep_att_se, NA)
	expect_true(all(is.na(res_high$catt_ses)))
	expect_true(all(is.na(res_high$catt_df$se)))
	expect_false(res_high$calc_ses)
})

test_that("gls gating (#304): q < 1 but gls = FALSE -> ret_se = NA, calc_ses = FALSE", {
	# A gls = FALSE degenerate fit has no oracle SE, so its degenerate SE is NA and
	# calc_ses is FALSE -- mirroring the normal path's `(q < 1) && gls`. calc_ses,
	# the SEs, and the variance components must all agree (else the C1 SE-consistency
	# validator trips at fit construction).
	args <- .make_helper_args()
	args$q <- 0.5
	args$gls <- FALSE
	res <- do.call(fetwfe:::.build_selected_out_result, args)
	expect_false(res$calc_ses)
	expect_identical(res$in_sample_att_se, NA)
	expect_identical(res$indep_att_se, NA)
	expect_true(all(is.na(res$catt_ses)))
	expect_true(all(is.na(res$catt_df$se)))
	# the default (gls absent) keeps the old q < 1 -> calc_ses = TRUE behavior.
	args_default <- .make_helper_args()
	args_default$q <- 0.5
	expect_true(
		do.call(fetwfe:::.build_selected_out_result, args_default)$calc_ses
	)
})

test_that("verbose = TRUE emits the supplied message_text; verbose = FALSE is silent", {
	args <- .make_helper_args()
	args$verbose <- TRUE
	args$message_text <- "No features selected (or only intercept); all treatment effects estimated to be 0."
	expect_message(
		do.call(fetwfe:::.build_selected_out_result, args),
		"No features selected"
	)

	args$verbose <- FALSE
	expect_no_message(do.call(fetwfe:::.build_selected_out_result, args))
})
