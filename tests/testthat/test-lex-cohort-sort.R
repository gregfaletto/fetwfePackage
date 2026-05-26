# Tests for the lexicographic-vs-numeric cohort sort fix.
#
# Background: prior to v1.9.3, .truncate_catt() (R/class_helpers.R:40)
# and .tidy_estimator_output() (R/broom_methods.R:57) both used
# order(catt$Cohort) on a character column. R's default order() on
# character is lexicographic, so for cohorts ≥ 10 the output came back
# as "10","11","12","2","3" rather than "2","3","10","11","12". For
# .truncate_catt(), the worst consequence is that print()/summary() with
# max_cohorts truncation dropped the *earliest* cohorts as if they were
# "later" alphabetically. See GitHub issue #53.

# Build a 11-row catt_df with mixed-digit-count cohort labels to exercise
# the bug. Without the fix, both order() outputs would be lex-sorted.
.make_catt_df_11 <- function() {
	out <- data.frame(
		cohort = as.character(2:12),
		estimate = seq(0.01, 0.11, length.out = 11),
		se = rep(0.005, 11),
		ci_low = seq(0.0, 0.10, length.out = 11),
		ci_high = seq(0.02, 0.12, length.out = 11),
		p_value = rep(0.05, 11),
		selected = rep(TRUE, 11),
		stringsAsFactors = FALSE
	)
	class(out) <- c("catt_df", "data.frame")
	out
}

test_that(".truncate_catt sorts cohorts numerically and drops the right cohort on truncation", {
	df <- .make_catt_df_11()

	# Truncate to 10 — should keep cohorts 2..11 (the earliest 10) and
	# drop "12". Under the pre-fix lex sort, the kept rows would be
	# "10","11","12","2","3","4","5","6","7","8" and the dropped cohort
	# would be "9" (wrong).
	res <- fetwfe:::.truncate_catt(
		df,
		max_cohorts = 10,
		order_by = "cohort"
	)
	expect_equal(res$cohort, as.character(2:11))
	expect_true(isTRUE(attr(res, "truncated")))
	expect_equal(attr(res, "n_discarded"), 1L)

	# Confirm the dropped cohort is "12" (the largest), not "9" (which
	# would be the pre-fix result).
	expect_false("12" %in% res$cohort)
	expect_true("9" %in% res$cohort)
})

test_that(".truncate_catt full-output sort order is numeric (no truncation)", {
	df <- .make_catt_df_11()
	res <- fetwfe:::.truncate_catt(
		df,
		max_cohorts = 20,
		order_by = "cohort"
	)
	expect_equal(res$cohort, as.character(2:12))
	expect_false(isTRUE(attr(res, "truncated")))
})

test_that("tidy.<class> sorts cohorts numerically when labels include >= 10", {
	skip_if_not_installed("broom")

	# Build a fully-validator-compliant fetwfe-classed object with a
	# catt_df containing cohorts >= 10. Most slots are NA/0/empty
	# placeholders; only catt_df, R, alpha, and the SE-consistency
	# slots need to be coherent. The constructor validator from #85
	# (called by .check_for_tidy from #86) requires the full slot
	# inventory; pre-#85 versions of this test got away with a
	# minimal mock.
	R_test <- 11L
	T_test <- 13L
	d_test <- 2L
	p_test <- 50L
	num_treats <- T_test * R_test - R_test * (R_test + 1L) / 2L # 88
	obj <- list(
		att_hat = 0.05,
		att_se = 0.01,
		att_p_value = 0.001,
		att_selected = TRUE,
		catt_hats = setNames(rep(0.01, R_test), as.character(2:12)),
		catt_ses = setNames(rep(0.005, R_test), as.character(2:12)),
		cohort_probs = rep(1 / (R_test + 5L), R_test), # < 1/(R+1) each so sum < 1
		cohort_probs_overall = rep(1 / (R_test + 5L), R_test),
		catt_df = .make_catt_df_11(),
		beta_hat = rep(0, p_test),
		treat_inds = seq_len(num_treats),
		treat_int_inds = (num_treats + 1L):p_test,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		lambda.max = 1,
		lambda.max_model_size = 1L,
		lambda.min = 0.01,
		lambda.min_model_size = p_test + 1L,
		lambda_star = 0.1,
		lambda_star_model_size = 5L,
		N = 100L,
		T = T_test,
		R = R_test,
		d = d_test,
		p = p_test,
		alpha = 0.05,
		indep_counts_used = FALSE,
		se_type = "default",
		y_mean = 0,
		response_col_name = "y",
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		internal = list(
			X_ints = matrix(0, 100L * T_test, p_test),
			y = rep(0, 100L * T_test),
			X_final = matrix(0, 100L * T_test, p_test),
			y_final = rep(0, 100L * T_test),
			theta_hat = rep(0, p_test + 1L),
			calc_ses = TRUE
		)
	)
	class(obj) <- "fetwfe"

	td <- broom::tidy(obj)
	cohort_rows <- td[td$term != "ATT", ]
	cohort_labels <- sub("^Cohort ", "", cohort_rows$term)

	# Numeric sort: "2","3",...,"12".
	expect_equal(cohort_labels, as.character(2:12))

	# Negative check: confirm this is NOT the lex-sort outcome
	# ("10","11","12","2","3",...).
	expect_false(identical(cohort_labels, sort(as.character(2:12))))
})
