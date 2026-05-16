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
	data.frame(
		Cohort = as.character(2:12),
		`Estimated TE` = seq(0.01, 0.11, length.out = 11),
		SE = rep(0.005, 11),
		ConfIntLow = seq(0.0, 0.10, length.out = 11),
		ConfIntHigh = seq(0.02, 0.12, length.out = 11),
		P_value = rep(0.05, 11),
		selected = rep(TRUE, 11),
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
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
	expect_equal(res$Cohort, as.character(2:11))
	expect_true(isTRUE(attr(res, "truncated")))
	expect_equal(attr(res, "n_discarded"), 1L)

	# Confirm the dropped cohort is "12" (the largest), not "9" (which
	# would be the pre-fix result).
	expect_false("12" %in% res$Cohort)
	expect_true("9" %in% res$Cohort)
})

test_that(".truncate_catt full-output sort order is numeric (no truncation)", {
	df <- .make_catt_df_11()
	res <- fetwfe:::.truncate_catt(
		df,
		max_cohorts = 20,
		order_by = "cohort"
	)
	expect_equal(res$Cohort, as.character(2:12))
	expect_false(isTRUE(attr(res, "truncated")))
})

test_that("tidy.<class> sorts cohorts numerically when labels include >= 10", {
	skip_if_not_installed("broom")

	# Build a minimal fetwfe-classed object with a catt_df containing
	# cohorts >= 10. This bypasses the full estimator pipeline; we only
	# need .tidy_estimator_output to be exercised.
	obj <- list(
		att_hat = 0.05,
		att_se = 0.01,
		att_p_value = 0.001,
		att_selected = TRUE,
		alpha = 0.05,
		catt_df = .make_catt_df_11(),
		N = 100,
		T = 13,
		R = 11,
		d = 2,
		p = 50
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
