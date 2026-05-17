# Tests for the shared S3-class helpers in R/class_helpers.R.
#
# These exercise the cross-class invariant that the per-class options
# `fetwfe.max_cohorts`, `etwfe.max_cohorts`, `betwfe.max_cohorts` remain
# independent after `.truncate_catt` was deduplicated into a single
# definition. The helper now takes `max_cohorts` as a required positional
# argument; the per-class option lives on each `print.<class>()` method's
# own `max_cohorts =` argument, and the three should not collide.

# Helper: build a 5-row catt_df with the columns print.<class> reads. The
# Cohort values are descending so `order_by = "cohort"` (the default for
# print methods) reorders them ascending; that lets the truncation behavior
# be observed deterministically.
.make_catt_df_5 <- function() {
	data.frame(
		Cohort = paste0("c", 5:1),
		`Estimated TE` = seq(0.1, 0.5, length.out = 5),
		SE = rep(0.05, 5),
		ConfIntLow = seq(0.0, 0.4, length.out = 5),
		ConfIntHigh = seq(0.2, 0.6, length.out = 5),
		P_value = rep(0.01, 5),
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}

# Helper: build a minimal classed object carrying only the fields each
# `print.<class>(show_internal = FALSE)` body touches. The exact field
# inventory was confirmed by reading each print method end-to-end during
# the plan-review pass.
.make_minimal_obj <- function(klass) {
	obj <- list(
		alpha = 0.05,
		att_hat = 0.3,
		att_se = 0.05,
		att_p_value = 0.01,
		# se_type must mirror the live `match.arg(c("default", "cluster"))`
		# choice in the estimator entry points; "standard" was a stale
		# mock value (issue #55).
		se_type = "default",
		catt_df = .make_catt_df_5(),
		N = 100,
		T = 6,
		R = 3,
		d = 2,
		p = 30
	)
	if (klass %in% c("fetwfe", "betwfe")) {
		obj$att_selected <- TRUE
		obj$lambda_star_model_size <- 5L
		obj$lambda_star <- 0.1
	}
	class(obj) <- klass
	obj
}

test_that("per-class max_cohorts options are independent across print methods", {
	prev <- options(
		fetwfe.max_cohorts = 1,
		etwfe.max_cohorts = 2,
		betwfe.max_cohorts = 3
	)
	on.exit(options(prev), add = TRUE)

	out_fetwfe <- paste(
		capture.output(print(.make_minimal_obj("fetwfe"))),
		collapse = "\n"
	)
	out_etwfe <- paste(
		capture.output(print(.make_minimal_obj("etwfe"))),
		collapse = "\n"
	)
	out_betwfe <- paste(
		capture.output(print(.make_minimal_obj("betwfe"))),
		collapse = "\n"
	)

	# Each print method should truncate to its own option's value: 5 - K rows
	# are discarded, so the suffix reads "... and (5 - K) more cohorts."
	expect_match(out_fetwfe, "and 4 more cohorts", fixed = TRUE)
	expect_match(out_etwfe, "and 3 more cohorts", fixed = TRUE)
	expect_match(out_betwfe, "and 2 more cohorts", fixed = TRUE)
})

test_that(".truncate_catt requires max_cohorts (no default)", {
	# Regression guard: after deduplication, .truncate_catt has no default
	# on `max_cohorts`. Calling without it should error.
	expect_error(
		fetwfe:::.truncate_catt(.make_catt_df_5()),
		'argument "max_cohorts" is missing'
	)
})
