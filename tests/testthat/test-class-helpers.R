# Tests for the shared S3-class helpers in R/class_helpers.R.
#
# These exercise the cross-class invariant that the per-class options
# `fetwfe.max_cohorts`, `etwfe.max_cohorts`, `betwfe.max_cohorts` remain
# independent after `.truncate_catt` was deduplicated into a single
# definition. The helper now takes `max_cohorts` as a required positional
# argument; the per-class option lives on each `print.<class>()` method's
# own `max_cohorts =` argument, and the three should not collide.

# Helper: build a 5-row catt_df with the columns print.<class> reads. The
# cohort values are descending so `order_by = "cohort"` (the default for
# print methods) reorders them ascending; that lets the truncation behavior
# be observed deterministically.
.make_catt_df_5 <- function() {
	out <- data.frame(
		cohort = paste0("c", 5:1),
		estimate = seq(0.1, 0.5, length.out = 5),
		se = rep(0.05, 5),
		ci_low = seq(0.0, 0.4, length.out = 5),
		ci_high = seq(0.2, 0.6, length.out = 5),
		p_value = rep(0.01, 5),
		stringsAsFactors = FALSE
	)
	class(out) <- c("catt_df", "data.frame")
	out
}

# Helper: build a printable classed object on a small simulated fixture,
# then overwrite `catt_df` with the 5-row stub above. Prior to #174 we
# used a sparse hand-rolled list with only the slots `print.<class>()`
# read; #174 removed the silent-swallow `tryCatch(eventStudy(x), ...)`
# in `print` so it now also exercises the full slot inventory via
# `eventStudy()`'s precondition (which calls `.validate_<class>()`).
# Sourcing the fixture once per call keeps the test focused on
# truncation behavior while satisfying the post-#174 contract.
# The fixture pins R = 5 because the validator's C4 contract requires
# `nrow(catt_df) == R`; the 5-row stub above is the load-bearing piece
# for the truncation assertions, so the fixture's R has to match.
.make_minimal_obj <- function(klass) {
	set.seed(2026)
	coefs <- genCoefs(R = 5, T = 7, d = 0, density = 0.5, eff_size = 1)
	sim <- simulateData(coefs, N = 80, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	obj <- switch(
		klass,
		fetwfe = fetwfeWithSimulatedData(sim, q = 0.5),
		etwfe = etwfeWithSimulatedData(sim),
		betwfe = betwfeWithSimulatedData(sim, q = 0.5),
		stop("Unknown klass: ", klass)
	)
	# Pre-#174 the mock pinned `R = 3 / N = 100 / lambda_star = 0.1`
	# etc.; only `catt_df` is load-bearing for truncation behavior, so
	# we keep the 5-row stub there and let the real fit supply the rest.
	obj$catt_df <- .make_catt_df_5()
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
