library(testthat)
library(fetwfe)

# Issue #439 (#366 part A2): coverage for the five rendering pieces shared by
# `.print_estimator_output()` and `.print_summary_estimator_output()` --
# `.se_qualifier()`, `.band_label()`, `.att_wald_ci()`, `.cat_preview_block()`
# and `.cat_model_details()`.
#
# WHY THIS FILE EXISTS. The consolidation is behavior-preserving, so none of
# these assertions goes red before the source change; that is what a guardrail
# is. Their value is measured instead by a mutation battery, and each group
# below names the mutation it defends against. Before #439 the suite could not
# distinguish a correct helper from a stubbed one: measured during plan review,
# `.band_label <- function(ci_type) "simultaneous"`,
# `.se_qualifier <- function(se_type) ""`, and swapping the G and d rows of the
# Model Details block EACH left the entire suite green.
#
# ---------------------------------------------------------------------------
# `fixed = TRUE` IS MANDATORY ON EVERY RENDERED-LITERAL ASSERTION HERE.
#
# `expect_match()` and `grepl()` default to `fixed = FALSE`, and three of the
# literals this file pins contain regex metacharacters. Measured:
#
#   grepl("[pointwise 95% CI]", "...(CATT) [simultaneous 95% CI]:")   TRUE (!)
#   grepl("SE (cluster-robust) =", "...(SE (cluster-robust) = 0,...") FALSE (!)
#   grepl("  ... + 7 more cohorts.", "  ... + 7 more cohorts.")       FALSE (!)
#
# `[pointwise 95% CI]` is a CHARACTER CLASS, so as a regex it matches almost any
# non-empty string -- including output from the `.band_label -> "simultaneous"`
# mutation it is the only guard against. `(cluster-robust)` is a CAPTURE GROUP,
# so those presence assertions fail against correct output while their paired
# absence assertions pass by accident. The `+` in the footer is a QUANTIFIER, so
# that assertion cannot match its own expected text.
#
# Absence is asserted as `expect_false(any(grepl(..., fixed = TRUE)))` rather
# than `expect_no_match`, so the `fixed` argument is impossible to drop by
# inattention. `tests/testthat/test-class-helpers.R` uses the same convention.
# ---------------------------------------------------------------------------
#
# Panel-data fixture helpers (generate_panel_data, ...) live in
# tests/testthat/helper-panel-fixture.R and are sourced by testthat first (#91).

.render_print <- function(x, ...) {
	paste(capture.output(print(x, ...)), collapse = "\n")
}

.render_summary <- function(x, ...) {
	paste(capture.output(print(summary(x, ...))), collapse = "\n")
}

.expect_absent <- function(blob, literal) {
	expect_false(any(grepl(literal, blob, fixed = TRUE)))
}

# ------------------------------------------------------------------------------
# Group 1: direct helper contracts.
#
# Defends against: any helper stubbed to a constant. These are the only
# assertions that see the helpers in isolation; everything below sees them
# through a call site, which is the distinction the mutation battery turns on.
# ------------------------------------------------------------------------------

test_that(".se_qualifier() returns the parenthetical, not a finished label", {
	expect_identical(fetwfe:::.se_qualifier("cluster"), " (cluster-robust)")
	expect_identical(fetwfe:::.se_qualifier("conservative"), " (conservative)")
	expect_identical(fetwfe:::.se_qualifier("default"), "")
	# Anything unrecognized falls through to the unqualified form rather than
	# erroring: `se_type` is validated at fit time (C8), so this is a rendering
	# helper's tolerance, not an input gate.
	expect_identical(fetwfe:::.se_qualifier("nonsense"), "")
	expect_identical(fetwfe:::.se_qualifier(NULL), "")

	# The leading space is load-bearing: the print caller composes
	# sprintf("  Std. Error%s: %.4f\n", ...), so dropping it would render
	# "Std. Error(cluster-robust):".
	expect_identical(
		sprintf("  Std. Error%s: %.4f\n", fetwfe:::.se_qualifier("cluster"), 1),
		"  Std. Error (cluster-robust): 1.0000\n"
	)
	# ...and the summary caller composes paste0("SE", .), with no space.
	expect_identical(
		paste0("SE", fetwfe:::.se_qualifier("cluster")),
		"SE (cluster-robust)"
	)
})

test_that(".band_label() takes the scalar and treats a missing slot as pointwise", {
	expect_identical(fetwfe:::.band_label("simultaneous"), "simultaneous")
	expect_identical(fetwfe:::.band_label("pointwise"), "pointwise")
	# NULL is the pre-1.16.0 no-slot case. Those objects' stored bounds really
	# are pointwise, so this is a compatibility behavior, not a fallback
	# default -- see the helper's roxygen and R/event_study.R's
	# `.resolve_event_study_ci_type()`, which encodes the same rule.
	expect_identical(fetwfe:::.band_label(NULL), "pointwise")
})

test_that(".att_wald_ci() computes the Wald interval and reads its alpha", {
	# Hand-computed rather than recomputed from the helper's own output: a
	# round-trip against `stats::qnorm` would hold by construction.
	expect_equal(
		fetwfe:::.att_wald_ci(2, 0.5, 0.05),
		c(2 - 1.959963984540054 * 0.5, 2 + 1.959963984540054 * 0.5),
		tolerance = 1e-12
	)
	# A DIFFERENT alpha must produce a different, wider interval. Without this
	# the helper could ignore `alpha` entirely.
	wide <- fetwfe:::.att_wald_ci(2, 0.5, 0.01)
	narrow <- fetwfe:::.att_wald_ci(2, 0.5, 0.05)
	expect_gt(wide[2] - wide[1], narrow[2] - narrow[1])
	expect_equal(
		narrow[2] - narrow[1],
		2 * stats::qnorm(1 - 0.05 / 2) * 0.5,
		tolerance = 1e-12
	)

	# Lower bound first.
	expect_lt(narrow[1], narrow[2])
	# NA `se` (a q >= 1 fit, where calc_ses is FALSE) propagates rather than
	# erroring -- the render paths print "[NA, NA]".
	expect_true(all(is.na(fetwfe:::.att_wald_ci(2, NA_real_, 0.05))))

	# The unnamed-return property, asserted with NAMED inputs. With unnamed
	# inputs this passes even with `unname()` deleted, and would be a vacuous
	# fixture: `unname()` is load-bearing at exactly one call site, the summary
	# one, where the inputs are `x$att["estimate"]` and `x$att["se"]`.
	expect_null(names(fetwfe:::.att_wald_ci(
		c(estimate = 0.5),
		c(se = 0.1),
		0.05
	)))
})

test_that(".cat_model_details() rejects a mis-named field instead of dropping a row", {
	ok <- list(N = 40, T = 6, G = 3, d = 1, p = 41)
	expect_silent(invisible(capture.output(
		fetwfe:::.cat_model_details(ok, show_lambda = FALSE)
	)))

	# `sprintf("%d", NULL)` returns character(0) and `cat(character(0))` emits
	# nothing, so without the name check this would print a Model Details block
	# silently missing the Units row.
	bad <- list(Ntotal = 40, T = 6, G = 3, d = 1, p = 41)
	expect_error(
		fetwfe:::.cat_model_details(bad, show_lambda = FALSE),
		"missing field\\(s\\): N"
	)
	# The two optional fields are required exactly when they are read.
	expect_error(
		fetwfe:::.cat_model_details(ok, show_lambda = TRUE),
		"missing field\\(s\\): model_size, lambda_star"
	)
	# Extra names are ignored -- this is what lets the summary caller pass
	# `model_info` (which also carries R / sig_eps_sq / sig_eps_c_sq) unadapted.
	extra <- c(ok, list(R = 3, sig_eps_sq = 1, sig_eps_c_sq = 0.5))
	out <- capture.output(fetwfe:::.cat_model_details(
		extra,
		show_lambda = FALSE
	))
	expect_true(any(grepl("  Units (N)           : 40", out, fixed = TRUE)))
})

test_that(".cat_preview_block() rejects a truncated frame with no discard count", {
	df <- data.frame(cohort = "2", estimate = 1)
	attr(df, "truncated") <- TRUE
	# `n_discarded` deliberately unset: without the guard, `sprintf(fmt, NULL)`
	# returns character(0) and the footer vanishes silently -- in the helper
	# whose entire purpose is that footer.
	expect_error(
		fetwfe:::.cat_preview_block(df, "H:\n", "  ... and %d more cohorts.\n"),
		"`n_discarded` is missing or not a scalar"
	)

	attr(df, "n_discarded") <- 4L
	out <- capture.output(
		fetwfe:::.cat_preview_block(df, "H:\n", "  ... and %d more cohorts.\n")
	)
	expect_true(any(grepl("  ... and 4 more cohorts.", out, fixed = TRUE)))
})

# ------------------------------------------------------------------------------
# Group 2: the label branches, end to end, on BOTH render paths.
#
# Defends against: `.band_label()` stubbed to "simultaneous" and
# `.se_qualifier()` stubbed to "" -- each of which left the whole suite green
# before this file existed, because all ten snapshot goldens are the "default"
# se_type and the "simultaneous" ci_type. Also defends against a branch swap
# inside `.se_qualifier()`, via the paired absence assertions.
# ------------------------------------------------------------------------------

.pdata_439 <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

.fit_439 <- function(...) {
	fetwfe(
		pdata = .pdata_439,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		verbose = FALSE,
		lambda_selection = "bic",
		...
	)
}

test_that("se_type = 'cluster' and ci_type = 'pointwise' are rendered on both paths", {
	fit <- .fit_439(se_type = "cluster", ci_type = "pointwise")

	out_print <- .render_print(fit)
	expect_match(out_print, "  Std. Error (cluster-robust): ", fixed = TRUE)
	expect_match(out_print, "[pointwise 95% CI]", fixed = TRUE)
	.expect_absent(out_print, "[simultaneous 95% CI]")
	.expect_absent(out_print, " (conservative)")

	out_summary <- .render_summary(fit)
	expect_match(out_summary, "(SE (cluster-robust) = ", fixed = TRUE)
	expect_match(out_summary, "[pointwise 95% CI]", fixed = TRUE)
	.expect_absent(out_summary, "[simultaneous 95% CI]")
	.expect_absent(out_summary, " (conservative)")
})

test_that("se_type = 'conservative' is rendered on both paths", {
	fit <- .fit_439(se_type = "conservative")

	out_print <- .render_print(fit)
	expect_match(out_print, "  Std. Error (conservative): ", fixed = TRUE)
	.expect_absent(out_print, " (cluster-robust)")
	# The default ci_type travels with this fit, so it also pins the other
	# branch of `.band_label()`.
	expect_match(out_print, "[simultaneous 95% CI]", fixed = TRUE)
	.expect_absent(out_print, "[pointwise 95% CI]")

	out_summary <- .render_summary(fit)
	expect_match(out_summary, "(SE (conservative) = ", fixed = TRUE)
	.expect_absent(out_summary, " (cluster-robust)")
	expect_match(out_summary, "[simultaneous 95% CI]", fixed = TRUE)
	.expect_absent(out_summary, "[pointwise 95% CI]")
})

test_that("the default se_type renders an unqualified Std. Error on both paths", {
	fit <- .fit_439()

	out_print <- .render_print(fit)
	expect_match(out_print, "  Std. Error: ", fixed = TRUE)
	.expect_absent(out_print, " (cluster-robust)")
	.expect_absent(out_print, " (conservative)")

	out_summary <- .render_summary(fit)
	expect_match(out_summary, "(SE = ", fixed = TRUE)
	.expect_absent(out_summary, " (cluster-robust)")
	.expect_absent(out_summary, " (conservative)")
})

test_that("a non-default alpha reaches BOTH rendered CI labels", {
	# Without this, `.att_wald_ci()` ignoring its `alpha` and hard-coding 0.05
	# would be byte-identical at both rendering call sites: measured, nothing
	# else in the repo renders a print or summary path at a non-default alpha.
	# The broom site is covered by test-broom-methods.R's conf.level test.
	fit <- .fit_439(alpha = 0.10)

	out_print <- .render_print(fit)
	expect_match(out_print, "90% CI:", fixed = TRUE)
	expect_match(out_print, "[simultaneous 90% CI]", fixed = TRUE)
	.expect_absent(out_print, "95% CI")

	out_summary <- .render_summary(fit)
	expect_match(out_summary, "90% CI = ", fixed = TRUE)
	expect_match(out_summary, "[simultaneous 90% CI]", fixed = TRUE)
	.expect_absent(out_summary, "95% CI")
})

# ------------------------------------------------------------------------------
# Group 3: Model Details on a fixture with all five dimensions distinct.
#
# Defends against: a G/d (or N/T, or any other) transposition in
# `.cat_model_details()`. Every snapshot fixture predating #439 had G == d == 2,
# so a G/d swap rendered identically and the whole suite stayed green.
# ------------------------------------------------------------------------------

.pdata_distinct_439 <- generate_panel_data(N = 40, T = 6, R = 3, seed = 123)

.fit_distinct_439 <- fetwfe(
	pdata = .pdata_distinct_439,
	time_var = "time",
	unit_var = "unit",
	treatment = "treatment",
	response = "y",
	covs = "cov1",
	verbose = FALSE,
	lambda_selection = "bic"
)

test_that("the distinct-dimension fixture really has five distinct dimensions", {
	# Asserted first and separately so that if a future change to
	# `generate_panel_data()` alters the realized cohort count, this fails as a
	# FIXTURE problem rather than silently degrading into a G == d fixture that
	# proves nothing.
	#
	# `expect_equal`, not `expect_identical`: `x$G` and `x$p` are doubles while
	# `x$N`, `x$T` and `x$d` are integers, so the identical form would fail on
	# type and read as a broken assertion.
	dims <- c(
		.fit_distinct_439$N,
		.fit_distinct_439$T,
		.fit_distinct_439$G,
		.fit_distinct_439$d,
		.fit_distinct_439$p
	)
	expect_equal(dims, c(40, 6, 3, 1, 41))
	expect_length(unique(dims), 5L)
})

test_that("Model Details maps each dimension to its own label on both paths", {
	expected <- c(
		"  Units (N)           : 40",
		"  Time periods (T)    : 6",
		"  Treated cohorts (G) : 3",
		"  Covariates (d)      : 1",
		"  Features (p)        : 41"
	)
	for (blob in list(
		.render_print(.fit_distinct_439),
		.render_summary(.fit_distinct_439)
	)) {
		for (line in expected) {
			expect_match(blob, line, fixed = TRUE)
		}
		# A G/d transposition would render these instead.
		.expect_absent(blob, "  Treated cohorts (G) : 1")
		.expect_absent(blob, "  Covariates (d)      : 3")
		# ...and an N/T transposition these.
		.expect_absent(blob, "  Units (N)           : 6")
		.expect_absent(blob, "  Time periods (T)    : 40")
	}
})

# ------------------------------------------------------------------------------
# Group 4: all FOUR truncation footers.
#
# Defends against: `.cat_preview_block()` ignoring its `more_fmt` argument. The
# four literals vary on two axes at once -- "and" vs "+" (print vs summary) and
# "cohorts" vs "event times" (CATT vs event study) -- so a helper that
# hard-coded any one of them renders plausible output at the wrong site.
#
# Before #439 only the two PRINT-side footers were pinned anywhere in the repo
# (test-class-helpers.R and test-event-study-in-display-138.R). The two
# summary-side literals were pinned nowhere and reached by nothing.
# ------------------------------------------------------------------------------

test_that("the print path emits both of its truncation footers", {
	out_cohorts <- .render_print(.fit_distinct_439, max_cohorts = 1)
	expect_match(out_cohorts, "  ... and 2 more cohorts.", fixed = TRUE)
	.expect_absent(out_cohorts, "  ... + 2 more cohorts.")
	.expect_absent(out_cohorts, "more event times")

	out_es <- .render_print(.fit_distinct_439, max_event_times = 1)
	expect_match(out_es, "  ... and 4 more event times.", fixed = TRUE)
	.expect_absent(out_es, "  ... + 4 more event times.")
	.expect_absent(out_es, "more cohorts")
})

test_that("the summary path emits both of its truncation footers", {
	# `summary.<class>()` exposes neither `max_cohorts` nor `max_event_times`,
	# and `.summary_estimator_output()` hard-codes 20 for both, so the only
	# public route to these two footers is a fit with more than 20 cohorts or
	# event times -- measured at ~10s. Setting the attributes on a real summary
	# object reaches them at no cost and still runs the REAL
	# `print.summary.fetwfe()` -> `.print_summary_estimator_output()` ->
	# `.cat_preview_block()` with the real `more_fmt` argument, which is the
	# thing under test.
	#
	# NOTE: these two objects are internally INCONSISTENT on purpose. The real
	# `.truncate_catt()` sets `truncated` and then subsets the frame; here the
	# frame is unsubsetted and the count is invented, so the table above the
	# footer disagrees with the footer. They can pin ONLY the footer literal --
	# do not extend them into any assertion relating the table to the footer.
	s_catt <- summary(.fit_distinct_439)
	attr(s_catt$catt, "truncated") <- TRUE
	attr(s_catt$catt, "n_discarded") <- 7L
	out_catt <- paste(capture.output(print(s_catt)), collapse = "\n")
	expect_match(out_catt, "  ... + 7 more cohorts.", fixed = TRUE)
	.expect_absent(out_catt, "  ... and 7 more cohorts.")
	.expect_absent(out_catt, "7 more event times")

	s_es <- summary(.fit_distinct_439)
	attr(s_es$event_study, "truncated") <- TRUE
	attr(s_es$event_study, "n_discarded") <- 4L
	out_es <- paste(capture.output(print(s_es)), collapse = "\n")
	expect_match(out_es, "  ... + 4 more event times.", fixed = TRUE)
	.expect_absent(out_es, "  ... and 4 more event times.")
	.expect_absent(out_es, "4 more cohorts")
})
