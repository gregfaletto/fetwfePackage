# Helpful-error S3 layer on the `catt_df` class (issue #136).
#
# As of fetwfe 1.11.0 the columns of `result$catt_df` were renamed from
# Title-Case to snake_case to match `eventStudy(result)`'s convention.
# To make migration loud rather than silent, every `catt_df` carries
# class c("catt_df", "data.frame") and the three S3 methods on
# `[[` / `$` / `[` intercept the six old column names and stop()
# with a migration message pointing to the new name.

library(testthat)
library(fetwfe)

# -------- Fixture: build a `catt_df`-classed data frame ------------------
.make_catt_df <- function() {
	out <- data.frame(
		cohort = c("2", "3", "4"),
		estimate = c(0.1, 0.2, 0.3),
		se = c(0.05, 0.06, 0.07),
		ci_low = c(0.0, 0.08, 0.16),
		ci_high = c(0.2, 0.32, 0.44),
		p_value = c(0.04, 0.01, 0.005),
		selected = c(TRUE, TRUE, TRUE),
		stringsAsFactors = FALSE
	)
	class(out) <- c("catt_df", "data.frame")
	out
}

# The 6 names that should fire the helpful-error layer.
.old_names <- c(
	"Cohort",
	"Estimated TE",
	"SE",
	"ConfIntLow",
	"ConfIntHigh",
	"P_value"
)

# Their new names; the migration message must surface the new name.
.new_names <- c(
	"cohort",
	"estimate",
	"se",
	"ci_low",
	"ci_high",
	"p_value"
)

# ------------------------------------------------------------------------
# 1) Regression guard: the rename-map covers exactly the historical names
# ------------------------------------------------------------------------
test_that(".catt_df_rename_map maps the historical names to the new names", {
	# The internal map is the single source of truth for the rename. If
	# this assertion ever fails, the rename map drifted from the historical
	# column set; downstream S3 dispatch and the helpful-error messages
	# (which interpolate the new name out of the map) will be inconsistent
	# with NEWS.md.
	expect_identical(names(fetwfe:::.catt_df_rename_map), .old_names)
	expect_identical(unname(fetwfe:::.catt_df_rename_map), .new_names)
})

# ------------------------------------------------------------------------
# 2) The class hierarchy is correct
# ------------------------------------------------------------------------
test_that("catt_df objects are data.frames AND carry the catt_df class", {
	df <- .make_catt_df()
	expect_true(is.data.frame(df))
	expect_s3_class(df, "catt_df")
	expect_s3_class(df, "data.frame")
})

# ------------------------------------------------------------------------
# 3) [[ on each old name fires the helpful error
# ------------------------------------------------------------------------
test_that("[[<old name>] fires helpful error on every renamed column", {
	df <- .make_catt_df()
	for (i in seq_along(.old_names)) {
		old <- .old_names[[i]]
		new <- .new_names[[i]]
		expect_error(
			df[[old]],
			sprintf("`%s` was renamed to `%s` in fetwfe 1.11.0", old, new),
			fixed = TRUE
		)
	}
})

# ------------------------------------------------------------------------
# 4) $ on each old name fires the helpful error
# ------------------------------------------------------------------------
test_that("$<old name> fires helpful error on every renamed column", {
	df <- .make_catt_df()
	for (i in seq_along(.old_names)) {
		old <- .old_names[[i]]
		new <- .new_names[[i]]
		expect_error(
			df[[old]],
			sprintf("`%s` was renamed to `%s` in fetwfe 1.11.0", old, new),
			fixed = TRUE
		)
		# Direct $ access too (needs eval(call) since some old names contain
		# spaces; `df$"Estimated TE"` is valid syntax).
		call <- substitute(df$X, list(X = as.name(old)))
		expect_error(
			eval(call),
			sprintf("`%s` was renamed to `%s` in fetwfe 1.11.0", old, new),
			fixed = TRUE
		)
	}
})

# ------------------------------------------------------------------------
# 5) [ on each old name fires the helpful error in the column-selector
#    position (both df[old] and df[, old] forms)
# ------------------------------------------------------------------------
test_that("[<old name>] (one-index column selection) fires helpful error", {
	df <- .make_catt_df()
	for (i in seq_along(.old_names)) {
		old <- .old_names[[i]]
		new <- .new_names[[i]]
		expect_error(
			df[old],
			sprintf("`%s` was renamed to `%s` in fetwfe 1.11.0", old, new),
			fixed = TRUE
		)
	}
})

test_that("[, <old name>] (two-index column selection) fires helpful error", {
	df <- .make_catt_df()
	for (i in seq_along(.old_names)) {
		old <- .old_names[[i]]
		new <- .new_names[[i]]
		expect_error(
			df[, old],
			sprintf("`%s` was renamed to `%s` in fetwfe 1.11.0", old, new),
			fixed = TRUE
		)
	}
})

# ------------------------------------------------------------------------
# 6) [[ / $ / [ on the new column names work as expected
# ------------------------------------------------------------------------
test_that("[[<new name>] returns the column unchanged", {
	df <- .make_catt_df()
	expect_equal(df[["estimate"]], c(0.1, 0.2, 0.3))
	expect_equal(df[["se"]], c(0.05, 0.06, 0.07))
	expect_equal(df[["cohort"]], c("2", "3", "4"))
	expect_equal(df[["p_value"]], c(0.04, 0.01, 0.005))
})

test_that("$<new name> returns the column unchanged", {
	df <- .make_catt_df()
	expect_equal(df$estimate, c(0.1, 0.2, 0.3))
	expect_equal(df$se, c(0.05, 0.06, 0.07))
	expect_equal(df$cohort, c("2", "3", "4"))
	expect_equal(df$p_value, c(0.04, 0.01, 0.005))
})

test_that("[<new name>] returns the column-frame unchanged", {
	df <- .make_catt_df()
	out <- df["estimate"]
	expect_s3_class(out, "data.frame")
	expect_equal(names(out), "estimate")
	expect_equal(out[[1]], c(0.1, 0.2, 0.3))
})

test_that("[, <new name>] returns the column vector unchanged", {
	df <- .make_catt_df()
	expect_equal(df[, "estimate"], c(0.1, 0.2, 0.3))
})

# ------------------------------------------------------------------------
# 7) Row-only access and other non-column-selector forms fall through
# ------------------------------------------------------------------------
test_that("row-only access does NOT fire the helpful error", {
	df <- .make_catt_df()
	out <- df[1, ]
	expect_s3_class(out, "data.frame")
	expect_equal(nrow(out), 1L)
	expect_equal(out$cohort, "2")
})

test_that("numeric column selection does NOT fire the helpful error", {
	df <- .make_catt_df()
	expect_equal(df[, 2], c(0.1, 0.2, 0.3))
	expect_equal(df[[2]], c(0.1, 0.2, 0.3))
})

test_that("multi-column new-name selection does NOT fire the helpful error", {
	df <- .make_catt_df()
	out <- df[, c("cohort", "estimate")]
	expect_named(out, c("cohort", "estimate"))
})

# ------------------------------------------------------------------------
# 8) print() falls through to print.data.frame (no custom print method)
# ------------------------------------------------------------------------
test_that("print(catt_df) dispatches to print.data.frame and includes the snake_case column names", {
	df <- .make_catt_df()
	out <- paste(capture.output(print(df)), collapse = "\n")
	expect_match(out, "cohort", fixed = TRUE)
	expect_match(out, "estimate", fixed = TRUE)
	expect_match(out, "p_value", fixed = TRUE)
	# Old names must not appear in the printed output.
	expect_false(grepl("Estimated TE", out, fixed = TRUE))
	expect_false(grepl("ConfIntLow", out, fixed = TRUE))
	expect_false(grepl("ConfIntHigh", out, fixed = TRUE))
})

# ------------------------------------------------------------------------
# 9) The result of fetwfe()$catt_df carries the class
# ------------------------------------------------------------------------
test_that("fetwfe()$catt_df carries class catt_df and the snake_case columns", {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	res <- fetwfeWithSimulatedData(sim, verbose = FALSE)
	expect_s3_class(res$catt_df, "catt_df")
	expect_s3_class(res$catt_df, "data.frame")
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
})

test_that("etwfe()$catt_df carries class catt_df and the snake_case columns (no `selected`)", {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	res <- etwfeWithSimulatedData(sim, verbose = FALSE)
	expect_s3_class(res$catt_df, "catt_df")
	expect_s3_class(res$catt_df, "data.frame")
	expect_named(
		res$catt_df,
		c("cohort", "estimate", "se", "ci_low", "ci_high", "p_value")
	)
})
