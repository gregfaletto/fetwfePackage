library(testthat)
library(fetwfe)

# Issue #28: formula interface for `covs`. The four estimators now accept
# `covs = ~ x1 + x2` (one-sided formula) in addition to the long-standing
# `covs = c("x1", "x2")` (character vector). Internal helper
# `.process_covs_input()` normalizes formula -> character vector; the
# rest of the package assumes the character form.

# ----------------------------------------------------------------------
# 1) Helper-level: character vectors are returned unchanged.
# ----------------------------------------------------------------------

test_that(".process_covs_input() returns character vectors unchanged", {
	expect_identical(
		fetwfe:::.process_covs_input(c("x1", "x2")),
		c("x1", "x2")
	)
	expect_identical(fetwfe:::.process_covs_input(c()), character(0))
	expect_identical(
		fetwfe:::.process_covs_input(character(0)),
		character(0)
	)
})

# ----------------------------------------------------------------------
# 2) Helper-level: NULL maps to character(0).
# ----------------------------------------------------------------------

test_that(".process_covs_input(NULL) returns character(0)", {
	expect_identical(fetwfe:::.process_covs_input(NULL), character(0))
})

# ----------------------------------------------------------------------
# 3) Helper-level: one-sided formulas are normalized to their term labels.
# ----------------------------------------------------------------------

test_that("one-sided formulas normalize to term labels", {
	expect_identical(
		fetwfe:::.process_covs_input(~x1),
		"x1"
	)
	expect_identical(
		fetwfe:::.process_covs_input(~ x1 + x2),
		c("x1", "x2")
	)
	expect_identical(
		fetwfe:::.process_covs_input(~ x1 + x2 + x3),
		c("x1", "x2", "x3")
	)
	# `~ 1` and `~ 0` are "no covariates" spellings -> character(0).
	expect_identical(fetwfe:::.process_covs_input(~1), character(0))
	expect_identical(fetwfe:::.process_covs_input(~0), character(0))
})

# ----------------------------------------------------------------------
# 4) Helper-level: rejections with structured messages.
# ----------------------------------------------------------------------

test_that("two-sided formulas are rejected", {
	expect_error(
		fetwfe:::.process_covs_input(y ~ x1 + x2),
		"one-sided"
	)
})

test_that("interactions and transformations are rejected with a clear pointer", {
	# Interactions: `x1:x2` produces term label `x1:x2` (with `:`)
	expect_error(
		fetwfe:::.process_covs_input(~ x1:x2),
		"additive bare variable names"
	)
	# `x1 * x2` expands to `x1 + x2 + x1:x2`; the `:` term triggers.
	expect_error(
		fetwfe:::.process_covs_input(~ x1 * x2),
		"additive bare variable names"
	)
	# I() transformations
	expect_error(
		fetwfe:::.process_covs_input(~ I(x1^2)),
		"additive bare variable names"
	)
})

test_that("non-character, non-formula input is rejected with class info", {
	expect_error(
		fetwfe:::.process_covs_input(123),
		"character vector.*one-sided formula"
	)
	expect_error(
		fetwfe:::.process_covs_input(list("x1")),
		"character vector.*one-sided formula"
	)
})

# ----------------------------------------------------------------------
# 5) End-to-end: formula and character forms produce identical fits on
#    the same simulated panel. Each of the four estimators tested.
# ----------------------------------------------------------------------

# Build a small simulated panel WITH d > 0 so covs actually matters.
.covs_setup <- function() {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(coefs, N = 60, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	# simulateData() produces pdata with covariate columns named "X1", "X2".
	# Extract for direct fits.
	dat
}

test_that("fetwfe(): formula and character forms produce identical results", {
	sim <- .covs_setup()
	df <- sim$pdata
	# fetwfeWithSimulatedData() goes through fetwfe() with the sim's
	# preset covs. Run direct fetwfe() for the comparison.
	covs_char <- c("cov1", "cov2")
	covs_form <- ~ cov1 + cov2
	res_char <- fetwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = covs_char,
		q = 0.5
	)
	res_form <- fetwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = covs_form,
		q = 0.5
	)
	expect_equal(res_char$att_hat, res_form$att_hat)
	expect_equal(res_char$att_se, res_form$att_se)
	expect_equal(res_char$catt_df$estimate, res_form$catt_df$estimate)
	expect_equal(res_char$catt_df$se, res_form$catt_df$se)
})

test_that("etwfe(): formula and character forms produce identical results", {
	sim <- .covs_setup()
	df <- sim$pdata
	res_char <- etwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2")
	)
	res_form <- etwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = ~ cov1 + cov2
	)
	expect_equal(res_char$att_hat, res_form$att_hat)
	expect_equal(res_char$att_se, res_form$att_se)
	expect_equal(res_char$catt_df$estimate, res_form$catt_df$estimate)
})

test_that("betwfe(): formula and character forms produce identical results", {
	sim <- .covs_setup()
	df <- sim$pdata
	res_char <- betwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		q = 0.5
	)
	res_form <- betwfe(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = ~ cov1 + cov2,
		q = 0.5
	)
	expect_equal(res_char$att_hat, res_form$att_hat)
	expect_equal(res_char$att_se, res_form$att_se)
})

test_that("twfeCovs(): formula and character forms produce identical results", {
	sim <- .covs_setup()
	df <- sim$pdata
	res_char <- twfeCovs(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2")
	)
	res_form <- twfeCovs(
		df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = ~ cov1 + cov2
	)
	expect_equal(res_char$att_hat, res_form$att_hat)
	expect_equal(res_char$att_se, res_form$att_se)
})

# ----------------------------------------------------------------------
# 6) End-to-end: bad formula at the estimator boundary is rejected before
#    any expensive computation.
# ----------------------------------------------------------------------

test_that("invalid covs formula is rejected at the estimator entry", {
	sim <- .covs_setup()
	df <- sim$pdata
	expect_error(
		fetwfe(
			df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = ~ cov1 * cov2
		),
		"additive bare variable names"
	)
	expect_error(
		fetwfe(
			df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = y ~ cov1 + cov2
		),
		"one-sided"
	)
})
