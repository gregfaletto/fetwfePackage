# Tests for the broom-package S3 methods registered in R/broom_methods.R.
#
# The generics `tidy` / `glance` / `augment` come from the `generics` package
# (which `broom` re-exports). The package itself only imports them; user code
# is expected to call them via `broom::` (or after `library(broom)`).
# Tests skip when `broom` isn't installed.

if (!requireNamespace("broom", quietly = TRUE)) {
	testthat::skip("broom not installed")
}

# ------------------------------------------------------------------------------
# Shared fixtures
# ------------------------------------------------------------------------------

# Reuse the small simulated regime from the existing class-method tests.
# Cohorts in the simulator's default setup adopt at times >= 2, so no unit is
# treated in period 1 and idCohorts() drops no rows; therefore
# nrow(sim$pdata) == res$N * res$T == nrow(res[$internal]$X_ints).
.simulated_setup <- function(seed = 20260515) {
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = seed
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	list(coefs = coefs, sim = sim)
}

.fetwfe_fixture <- function() {
	setup <- .simulated_setup()
	fetwfeWithSimulatedData(setup$sim)
}

.etwfe_fixture <- function() {
	setup <- .simulated_setup()
	etwfeWithSimulatedData(setup$sim)
}

.betwfe_fixture <- function() {
	setup <- .simulated_setup()
	betwfeWithSimulatedData(setup$sim)
}

# ------------------------------------------------------------------------------
# tidy.<class>
# ------------------------------------------------------------------------------

test_that("tidy.fetwfe returns broom-schema columns with selected", {
	res <- .fetwfe_fixture()
	td <- broom::tidy(res)
	expect_s3_class(td, "data.frame")
	expect_named(
		td,
		c(
			"term",
			"estimate",
			"std.error",
			"statistic",
			"p.value",
			"conf.low",
			"conf.high",
			"selected"
		)
	)
	expect_equal(nrow(td), res$R + 1L)
	expect_equal(td$term[1], "ATT")
	expect_match(td$term[-1], "^Cohort ")
	# Cohort rows are sorted ascending by cohort label.
	cohort_labels <- sub("^Cohort ", "", td$term[-1])
	expect_equal(cohort_labels, sort(cohort_labels))
})

test_that("tidy.etwfe has no selected column (no regularized selection)", {
	res <- .etwfe_fixture()
	td <- broom::tidy(res)
	expect_false("selected" %in% names(td))
	expect_named(
		td,
		c(
			"term",
			"estimate",
			"std.error",
			"statistic",
			"p.value",
			"conf.low",
			"conf.high"
		)
	)
	expect_equal(nrow(td), res$R + 1L)
})

test_that("tidy.betwfe has selected column", {
	res <- .betwfe_fixture()
	td <- broom::tidy(res)
	expect_true("selected" %in% names(td))
	expect_equal(nrow(td), res$R + 1L)
})

test_that("broom::tidy(..., conf.int = FALSE) omits CI columns", {
	res <- .fetwfe_fixture()
	td <- broom::tidy(res, conf.int = FALSE)
	expect_false("conf.low" %in% names(td))
	expect_false("conf.high" %in% names(td))
	# Other broom columns remain.
	expect_true(all(
		c("term", "estimate", "std.error", "statistic", "p.value") %in%
			names(td)
	))
})

test_that("broom::tidy() respects custom conf.level", {
	res <- .fetwfe_fixture()
	td_default <- broom::tidy(res, conf.level = 0.95)
	td_wide <- broom::tidy(res, conf.level = 0.99)
	# 99% CI must be at least as wide as 95% CI for every row.
	expect_true(all(
		(td_wide$conf.high - td_wide$conf.low) >=
			(td_default$conf.high - td_default$conf.low) - 1e-10
	))
})

test_that("broom::tidy() propagates NA std.error to statistic and p.value", {
	res <- .fetwfe_fixture()
	# Force a row to have NA SE and verify propagation.
	res$catt_df$SE[1] <- NA_real_
	res$catt_df$P_value[1] <- NA_real_
	td <- broom::tidy(res)
	expect_true(is.na(td$std.error[2]))
	expect_true(is.na(td$statistic[2]))
	expect_true(is.na(td$p.value[2]))
})

# ------------------------------------------------------------------------------
# glance.<class>
# ------------------------------------------------------------------------------

test_that("glance.fetwfe returns 13 columns including lambda_star", {
	res <- .fetwfe_fixture()
	gl <- broom::glance(res)
	expect_s3_class(gl, "data.frame")
	expect_equal(nrow(gl), 1L)
	expect_named(
		gl,
		c(
			"nobs",
			"n_units",
			"n_periods",
			"n_cohorts",
			"n_covs",
			"n_features",
			"lambda_star",
			"lambda_star_model_size",
			"sig_eps_sq",
			"sig_eps_c_sq",
			"alpha",
			"se_type",
			"indep_counts_used"
		)
	)
	expect_equal(gl$nobs, res$N * res$T)
	expect_equal(gl$n_units, res$N)
	expect_equal(gl$n_periods, res$T)
})

test_that("glance.etwfe omits lambda_star columns (11 columns)", {
	res <- .etwfe_fixture()
	gl <- broom::glance(res)
	expect_equal(nrow(gl), 1L)
	expect_false("lambda_star" %in% names(gl))
	expect_false("lambda_star_model_size" %in% names(gl))
	expect_named(
		gl,
		c(
			"nobs",
			"n_units",
			"n_periods",
			"n_cohorts",
			"n_covs",
			"n_features",
			"sig_eps_sq",
			"sig_eps_c_sq",
			"alpha",
			"se_type",
			"indep_counts_used"
		)
	)
})

test_that("glance.betwfe matches fetwfe schema (13 columns)", {
	res <- .betwfe_fixture()
	gl <- broom::glance(res)
	expect_equal(ncol(gl), 13L)
	expect_true("lambda_star" %in% names(gl))
})

# ------------------------------------------------------------------------------
# augment.<class>
# ------------------------------------------------------------------------------

test_that("augment.fetwfe returns data + .fitted + .resid with reconstruction", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	aug <- broom::augment(res, data = setup$sim$pdata, response = "y")
	expect_s3_class(aug, "data.frame")
	expect_equal(nrow(aug), nrow(setup$sim$pdata))
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	expect_true(is.numeric(aug$.fitted))
	expect_true(is.numeric(aug$.resid))
	# .fitted + .resid reconstructs the original (uncentered) response.
	expect_equal(aug$.fitted + aug$.resid, aug$y, tolerance = 1e-8)
})

test_that("augment errors without data argument", {
	res <- .fetwfe_fixture()
	expect_error(broom::augment(res), "data")
})

test_that("augment errors without response argument", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	expect_error(broom::augment(res, data = setup$sim$pdata), "response")
})

test_that("augment errors on row-count mismatch with first-period hint", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	bad_data <- setup$sim$pdata[1:10, , drop = FALSE]
	expect_error(
		broom::augment(res, data = bad_data, response = "y"),
		"first time period"
	)
})

test_that("augment errors when response column not found", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	expect_error(
		broom::augment(res, data = setup$sim$pdata, response = "nonexistent"),
		"not found"
	)
})

test_that("augment.etwfe and augment.betwfe work with flat X_ints slot", {
	setup <- .simulated_setup()
	res_etwfe <- etwfeWithSimulatedData(setup$sim)
	aug_etwfe <- broom::augment(
		res_etwfe,
		data = setup$sim$pdata,
		response = "y"
	)
	expect_equal(
		aug_etwfe$.fitted + aug_etwfe$.resid,
		aug_etwfe$y,
		tolerance = 1e-8
	)

	res_betwfe <- betwfeWithSimulatedData(setup$sim)
	aug_betwfe <- broom::augment(
		res_betwfe,
		data = setup$sim$pdata,
		response = "y"
	)
	expect_equal(
		aug_betwfe$.fitted + aug_betwfe$.resid,
		aug_betwfe$y,
		tolerance = 1e-8
	)
})

# ------------------------------------------------------------------------------
# tidy.fetwfe_event_study and tidy.FETWFE_tes
# ------------------------------------------------------------------------------

test_that("tidy.fetwfe_event_study returns broom-schema columns", {
	res <- .fetwfe_fixture()
	es <- event_study(res)
	td <- broom::tidy(es)
	expect_s3_class(td, "data.frame")
	expect_named(
		td,
		c(
			"term",
			"event_time",
			"n_cohorts",
			"estimate",
			"std.error",
			"statistic",
			"p.value",
			"conf.low",
			"conf.high"
		)
	)
	expect_equal(nrow(td), nrow(es))
	expect_match(td$term, "^e\\d+$")
})

test_that("tidy.FETWFE_tes has NA_real_ for SE columns and integer cohort terms", {
	setup <- .simulated_setup()
	tes <- getTes(setup$coefs)
	td <- broom::tidy(tes)
	expect_equal(nrow(td), length(tes$actual_cohort_tes) + 1L)
	expect_equal(td$term[1], "ATT_true")
	expect_match(td$term[-1], "^Cohort \\d+$")
	# All SE-related columns are numeric NA, not logical NA.
	expect_true(is.numeric(td$std.error))
	expect_true(is.numeric(td$p.value))
	expect_true(all(is.na(td$std.error)))
	expect_true(all(is.na(td$statistic)))
	expect_true(all(is.na(td$p.value)))
	expect_true(all(is.na(td$conf.low)))
	expect_true(all(is.na(td$conf.high)))
	# Estimates are the real true ATT + cohort effects.
	expect_equal(td$estimate[1], tes$att_true)
	expect_equal(td$estimate[-1], tes$actual_cohort_tes)
})

# ------------------------------------------------------------------------------
# Dispatch
# ------------------------------------------------------------------------------

test_that("broom::tidy / broom::glance / broom::augment dispatch to package methods", {
	skip_if_not_installed("broom")
	res <- .fetwfe_fixture()
	# `broom::tidy(res)` should hit `tidy.fetwfe` via S3 dispatch.
	td_via_broom <- broom::tidy(res)
	td_direct <- broom::tidy(res)
	expect_identical(td_via_broom, td_direct)
})
