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
	# Cohort rows are sorted ascending by NUMERIC cohort label. Note
	# `sort()` on character is lexicographic, which would give the wrong
	# answer for cohort labels ≥ 10 (see GitHub #53, fixed in v1.9.3 by
	# the composite numeric-first sort key). Assert numeric ordering
	# explicitly so this test catches any future regression on panels
	# with ≥ 10 cohorts.
	cohort_labels <- sub("^Cohort ", "", td$term[-1])
	expect_equal(as.numeric(cohort_labels), sort(as.numeric(cohort_labels)))
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
	res$catt_df$se[1] <- NA_real_
	res$catt_df$p_value[1] <- NA_real_
	td <- broom::tidy(res)
	expect_true(is.na(td$std.error[2]))
	expect_true(is.na(td$statistic[2]))
	expect_true(is.na(td$p.value[2]))
})

# ------------------------------------------------------------------------------
# glance.<class>
# ------------------------------------------------------------------------------

test_that("glance.fetwfe returns 16 columns including lambda_star and lambda_selection", {
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
			# v1.13.0 (#164): lambda-selection method provenance.
			"lambda_selection",
			"cv_folds",
			"cv_seed",
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

test_that("glance.betwfe matches fetwfe schema (16 columns)", {
	res <- .betwfe_fixture()
	gl <- broom::glance(res)
	expect_equal(ncol(gl), 16L)
	expect_true("lambda_star" %in% names(gl))
	expect_true("lambda_selection" %in% names(gl))
})

# ------------------------------------------------------------------------------
# augment.<class>
# ------------------------------------------------------------------------------

test_that("fitted objects carry y_mean and response_col_name slots", {
	setup <- .simulated_setup()
	for (fit_fn in list(
		fetwfeWithSimulatedData,
		etwfeWithSimulatedData,
		betwfeWithSimulatedData
	)) {
		res <- fit_fn(setup$sim)
		expect_true(is.numeric(res$y_mean))
		expect_equal(length(res$y_mean), 1L)
		expect_equal(res$y_mean, mean(setup$sim$pdata$y))
		expect_identical(res$response_col_name, "y")
	}
})

test_that("augment.fetwfe returns data + .fitted + .resid with reconstruction", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	aug <- broom::augment(res, data = setup$sim$pdata)
	expect_s3_class(aug, "data.frame")
	expect_equal(nrow(aug), nrow(setup$sim$pdata))
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	expect_true(is.numeric(aug$.fitted))
	expect_true(is.numeric(aug$.resid))
	# NOTE: the round-trip identity `.fitted + .resid == y` is tautological
	# — `.resid` is defined as `y - .fitted`, so this holds by construction
	# even when `.fitted` is glued to the wrong rows (the augment row-order
	# bug fixed in PR #48). The load-bearing correctness check for row
	# alignment lives in the "augment aligns rows to X_ints regardless of
	# input row order" block below. The assertion here is kept as a guard
	# against NaN propagation or factor coercion accidents, not as a
	# correctness check. See issue #57.
	expect_equal(aug$.fitted + aug$.resid, aug$y, tolerance = 1e-8)
})

test_that("augment errors without data argument", {
	res <- .fetwfe_fixture()
	expect_error(broom::augment(res), "data")
})

test_that("augment errors when fitted object lacks response_col_name slot", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	# Simulate a fitted object from a pre-1.9.0 dev build (no slot).
	res$response_col_name <- NULL
	expect_error(
		broom::augment(res, data = setup$sim$pdata),
		"response_col_name"
	)
})

test_that("augment errors when response column is missing from data", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	# Drop the response column from `data` so the lookup fails.
	data_without_y <- setup$sim$pdata
	data_without_y$y <- NULL
	expect_error(
		broom::augment(res, data = data_without_y),
		"response column 'y'"
	)
})

test_that("augment errors when data row count cannot be reconciled with fit", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	# 10 arbitrary rows from the panel won't form a balanced panel idCohorts
	# can process; the error propagates out of idCohorts.
	bad_data <- setup$sim$pdata[1:10, , drop = FALSE]
	expect_error(broom::augment(res, data = bad_data))
})

test_that("augment aligns rows to X_ints regardless of input row order", {
	# Latent-bug guard: if the user's `data` is in a different row order
	# from what the estimator sorted into when building X_ints, augment must
	# sort `data` to (unit, time) order so .fitted[i] corresponds to the
	# same observation as data[i, ]. Without this sort, .fitted gets glued
	# to the wrong rows and downstream interpretation is silently wrong.
	# (The `.fitted + .resid == y` round-trip holds by construction even
	# under misalignment, so a stronger check is needed.)
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	aug_orig <- broom::augment(res, data = setup$sim$pdata)
	set.seed(99)
	shuf_idx <- sample(nrow(setup$sim$pdata))
	aug_shuf <- broom::augment(res, data = setup$sim$pdata[shuf_idx, ])
	# Match shuffled output back to original by (unit, time) and check
	# .fitted values are identical.
	key_o <- paste(aug_orig$unit, aug_orig$time, sep = "_")
	key_s <- paste(aug_shuf$unit, aug_shuf$time, sep = "_")
	m <- match(key_o, key_s)
	expect_equal(aug_orig$.fitted, aug_shuf$.fitted[m])
})

test_that("augment works with factor covariates in the user's data", {
	# Factor covariates: the estimator expands them to dummies internally
	# via processFactors() during fitting; the user's pdata keeps the
	# original factor column. augment should produce correct fitted values
	# (computed from the pre-expanded X_ints) and return the panel with
	# the factor column intact.
	setup <- .simulated_setup()
	pdata <- setup$sim$pdata
	# Random factor with enough variation across units that the dummies
	# survive processCovs (which drops covariates whose first-period value
	# is constant across all units).
	set.seed(13)
	unit_to_group <- setNames(
		sample(c("a", "b", "c"), length(unique(pdata$unit)), replace = TRUE),
		unique(pdata$unit)
	)
	pdata$group <- factor(unit_to_group[pdata$unit])
	res <- fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2", "group"),
		q = 0.5,
		verbose = FALSE
	)
	# Stored covs slot keeps the ORIGINAL covs (pre-expansion), not the
	# dummy column names — so idCohorts in augment can find the right cols.
	expect_identical(res$covs, c("cov1", "cov2", "group"))
	aug <- broom::augment(res, data = pdata)
	expect_equal(nrow(aug), nrow(pdata))
	# Round-trip holds (sanity, not row-alignment proof — see preceding test).
	expect_equal(aug$.fitted + aug$.resid, aug$y, tolerance = 1e-8)
	# Factor column is preserved in the output.
	expect_true("group" %in% names(aug))
	expect_s3_class(aug$group, "factor")
})

test_that("augment auto-trims first-period-treated units when present in data", {
	setup <- .simulated_setup()
	# Take the simulated panel and mark one never-treated unit as treated at
	# time 1, simulating a real-data panel that needs idCohorts trimming.
	# The estimator will drop this unit at fit time; augment should drop the
	# same rows when given the original (pre-trim) panel.
	pdata <- setup$sim$pdata
	# Pick a unit that's never treated in the simulator output.
	never_treated_units <- unique(
		pdata$unit[
			vapply(
				unique(pdata$unit),
				function(u) {
					all(pdata$treatment[pdata$unit == u] == 0)
				},
				logical(1)
			)
		]
	)
	skip_if(
		length(never_treated_units) == 0,
		"no never-treated units in fixture"
	)
	u <- never_treated_units[1]
	# Treatment must be absorbing (once treated, always treated), so set it
	# to 1 for all of this unit's periods. The estimator will detect that
	# treatment began at time 1 and drop the unit.
	pdata$treatment[pdata$unit == u] <- 1L

	# Fit on the modified panel — fetwfe internally drops the period-1-treated
	# unit, so X_ints has (N - 1) * T rows.
	expect_warning(
		res <- fetwfe(
			pdata = pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2"),
			q = 0.5,
			verbose = FALSE
		),
		"treated in the first time period"
	)
	# Passing the original (un-trimmed) panel — augment should auto-trim.
	expect_warning(
		aug <- broom::augment(res, data = pdata),
		"treated in the first time period"
	)
	expect_equal(nrow(aug), res$N * res$T)
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	# Dropped unit should not appear in the augmented output.
	expect_false(u %in% aug$unit)
	# Round-trip property still holds on the trimmed panel.
	expect_equal(aug$.fitted + aug$.resid, aug$y, tolerance = 1e-8)
})

# ------------------------------------------------------------------------------
# Issue #116 Gap 3: the auto-trim test above asserts only the round-trip
# identity .fitted + .resid == y, which holds by construction (.resid is
# defined as y - .fitted) regardless of whether the rows are correctly
# aligned. The row-order-invariance test at line ~286 is the genuine
# alignment lock, but it runs on an UNTRIMMED panel. This test is the
# intersection: it runs the permutation-invariance check on an auto-trimmed
# fit, so a regression in augment's (unit, time) re-sort along the trim path
# is caught. It deliberately does NOT assert the .fitted + .resid == y
# tautology.
# ------------------------------------------------------------------------------
test_that("augment row order is invariant after auto-trim of first-period-treated units", {
	setup <- .simulated_setup()
	pdata <- setup$sim$pdata
	# Mark one never-treated unit as treated at time 1, so the estimator and
	# augment both auto-drop it (same pattern as the auto-trim test above).
	never_treated_units <- unique(
		pdata$unit[
			vapply(
				unique(pdata$unit),
				function(u) {
					all(pdata$treatment[pdata$unit == u] == 0)
				},
				logical(1)
			)
		]
	)
	skip_if(
		length(never_treated_units) == 0,
		"no never-treated units in fixture"
	)
	u <- never_treated_units[1]
	# Treatment is absorbing: set it to 1 for all of this unit's periods so
	# the estimator detects first-period treatment and drops the unit.
	pdata$treatment[pdata$unit == u] <- 1L

	# Fit on the modified panel; the trim warns, so suppress.
	res <- suppressWarnings(fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		q = 0.5,
		verbose = FALSE
	))

	# augment the original panel and a row-shuffled copy of the SAME panel;
	# augment re-warns via its internal idCohorts(), so suppress both.
	aug_orig <- suppressWarnings(broom::augment(res, data = pdata))
	set.seed(116)
	shuf_idx <- sample(nrow(pdata))
	aug_shuf <- suppressWarnings(
		broom::augment(res, data = pdata[shuf_idx, ])
	)

	# Match the shuffled output back to the original by (unit, time) and
	# require .fitted values to be identical -- the genuine row-order lock.
	key_o <- paste(aug_orig$unit, aug_orig$time, sep = "_")
	key_s <- paste(aug_shuf$unit, aug_shuf$time, sep = "_")
	m <- match(key_o, key_s)
	expect_false(anyNA(m))
	expect_equal(aug_orig$.fitted, aug_shuf$.fitted[m])
	# The auto-dropped unit must be absent from both outputs, and the row
	# count must equal the trimmed-fit dimensions.
	expect_false(u %in% aug_orig$unit)
	expect_false(u %in% aug_shuf$unit)
	expect_equal(nrow(aug_orig), res$N * res$T)
	expect_equal(nrow(aug_shuf), res$N * res$T)
})

test_that("augment.etwfe and augment.betwfe work with flat X_ints slot", {
	setup <- .simulated_setup()
	res_etwfe <- etwfeWithSimulatedData(setup$sim)
	aug_etwfe <- broom::augment(res_etwfe, data = setup$sim$pdata)
	expect_equal(
		aug_etwfe$.fitted + aug_etwfe$.resid,
		aug_etwfe$y,
		tolerance = 1e-8
	)

	res_betwfe <- betwfeWithSimulatedData(setup$sim)
	aug_betwfe <- broom::augment(res_betwfe, data = setup$sim$pdata)
	expect_equal(
		aug_betwfe$.fitted + aug_betwfe$.resid,
		aug_betwfe$y,
		tolerance = 1e-8
	)
})

# ------------------------------------------------------------------------------
# tidy.eventStudy and tidy.FETWFE_tes
# ------------------------------------------------------------------------------

test_that("tidy.eventStudy returns broom-schema columns", {
	res <- .fetwfe_fixture()
	es <- eventStudy(res)
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

test_that("tidy.eventStudy respects conf.int = FALSE", {
	res <- .fetwfe_fixture()
	es <- eventStudy(res)
	td <- broom::tidy(es, conf.int = FALSE)
	expect_false("conf.low" %in% names(td))
	expect_false("conf.high" %in% names(td))
})

test_that("tidy.eventStudy recomputes CIs at a custom conf.level", {
	res <- .fetwfe_fixture()
	es <- eventStudy(res)
	td_95 <- broom::tidy(es, conf.level = 0.95)
	td_99 <- broom::tidy(es, conf.level = 0.99)
	# Sharp check: the CI half-width is `qnorm(1 - alpha/2) * std.error`,
	# so the difference between 99% and 95% high bounds is
	# `(qnorm(0.995) - qnorm(0.975)) * std.error` at every row. The pre-
	# #178 parametric `width_99 >= width_95 - 1e-10` assertion held even
	# if `conf.level` were ignored entirely (two equal widths satisfy
	# `>=`); the sharp form pins the exact formula. `td_95$std.error`
	# is used because it equals `td_99$std.error` (`std.error` doesn't
	# depend on `conf.level`).
	z_diff <- stats::qnorm(0.995) - stats::qnorm(0.975)
	expect_equal(
		td_99$conf.high - td_95$conf.high,
		z_diff * td_95$std.error
	)
	expect_equal(
		td_95$conf.low - td_99$conf.low,
		z_diff * td_95$std.error
	)
})

test_that("tidy.eventStudy localizes error when required columns are missing", {
	# Parallels the guard PR #150 added to tidy.cohortStudy(). If a user
	# mutated the `eventStudy` frame to drop a required column, the call
	# should report which column is missing rather than producing a
	# cryptic "differing number of rows" error from `data.frame()`.
	res <- .fetwfe_fixture()
	es <- eventStudy(res)
	# Drop a single required column.
	es_broken <- es
	es_broken$se <- NULL
	expect_error(
		broom::tidy(es_broken),
		"missing required columns: se"
	)
	# Drop two required columns; both are listed in the error message.
	es_broken2 <- es
	es_broken2$se <- NULL
	es_broken2$p_value <- NULL
	expect_error(
		broom::tidy(es_broken2),
		"missing required columns: se, p_value"
	)
	# CI columns are NOT required (this method computes its own CIs from
	# estimate +/- z * se). Dropping them should NOT fire the guard.
	es_no_ci <- es
	es_no_ci$ci_low <- NULL
	es_no_ci$ci_high <- NULL
	expect_silent(broom::tidy(es_no_ci))
})

test_that("FETWFE_tes carries cohort_times slot (simulator convention)", {
	setup <- .simulated_setup()
	tes <- getTes(setup$coefs)
	expect_true(!is.null(tes$cohort_times))
	# Simulator convention: cohort r adopts at calendar time r + 1.
	expect_equal(tes$cohort_times, as.integer(seq_len(tes$R) + 1L))
})

test_that("FETWFE_tes cohort_times agree with simulateData()'s actual cohort assignment (assertion guarding the convention)", {
	# Cross-class invariant: `getTes(coefs)$cohort_times` is supposed to label
	# cohorts the same way `fetwfe(simulateData(coefs))$catt_df$cohort` does,
	# so `tidy.FETWFE_tes` rows align with `tidy.<estimator>` rows on the
	# same simulated panel. If `simulateData()` ever changes its cohort-
	# assignment scheme (e.g., adopts at non-sequential times), this test
	# fails and `getTes()` needs to be updated to derive `cohort_times` from
	# whatever new convention the simulator uses.
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	tes <- getTes(setup$coefs)
	# `catt_df$cohort` is character ("2", "3", ...); coerce to integer for
	# comparison.
	cohort_labels_from_fit <- as.integer(res$catt_df$cohort)
	expect_equal(cohort_labels_from_fit, tes$cohort_times)
})

test_that("tidy.FETWFE_tes labels cohorts by adoption time (matches estimator tidy)", {
	setup <- .simulated_setup()
	tes <- getTes(setup$coefs)
	td <- broom::tidy(tes)
	# Cohort labels should be the adoption times, not 1..R indices.
	expect_equal(td$term[-1], paste0("Cohort ", tes$cohort_times))
})

test_that("tidy.FETWFE_tes has NA_real_ for SE columns and Cohort <time> terms", {
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

test_that("broom S3 methods are registered AND dispatch correctly for fetwfe / etwfe / betwfe", {
	skip_if_not_installed("broom")
	# Replaces a prior tautology (calling `broom::tidy(res)` twice and
	# asserting identity) that passed even when no S3 method was
	# registered. Issue #57.
	#
	# This test exercises actual S3 dispatch through the broom generics
	# rather than just `getS3method()` lookup. The distinction matters:
	# `getS3method(generic, class, optional = TRUE)` finds a method by
	# name (`paste(generic, class, sep = ".")`) in the package namespace
	# even when the `S3method(generic, class)` line is missing from
	# NAMESPACE — so it does NOT catch missing-registration regressions.
	# Calling `broom::tidy(res)` and checking the result, by contrast,
	# routes through R's S3 dispatch table, which DOES use the NAMESPACE
	# registration. Mutation test: remove an `S3method(...)` line from
	# NAMESPACE; the corresponding `broom::tidy/glance/augment` call
	# below falls through to the broom default (returns a character
	# error message instead of a data.frame), and the `expect_s3_class`
	# assertion fails.
	setup <- .simulated_setup()
	fits <- list(
		fetwfe = fetwfeWithSimulatedData(setup$sim),
		etwfe = etwfeWithSimulatedData(setup$sim),
		betwfe = betwfeWithSimulatedData(setup$sim)
	)
	for (cls in names(fits)) {
		res <- fits[[cls]]
		expect_s3_class(broom::tidy(res), "data.frame")
		expect_s3_class(broom::glance(res), "data.frame")
		expect_s3_class(
			broom::augment(res, data = setup$sim$pdata),
			"data.frame"
		)
	}
	# Belt-and-suspenders: also confirm the method name lookups succeed
	# (catches renames even when registration is intact).
	for (cls in c("fetwfe", "etwfe", "betwfe")) {
		for (gen in c("tidy", "glance", "augment")) {
			expect_false(
				is.null(getS3method(gen, cls, optional = TRUE)),
				info = paste0("Missing S3 method: ", gen, ".", cls)
			)
		}
	}
	# event-study + cohort-level tidy methods.
	es <- eventStudy(fits$fetwfe)
	expect_s3_class(broom::tidy(es), "data.frame")
	expect_false(
		is.null(getS3method("tidy", "eventStudy", optional = TRUE))
	)
	expect_false(
		is.null(getS3method("tidy", "FETWFE_tes", optional = TRUE))
	)
})

# ------------------------------------------------------------------------------
# tidy.FETWFE_tes broom-convention args (#84 item 14)
# ------------------------------------------------------------------------------
test_that("tidy.FETWFE_tes accepts conf.int and conf.level (#84 item 14)", {
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 42
	)
	tes <- getTes(coefs)
	expect_s3_class(tes, "FETWFE_tes")

	# Default behavior: backward-compatible. CI columns included with NA values.
	td <- broom::tidy(tes)
	expect_true("conf.low" %in% names(td))
	expect_true("conf.high" %in% names(td))
	expect_true(all(is.na(td$conf.low)))
	expect_true(all(is.na(td$conf.high)))

	# conf.int = FALSE: CI columns dropped entirely.
	td_no_ci <- broom::tidy(tes, conf.int = FALSE)
	expect_false("conf.low" %in% names(td_no_ci))
	expect_false("conf.high" %in% names(td_no_ci))

	# conf.level accepted (validated but unused; output unchanged regardless).
	td_90 <- broom::tidy(tes, conf.int = TRUE, conf.level = 0.9)
	expect_identical(td, td_90)

	# conf.level is also ignored under conf.int = FALSE.
	td_no_ci_99 <- broom::tidy(tes, conf.int = FALSE, conf.level = 0.99)
	expect_identical(td_no_ci, td_no_ci_99)

	# Invalid conf.level errors via stopifnot.
	expect_error(broom::tidy(tes, conf.level = 1.5))
	expect_error(broom::tidy(tes, conf.level = 0))
})

# ------------------------------------------------------------------------------
# Item 5 (issue #84): augment.<class> dispatches successfully with a
# se_type = "cluster" fit. augment() only uses `beta_hat` (SE-independent),
# but the dispatch path with a cluster-SE fit hadn't been exercised before
# this test — the per-class augment.* methods route through the same
# `.augment_estimator_output` helper regardless of the SE route, so the
# correctness invariants (row alignment, `.fitted + .resid == y`) carry over.
# This test locks the dispatch path against future regressions that might
# accidentally couple augment() to the SE route.
# ------------------------------------------------------------------------------
test_that("augment.fetwfe dispatches with se_type = 'cluster' fit", {
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim, se_type = "cluster")
	aug <- broom::augment(res, data = setup$sim$pdata)
	expect_s3_class(aug, "data.frame")
	expect_equal(nrow(aug), nrow(setup$sim$pdata))
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	expect_true(is.numeric(aug$.fitted))
	expect_true(is.numeric(aug$.resid))
	expect_true(all(is.finite(aug$.fitted)))
})

# ------------------------------------------------------------------------------
# Item 6 (issue #84): augment.<class> dispatches successfully when the fit
# was constructed with `indep_counts` supplied. `fetwfeWithSimulatedData()`
# always passes the simulator's `indep_counts` through (the sim object
# carries one), so this fixture automatically exercises the indep-counts
# code path; the test pins augment() against any future regression that
# might cause indep-counts-bearing fits to lack the augment metadata slots.
# ------------------------------------------------------------------------------
test_that("augment.fetwfe dispatches when indep_counts was supplied", {
	setup <- .simulated_setup()
	# Sanity: the simulator's indep_counts is a valid integer vector that
	# fetwfeWithSimulatedData() will pass through.
	expect_true(is.numeric(setup$sim$indep_counts))
	expect_true(all(setup$sim$indep_counts >= 0))

	res <- fetwfeWithSimulatedData(setup$sim)
	# The fit carries the indep-counts-used flag set to TRUE.
	expect_true(isTRUE(res$indep_counts_used))

	aug <- broom::augment(res, data = setup$sim$pdata)
	expect_s3_class(aug, "data.frame")
	expect_equal(nrow(aug), nrow(setup$sim$pdata))
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	expect_true(is.numeric(aug$.fitted))
	expect_true(all(is.finite(aug$.fitted)))
})
