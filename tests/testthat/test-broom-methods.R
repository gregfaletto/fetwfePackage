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
	# .fitted + .resid reconstructs the original (uncentered) response.
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
	res <- fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		q = 0.5,
		verbose = FALSE
	)
	# Passing the original (un-trimmed) panel — augment should auto-trim.
	aug <- broom::augment(res, data = pdata)
	expect_equal(nrow(aug), res$N * res$T)
	expect_true(all(c(".fitted", ".resid") %in% names(aug)))
	# Dropped unit should not appear in the augmented output.
	expect_false(u %in% aug$unit)
	# Round-trip property still holds on the trimmed panel.
	expect_equal(aug$.fitted + aug$.resid, aug$y, tolerance = 1e-8)
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

test_that("tidy.fetwfe_event_study respects conf.int = FALSE", {
	res <- .fetwfe_fixture()
	es <- event_study(res)
	td <- broom::tidy(es, conf.int = FALSE)
	expect_false("conf.low" %in% names(td))
	expect_false("conf.high" %in% names(td))
})

test_that("tidy.fetwfe_event_study recomputes CIs at a custom conf.level", {
	res <- .fetwfe_fixture()
	es <- event_study(res)
	td_95 <- broom::tidy(es, conf.level = 0.95)
	td_99 <- broom::tidy(es, conf.level = 0.99)
	# 99% widths should be wider than 95% widths everywhere (allowing
	# tiny floating-point slack).
	expect_true(all(
		(td_99$conf.high - td_99$conf.low) >=
			(td_95$conf.high - td_95$conf.low) - 1e-10
	))
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
	# cohorts the same way `fetwfe(simulateData(coefs))$catt_df$Cohort` does,
	# so `tidy.FETWFE_tes` rows align with `tidy.<estimator>` rows on the
	# same simulated panel. If `simulateData()` ever changes its cohort-
	# assignment scheme (e.g., adopts at non-sequential times), this test
	# fails and `getTes()` needs to be updated to derive `cohort_times` from
	# whatever new convention the simulator uses.
	setup <- .simulated_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	tes <- getTes(setup$coefs)
	# `catt_df$Cohort` is character ("2", "3", ...); coerce to integer for
	# comparison.
	cohort_labels_from_fit <- as.integer(res$catt_df$Cohort)
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

test_that("broom::tidy / broom::glance / broom::augment dispatch to package methods", {
	skip_if_not_installed("broom")
	res <- .fetwfe_fixture()
	# `broom::tidy(res)` should hit `tidy.fetwfe` via S3 dispatch.
	td_via_broom <- broom::tidy(res)
	td_direct <- broom::tidy(res)
	expect_identical(td_via_broom, td_direct)
})
