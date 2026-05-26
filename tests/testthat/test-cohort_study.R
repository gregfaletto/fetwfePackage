# Tests for cohortStudy() and tidy.cohortStudy() — issue #137.
#
# cohortStudy(result) is a function-style accessor parallel to
# eventStudy(): it returns result$catt_df with class re-assigned to
# c("cohortStudy", "catt_df", "data.frame"). The new "cohortStudy" class
# dispatches the broom tidier (tidy.cohortStudy); the "catt_df" class
# preserves the #136 helpful-error layer for pre-1.11.0 column-name
# access; the "data.frame" base preserves standard df ops.

library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------
# Shared fixture: a small simulated panel exercising all four estimator
# entry points. Cohorts in the simulator's default setup adopt at times
# >= 2, so no unit is treated in period 1 and idCohorts() drops no rows.
# ------------------------------------------------------------------------
.cs_setup <- function(seed = 20260525) {
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

# ------------------------------------------------------------------------
# 1) Class hierarchy: cohortStudy returns the chained class vector.
# ------------------------------------------------------------------------
test_that("cohortStudy returns class c('cohortStudy', 'catt_df', 'data.frame')", {
	setup <- .cs_setup()
	for (fit_fn in list(
		fetwfeWithSimulatedData,
		etwfeWithSimulatedData,
		betwfeWithSimulatedData
	)) {
		res <- fit_fn(setup$sim)
		cs <- cohortStudy(res)
		expect_s3_class(cs, "cohortStudy")
		expect_s3_class(cs, "catt_df")
		expect_s3_class(cs, "data.frame")
		# Class vector is exactly the three-element chain in the right order.
		expect_identical(
			class(cs),
			c("cohortStudy", "catt_df", "data.frame")
		)
	}
})

# ------------------------------------------------------------------------
# 2) Column shape per estimator: fetwfe / betwfe include `selected`;
#    etwfe / twfeCovs do not.
# ------------------------------------------------------------------------
test_that("cohortStudy column shape: fetwfe and betwfe include `selected`", {
	setup <- .cs_setup()
	res_f <- fetwfeWithSimulatedData(setup$sim)
	cs_f <- cohortStudy(res_f)
	expect_named(
		cs_f,
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

	res_b <- betwfeWithSimulatedData(setup$sim)
	cs_b <- cohortStudy(res_b)
	expect_named(
		cs_b,
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

test_that("cohortStudy column shape: etwfe has no `selected`", {
	setup <- .cs_setup()
	res_e <- etwfeWithSimulatedData(setup$sim)
	cs_e <- cohortStudy(res_e)
	expect_named(
		cs_e,
		c("cohort", "estimate", "se", "ci_low", "ci_high", "p_value")
	)
})

# ------------------------------------------------------------------------
# 3) Content identity: cohortStudy(result) is a passthrough on
#    result$catt_df modulo class. Values, column count, and row count
#    are unchanged.
# ------------------------------------------------------------------------
test_that("cohortStudy preserves content of result$catt_df", {
	setup <- .cs_setup()
	for (fit_fn in list(
		fetwfeWithSimulatedData,
		etwfeWithSimulatedData,
		betwfeWithSimulatedData
	)) {
		res <- fit_fn(setup$sim)
		cs <- cohortStudy(res)
		# Identical row count.
		expect_equal(nrow(cs), nrow(res$catt_df))
		# Identical column names.
		expect_equal(names(cs), names(res$catt_df))
		# Identical numeric content per column (stripping classes for the
		# comparison, since cs carries an additional `cohortStudy` class).
		expect_equal(unclass(cs), unclass(res$catt_df))
	}
})

# ------------------------------------------------------------------------
# 4) Column types are well-formed.
# ------------------------------------------------------------------------
test_that("cohortStudy column types are well-formed", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	# `cohort` is a character label (stored as character throughout the
	# package; see fetwfe()'s catt_df spec). Numeric coercion succeeds.
	expect_true(is.character(cs$cohort) || is.numeric(cs$cohort))
	expect_false(anyNA(suppressWarnings(as.numeric(cs$cohort))))
	# Numeric columns are numeric (may be NA for selected-out cohorts).
	expect_true(is.numeric(cs$estimate))
	expect_true(is.numeric(cs$se))
	expect_true(is.numeric(cs$ci_low))
	expect_true(is.numeric(cs$ci_high))
	expect_true(is.numeric(cs$p_value))
	# `selected` is logical.
	expect_true(is.logical(cs$selected))
})

# ------------------------------------------------------------------------
# 5) tidy(cohortStudy(result)): broom-shape translation.
# ------------------------------------------------------------------------
test_that("tidy.cohortStudy produces broom-shape output (fetwfe with selected)", {
	skip_if_not_installed("broom")
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	td <- broom::tidy(cs)
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
	expect_equal(nrow(td), nrow(cs))
	# term = paste0("cohort_", cohort)
	expect_equal(td$term, paste0("cohort_", cs$cohort))
	# estimate / std.error / conf.low / conf.high / p.value pass through.
	expect_equal(td$estimate, cs$estimate)
	expect_equal(td$std.error, cs$se)
	expect_equal(td$conf.low, cs$ci_low)
	expect_equal(td$conf.high, cs$ci_high)
	expect_equal(td$p.value, cs$p_value)
	# statistic = estimate / std.error when se > 0; NA when se == 0 or NA.
	# Mirrors the same NA-propagation idiom used by tidy.eventStudy().
	for (i in seq_len(nrow(td))) {
		if (is.na(cs$se[i]) || cs$se[i] == 0) {
			expect_true(is.na(td$statistic[i]))
		} else {
			expect_equal(td$statistic[i], cs$estimate[i] / cs$se[i])
		}
	}
})

test_that("tidy.cohortStudy omits `selected` for etwfe (no selection step)", {
	skip_if_not_installed("broom")
	setup <- .cs_setup()
	res <- etwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	td <- broom::tidy(cs)
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
})

test_that("tidy.cohortStudy is registered and dispatches via broom::tidy", {
	skip_if_not_installed("broom")
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	# Method-name lookup succeeds (registration in NAMESPACE).
	expect_false(
		is.null(getS3method("tidy", "cohortStudy", optional = TRUE))
	)
	# Actual dispatch through the broom generic produces a data.frame.
	# A missing S3method() line would route to the broom default and
	# return a character error message — this is the load-bearing check.
	cs <- cohortStudy(res)
	expect_s3_class(broom::tidy(cs), "data.frame")
})

# ------------------------------------------------------------------------
# 6) Helpful-error layer from #136 still fires through the chained class.
#    The cohortStudy class vector includes "catt_df", so [[, $, [ all
#    intercept pre-1.11.0 Title-Case column names.
# ------------------------------------------------------------------------
test_that("cohortStudy()$Cohort (old name) fires the migration error", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	expect_error(
		cs[["Cohort"]],
		"`Cohort` was renamed to `cohort` in fetwfe 1.11.0",
		fixed = TRUE
	)
	expect_error(
		cs[["Estimated TE"]],
		"`Estimated TE` was renamed to `estimate` in fetwfe 1.11.0",
		fixed = TRUE
	)
	expect_error(
		cs$SE,
		"`SE` was renamed to `se` in fetwfe 1.11.0",
		fixed = TRUE
	)
	expect_error(
		cs[, "P_value"],
		"`P_value` was renamed to `p_value` in fetwfe 1.11.0",
		fixed = TRUE
	)
})

test_that("cohortStudy()$cohort (new name) works without firing the migration error", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	# New-name access succeeds and returns the column unchanged.
	expect_equal(cs$cohort, res$catt_df$cohort)
	expect_equal(cs[["estimate"]], res$catt_df$estimate)
})

# ------------------------------------------------------------------------
# 7) Input-validation: error on non-fitted-class input.
# ------------------------------------------------------------------------
test_that("cohortStudy() errors on non-fitted-class input", {
	expect_error(
		cohortStudy(list(foo = 1)),
		"fetwfe.*etwfe.*betwfe.*twfeCovs"
	)
	expect_error(
		cohortStudy("not a fit"),
		"fetwfe.*etwfe.*betwfe.*twfeCovs"
	)
	expect_error(
		cohortStudy(data.frame(x = 1)),
		"fetwfe.*etwfe.*betwfe.*twfeCovs"
	)
	expect_error(
		cohortStudy(NULL),
		"fetwfe.*etwfe.*betwfe.*twfeCovs"
	)
})

# ------------------------------------------------------------------------
# 8) Input-validation: error when the fitted object has no catt_df slot.
# ------------------------------------------------------------------------
test_that("cohortStudy() errors when catt_df slot is missing", {
	# Construct a synthetic object that passes the class check but has no
	# catt_df slot. The post-construction check should fire.
	fake <- list(att_hat = 0)
	class(fake) <- "fetwfe"
	expect_error(
		cohortStudy(fake),
		"no `catt_df` slot",
		fixed = TRUE
	)

	# Also fires when catt_df is explicitly NULL on an otherwise-valid fit.
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	res$catt_df <- NULL
	expect_error(
		cohortStudy(res),
		"no `catt_df` slot",
		fixed = TRUE
	)
})

# ------------------------------------------------------------------------
# 9) cohortStudy works on twfeCovs fits.
# ------------------------------------------------------------------------
test_that("cohortStudy works on a twfeCovs fit", {
	setup <- .cs_setup()
	res_t <- suppressWarnings(twfeCovsWithSimulatedData(setup$sim))
	cs_t <- cohortStudy(res_t)
	expect_s3_class(cs_t, "cohortStudy")
	expect_s3_class(cs_t, "catt_df")
	expect_s3_class(cs_t, "data.frame")
	# twfeCovs is an OLS-style estimator: no `selected` column.
	expect_named(
		cs_t,
		c("cohort", "estimate", "se", "ci_low", "ci_high", "p_value")
	)
	expect_equal(nrow(cs_t), res_t$R)
})

# ------------------------------------------------------------------------
# 10) Standard data.frame operations still work (print, head, nrow,
#     dim, dplyr::filter-like row selection by logical index).
# ------------------------------------------------------------------------
test_that("standard data.frame operations work on cohortStudy output", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	# nrow / ncol / dim.
	expect_equal(nrow(cs), 3L)
	expect_equal(ncol(cs), 7L)
	expect_equal(dim(cs), c(3L, 7L))
	# head() falls through to head.data.frame.
	expect_s3_class(head(cs, 2L), "data.frame")
	expect_equal(nrow(head(cs, 2L)), 2L)
	# Logical row selection (the form `cs[cs$selected, ]`) does not fire
	# the helpful-error layer (the row selector is logical, not a
	# Title-Case name).
	sub <- cs[cs$selected, ]
	expect_s3_class(sub, "data.frame")
	expect_true(nrow(sub) <= nrow(cs))
})

# ------------------------------------------------------------------------
# 11) Robustness guards: cohortStudy() rejects a non-data.frame catt_df
#     slot, and tidy.cohortStudy() localizes the error when required
#     columns are missing. Both paths require manual user mutation of
#     the fitted object's slot; no realistic estimator path produces
#     them, but the guards localize an otherwise cryptic failure.
# ------------------------------------------------------------------------
test_that("cohortStudy() rejects a non-data.frame catt_df slot", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	# Mutate the slot to a list; verify the helpful stop fires.
	res$catt_df <- list(estimate = 1, se = 0.1)
	expect_error(
		cohortStudy(res),
		"`catt_df` slot is not a data.frame"
	)
	# Mutate to an atomic vector; same guard fires.
	res$catt_df <- c(1, 2, 3)
	expect_error(
		cohortStudy(res),
		"`catt_df` slot is not a data.frame"
	)
})

test_that("tidy.cohortStudy() localizes error when required columns are missing", {
	setup <- .cs_setup()
	res <- fetwfeWithSimulatedData(setup$sim)
	cs <- cohortStudy(res)
	# Drop a required column; the tidy translation should report which
	# column is missing rather than producing a cryptic "differing number
	# of rows" error from `data.frame()`.
	cs_broken <- cs
	cs_broken$ci_low <- NULL
	expect_error(
		broom::tidy(cs_broken),
		"missing required columns: ci_low"
	)
	# Drop two required columns; both are listed in the error message.
	cs_broken2 <- cs
	cs_broken2$ci_low <- NULL
	cs_broken2$ci_high <- NULL
	expect_error(
		broom::tidy(cs_broken2),
		"missing required columns: ci_low, ci_high"
	)
})
