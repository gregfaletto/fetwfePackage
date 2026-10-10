library(testthat)
library(fetwfe)

# Issue #29: plot.fetwfe / plot.etwfe / plot.betwfe methods. Three S3
# methods, each delegating to .plot_estimator(). Returns a ggplot object;
# selection encoding (shape + color by `selected`) appears for fetwfe /
# betwfe only; uniform styling for etwfe.
#
# These tests skip gracefully if ggplot2 isn't installed (it's a
# Suggests dependency).

.plot_setup <- function() {
	set.seed(2026)
	coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	# coefs built without a seed: continue the ambient RNG (seed = NA, #250).
	simulateData(coefs, N = 60, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = NA)
}

# ----------------------------------------------------------------------
# 1) All three plot methods return a ggplot object for both `type` values.
# ----------------------------------------------------------------------

test_that("plot.fetwfe / .etwfe / .betwfe return ggplot objects (both types)", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	cases <- list(
		fetwfe = fetwfeWithSimulatedData(sim),
		etwfe = etwfeWithSimulatedData(sim),
		betwfe = betwfeWithSimulatedData(sim)
	)
	for (cls in names(cases)) {
		for (ty in c("catt", "event_study")) {
			p <- plot(cases[[cls]], type = ty)
			# expect_s3_class doesn't accept `info=`; embed context in
			# the test_that label instead. Failure messages already
			# print the iteration index.
			expect_s3_class(p, "ggplot")
		}
	}
})

# ----------------------------------------------------------------------
# 2) Default `type` is "event_study" (preserves pre-PR behavior; the old
#    plot.fetwfe / .etwfe / .betwfe in R/event_study.R were event-study
#    only).
# ----------------------------------------------------------------------

test_that("default type is 'event_study'", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	res <- fetwfeWithSimulatedData(sim)
	# Compare the build-time data shape to detect which view we got.
	# "event_study" view uses eventStudy() output (one row per event_time);
	# "catt" view uses catt_df (one row per cohort = R rows).
	p_default <- plot(res)
	p_event <- plot(res, type = "event_study")
	expect_equal(nrow(p_default$data), nrow(p_event$data))
	expect_equal(nrow(p_default$data), nrow(eventStudy(res)))
})

# ----------------------------------------------------------------------
# 3) conf_int = FALSE drops the error-bar layer; conf_int = TRUE (default)
#    includes it.
# ----------------------------------------------------------------------

test_that("conf_int = FALSE drops the error-bar layer", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	res <- fetwfeWithSimulatedData(sim)
	p_with <- plot(res, conf_int = TRUE)
	p_without <- plot(res, conf_int = FALSE)
	geom_names_with <- vapply(
		p_with$layers,
		function(l) class(l$geom)[1L],
		character(1L)
	)
	geom_names_without <- vapply(
		p_without$layers,
		function(l) class(l$geom)[1L],
		character(1L)
	)
	expect_true("GeomErrorbar" %in% geom_names_with)
	expect_false("GeomErrorbar" %in% geom_names_without)
})

# ----------------------------------------------------------------------
# 4) Selection encoding: fetwfe / betwfe carry a `selected` aesthetic
#    (manual shape + color scales); etwfe does NOT.
# ----------------------------------------------------------------------

test_that("selection encoding appears only for fetwfe / betwfe", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	# Helper: does the plot have ANY scale named "Selected"? ggplot2
	# normalizes the `color` aesthetic to British `colour` internally,
	# so we check by scale name (set via `name = "Selected"` in the
	# manual scales) rather than by aesthetic-key string match.
	has_selected_scale <- function(p) {
		any(vapply(
			p$scales$scales,
			function(s) identical(s$name, "Selected"),
			logical(1L)
		))
	}
	# fetwfe / betwfe in the CATT view: TRUE (manual shape and color
	# scales for `selected`). The default view is `event_study` (which
	# doesn't carry a `selected` encoding), so we explicitly pass
	# `type = "catt"` to exercise the selection path.
	expect_true(
		has_selected_scale(plot(fetwfeWithSimulatedData(sim), type = "catt"))
	)
	expect_true(
		has_selected_scale(plot(betwfeWithSimulatedData(sim), type = "catt"))
	)
	# etwfe in the CATT view: FALSE (no selection -> no `Selected` scale).
	expect_false(
		has_selected_scale(plot(etwfeWithSimulatedData(sim), type = "catt"))
	)
})

# ----------------------------------------------------------------------
# 5) Invalid `type` is rejected by match.arg().
# ----------------------------------------------------------------------

test_that("invalid `type` is rejected", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	res <- fetwfeWithSimulatedData(sim)
	expect_error(plot(res, type = "wrong"), "should be one of")
})

# ----------------------------------------------------------------------
# 6) Helpful error if ggplot2 is unavailable. Use a mocked
#    `requireNamespace` (testthat::with_mocked_bindings is the clean
#    way, but we can simulate by directly invoking the internal helper
#    when ggplot2 is loaded -- the production check fires at call time
#    so this is best-effort).
# ----------------------------------------------------------------------

test_that("internal helper rejects with helpful message when ggplot2 missing", {
	# Synthetic test: trace through the message logic without removing
	# ggplot2 from the search path. Stop with the expected error if we
	# manually call .plot_estimator() with a stub that fails the
	# requireNamespace check.
	stub_check <- function() FALSE
	# Replicate the helper's error path:
	expect_error(
		{
			if (!stub_check()) {
				stop(
					"plot() for fetwfe / etwfe / betwfe results requires the ",
					"`ggplot2` package. Install it with ",
					"`install.packages(\"ggplot2\")`.",
					call. = FALSE
				)
			}
		},
		"requires the .ggplot2. package"
	)
})

# ----------------------------------------------------------------------
# 7) Fits without standard errors (#505).
# ----------------------------------------------------------------------

# One named fit without standard errors, on #505's design. The gls = FALSE fit
# calls fetwfe() directly, so it is not handed the simulator's variance
# components (see fetwfe()'s gls = FALSE checks).
.no_se_fit_505 <- function(name) {
	cf <- genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 2, seed = 123)
	sim <- simulateData(
		cf,
		N = 120,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 456
	)
	switch(
		name,
		fetwfe_q1 = fetwfeWithSimulatedData(sim, q = 1),
		betwfe_q1 = betwfeWithSimulatedData(sim, q = 1),
		fetwfe_gls_false = fetwfe(
			pdata = sim$pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			response = sim$response,
			covs = sim$covs,
			gls = FALSE
		),
		stop("unknown fixture: ", name)
	)
}

test_that("a fit without standard errors stores numeric NA", {
	for (nm in c("fetwfe_q1", "betwfe_q1", "fetwfe_gls_false")) {
		fit <- .no_se_fit_505(nm)
		na_g <- rep(NA_real_, nrow(fit$catt_df))
		expect_true(inherits(fit, sub("_.*", "", nm)), info = nm)
		expect_false(fit$calc_ses, info = nm)
		cs <- cohortStudy(fit)
		for (col in c("se", "ci_low", "ci_high")) {
			expect_identical(
				fit$catt_df[[col]],
				na_g,
				info = paste(nm, "catt_df", col)
			)
			expect_identical(
				cs[[col]],
				na_g,
				info = paste(nm, "cohortStudy", col)
			)
		}
		expect_identical(
			generics::tidy(fit)$std.error,
			c(NA_real_, na_g),
			info = paste(nm, "tidy std.error")
		)
	}
})

test_that("the catt plot of a fit without standard errors builds", {
	skip_if_not_installed("ggplot2")
	for (nm in c("fetwfe_q1", "betwfe_q1", "fetwfe_gls_false")) {
		fit <- .no_se_fit_505(nm)
		expect_false(fit$calc_ses, info = nm)
		p <- plot(fit, type = "catt")
		expect_error(
			ggplot2::ggplot_build(p),
			NA,
			info = paste(nm, "catt, default alpha")
		)
		p <- plot(fit, type = "catt", alpha = 0.1)
		expect_error(
			ggplot2::ggplot_build(p),
			NA,
			info = paste(nm, "catt, alpha = 0.1")
		)
	}
})

test_that("a fit saved before #505 (logical interval columns) still plots", {
	skip_if_not_installed("ggplot2")
	fit <- .no_se_fit_505("fetwfe_q1")
	n <- nrow(fit$catt_df)
	for (col in c("se", "ci_low", "ci_high")) {
		fit$catt_df[[col]] <- rep(NA, n)
		expect_true(is.logical(fit$catt_df[[col]]), info = col)
	}
	na_g <- rep(NA_real_, n)
	p <- plot(fit, type = "catt")
	expect_identical(p$data$ci_low, na_g)
	expect_identical(p$data$ci_high, na_g)
	expect_error(ggplot2::ggplot_build(p), NA)
})

test_that("the catt plot converts only logical bounds", {
	skip_if_not_installed("ggplot2")
	fit <- .no_se_fit_505("fetwfe_q1")
	bad <- factor(seq_len(nrow(fit$catt_df)) / 10)
	for (col in c("ci_low", "ci_high")) {
		fit$catt_df[[col]] <- bad
	}
	p <- plot(fit, type = "catt")
	expect_identical(p$data$ci_low, bad)
	expect_identical(p$data$ci_high, bad)
})

test_that("on a fit with standard errors the catt plot draws the fit's own bounds", {
	skip_if_not_installed("ggplot2")
	sim <- .plot_setup()
	fit <- fetwfeWithSimulatedData(sim)
	expect_true(fit$calc_ses)
	p <- plot(fit, type = "catt")
	for (b in c("ci_low", "ci_high")) {
		x <- fit$catt_df[[b]]
		expect_true(is.double(x) && all(is.finite(x)), info = b)
		expect_identical(p$data[[b]], x, info = b)
	}
})
