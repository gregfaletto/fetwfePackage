# Tests for the `ci_type` fit-time argument and the simultaneous-by-default
# confidence bands (#197). The cohort-family / event-study-family band math is
# the #192 `simultaneousCIs()` machinery (correctness anchored by
# `test-simultaneous-cis.R`'s 500-rep MC coverage test + the
# `sqrt(diag(Sigma)) == se` identity); these tests lock the fit-time WIRING,
# the backward-compat `ci_type = "pointwise"` reproduction, the graceful
# degradation, the K = 1 bypass, the event-study NA-preservation blocker fix,
# and the Option-B tidy pass-through. Seeds are pinned; no PASS total pinned.

# Shared consecutive-cohort fixture (no degenerate event time): used by most
# tests. R = 4 so the cohort family has K = 4 >= 2 non-degenerate cohorts when
# eff_size is large enough that selection keeps them.
make_ci_panel <- function(
	seed = 101,
	R = 4,
	T = 6,
	d = 2,
	N = 200,
	eff_size = 2
) {
	coefs <- genCoefs(
		G = R,
		T = T,
		d = d,
		density = 0.5,
		eff_size = eff_size,
		seed = seed
	)
	simulateData(coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

# Scattered / short-overlap panel whose LAST event time is degenerate
# (cohorts adopting at calendar 3, 5, 6 + never-treated, T = 7 -> event time 5
# has no contributing cohort). Used by the blocker regression test (Test 11).
make_scattered_degenerate_panel <- function(n_per = 40, T = 7, seed = 909) {
	set.seed(seed)
	adopt <- c(3L, 5L, 6L, NA_integer_) # NA = never-treated
	rows <- vector("list", 0L)
	uid <- 0L
	for (a in adopt) {
		for (i in seq_len(n_per)) {
			uid <- uid + 1L
			for (t in seq_len(T)) {
				treated <- if (!is.na(a) && t >= a) 1L else 0L
				y <- stats::rnorm(1) + 0.5 * treated + 0.01 * uid
				rows[[length(rows) + 1L]] <- data.frame(
					unit = paste0("u", uid),
					year = t,
					treated = treated,
					y = y,
					stringsAsFactors = FALSE
				)
			}
		}
	}
	do.call(rbind, rows)
}

# AR(1)-corrupted panel for the se_type = "cluster" matrix cells: simulateData()
# produces F1-conforming (independent) shocks, so a clean panel's cluster band ~
# default. To make the cluster sandwich bite (and keep >= 2 finite-SE cohorts so
# the strict simultaneous > pointwise widening is assertable) we add a per-unit
# AR(1) shock on top of the simulated response. Mirrors test-fetwfe.R's
# make_ar1_panel but at R = 4 (the file's K = 4 convention).
make_ci_ar1_panel <- function(seed = 7, N = 150, rho = 0.85, sd_e = 1) {
	set.seed(seed)
	sim_coefs <- genCoefs(G = 4, T = 6, d = 2, density = 0.5, eff_size = 2)
	sim <- simulateData(sim_coefs, N = N, sig_eps_sq = 0.5, sig_eps_c_sq = 0.5)
	pdata <- sim$pdata
	pdata <- pdata[order(pdata$unit, pdata$time), ]
	units <- unique(pdata$unit)
	T_ <- length(unique(pdata$time))
	set.seed(seed + 100)
	for (u in units) {
		idx <- which(pdata$unit == u)
		ar <- numeric(T_)
		ar[1] <- stats::rnorm(1, sd = sd_e / sqrt(1 - rho^2))
		for (t in 2:T_) {
			ar[t] <- rho * ar[t - 1] + stats::rnorm(1, sd = sd_e)
		}
		pdata[idx, sim$response] <- pdata[idx, sim$response] + ar
	}
	sim$pdata <- pdata
	sim
}

# Smoke-level assertion shared by the interaction-matrix cells (Tests 14-19):
# the fit ran under the simultaneous default, the slot is set, and catt_df is
# well-formed (finite-SE rows bracket the estimate; degenerate rows are NA or 0
# per the established convention). Where >= 2 finite-SE cohorts survive, the
# simultaneous band is wider-or-equal vs a pointwise refit, and a pointwise
# refit reverts the bounds to estimate +/- qnorm(1 - alpha/2)*se. Behavioral
# (not well-formedness-only): the widen + revert checks are the locks.
expect_simultaneous_matrix_cell <- function(
	fit_s,
	fit_p,
	alpha = 0.05,
	label = ""
) {
	testthat::expect_identical(fit_s$ci_type, "simultaneous")
	testthat::expect_identical(fit_p$ci_type, "pointwise")
	z <- stats::qnorm(1 - alpha / 2)
	cd_s <- fit_s$catt_df
	cd_p <- fit_p$catt_df

	# (2) well-formed: finite-SE rows bracket the estimate; degenerate rows are
	# NA or 0 (the .build_selected_out_result() fallback writes 0).
	fin <- is.finite(cd_s$se) & cd_s$se > 0
	testthat::expect_true(all(
		cd_s$ci_low[fin] <= cd_s$estimate[fin] &
			cd_s$estimate[fin] <= cd_s$ci_high[fin]
	))
	deg <- !fin
	testthat::expect_true(all(
		is.na(cd_s$ci_low[deg]) | cd_s$ci_low[deg] == 0
	))

	# se unchanged across ci_type (only the multiplier differs).
	testthat::expect_equal(cd_s$se, cd_p$se)

	# (3) >= 2 finite-SE cohorts -> simultaneous wider-or-equal vs pointwise.
	if (sum(fin) >= 2L) {
		w_s <- (cd_s$ci_high - cd_s$ci_low)[fin]
		w_p <- (cd_p$ci_high - cd_p$ci_low)[fin]
		testthat::expect_true(all(w_s >= w_p - 1e-9))
	}

	# (4) pointwise refit reverts to estimate +/- z*se on finite-SE rows.
	if (any(fin)) {
		testthat::expect_equal(
			cd_p$ci_low[fin],
			(cd_p$estimate - z * cd_p$se)[fin],
			tolerance = 1e-12
		)
		testthat::expect_equal(
			cd_p$ci_high[fin],
			(cd_p$estimate + z * cd_p$se)[fin],
			tolerance = 1e-12
		)
	}
	invisible(NULL)
}

# ------------------------------------------------------------------------------
# Test 1: Default is simultaneous; bounds are wider than pointwise.
# ------------------------------------------------------------------------------
test_that("default ci_type is simultaneous and widens the catt_df bounds", {
	sim <- make_ci_panel()
	fit <- fetwfeWithSimulatedData(sim)
	expect_identical(fit$ci_type, "simultaneous")

	z <- stats::qnorm(1 - fit$alpha / 2)
	cd <- fit$catt_df
	# Restrict to non-degenerate cohorts (finite, positive SE) before asserting
	# strict widening: a selected-out cohort has se = 0 and an unchanged
	# point-mass interval.
	fin <- is.finite(cd$se) & cd$se > 0
	expect_gte(sum(fin), 2L) # fixture must keep >= 2 non-degenerate cohorts
	pw_low <- cd$estimate - z * cd$se
	pw_high <- cd$estimate + z * cd$se
	sim_width <- (cd$ci_high - cd$ci_low)[fin]
	pw_width <- (pw_high - pw_low)[fin]
	expect_true(all(sim_width > pw_width))

	es <- eventStudy(fit)
	es_fin <- is.finite(es$se) & es$se > 0
	expect_gte(sum(es_fin), 2L)
	es_pw_width <- 2 * z * es$se[es_fin]
	es_sim_width <- (es$ci_high - es$ci_low)[es_fin]
	expect_true(all(es_sim_width > es_pw_width))
})

# ------------------------------------------------------------------------------
# Test 2: ci_type = "pointwise" reproduces the EXACT pre-#197 bounds.
# THE backward-compat guard. Because getCohortATTsFinal() / eventStudy() are
# untouched on the pointwise path, the direct `estimate +/- z*se` formula IS
# the pre-#197 baseline -- so (a) doubles as the byte-for-byte lock (b).
# ------------------------------------------------------------------------------
test_that("ci_type = 'pointwise' reproduces the exact pre-#197 catt_df bounds", {
	sim <- make_ci_panel()
	fit <- fetwfeWithSimulatedData(sim, ci_type = "pointwise")
	expect_identical(fit$ci_type, "pointwise")

	z <- stats::qnorm(1 - fit$alpha / 2)
	cd <- fit$catt_df
	expect_equal(cd$ci_low, cd$estimate - z * cd$se, tolerance = 1e-12)
	expect_equal(cd$ci_high, cd$estimate + z * cd$se, tolerance = 1e-12)

	es <- eventStudy(fit, ci_type = "pointwise")
	expect_equal(es$ci_low, es$estimate - z * es$se, tolerance = 1e-12)
	expect_equal(es$ci_high, es$estimate + z * es$se, tolerance = 1e-12)
})

# ------------------------------------------------------------------------------
# Test 3: Simultaneous strictly wider than pointwise (the SUBSTANTIVE
# correctness-adjacent check; would FAIL on a ~50%-miscalibrated band). Gated
# on >= 2 non-degenerate cohorts (else the worker's K=1 bypass makes widths
# equal). Mirrors #192 Test 2's pointwise < simultaneous ordering.
# ------------------------------------------------------------------------------
test_that("simultaneous bands are strictly wider than pointwise (>= 2 cohorts)", {
	sim <- make_ci_panel(seed = 202, R = 5, T = 6, eff_size = 3)
	fit_s <- fetwfeWithSimulatedData(sim)
	fit_p <- fetwfeWithSimulatedData(sim, ci_type = "pointwise")

	cd_s <- fit_s$catt_df
	cd_p <- fit_p$catt_df
	fin <- is.finite(cd_s$se) & cd_s$se > 0
	expect_gte(sum(fin), 2L)
	w_s <- (cd_s$ci_high - cd_s$ci_low)[fin]
	w_p <- (cd_p$ci_high - cd_p$ci_low)[fin]
	expect_true(all(w_s > w_p))

	es_s <- eventStudy(fit_s)
	es_p <- eventStudy(fit_p)
	es_fin <- is.finite(es_s$se) & es_s$se > 0
	expect_gte(sum(es_fin), 2L)
	ew_s <- (es_s$ci_high - es_s$ci_low)[es_fin]
	ew_p <- (es_p$ci_high - es_p$ci_low)[es_fin]
	expect_true(all(ew_s > ew_p))
})

# ------------------------------------------------------------------------------
# Test 4: se / p_value / selected / estimate / overall-ATT CI are identical
# across ci_type. Bounds the blast radius to the two interval-bound columns.
# ------------------------------------------------------------------------------
test_that("only ci_low/ci_high differ across ci_type; se/p_value/selected/ATT identical", {
	sim <- make_ci_panel()
	fit_s <- fetwfeWithSimulatedData(sim)
	fit_p <- fetwfeWithSimulatedData(sim, ci_type = "pointwise")

	expect_equal(fit_s$catt_df$se, fit_p$catt_df$se)
	expect_equal(fit_s$catt_df$p_value, fit_p$catt_df$p_value)
	expect_identical(fit_s$catt_df$selected, fit_p$catt_df$selected)
	expect_equal(fit_s$catt_df$estimate, fit_p$catt_df$estimate)

	expect_equal(fit_s$att_hat, fit_p$att_hat)
	expect_equal(fit_s$att_se, fit_p$att_se)
	expect_equal(fit_s$att_p_value, fit_p$att_p_value)

	# Event-study se identical (only bounds differ).
	es_s <- eventStudy(fit_s)
	es_p <- eventStudy(fit_p)
	expect_equal(es_s$se, es_p$se)
	expect_equal(es_s$estimate, es_p$estimate)
})

# ------------------------------------------------------------------------------
# Test 5: ROW-ALIGNMENT + WIRING LOCK (NOT a correctness check). The fit-time
# catt_df bounds route to the #192 cohort worker and are column/row-aligned.
# Both sides route through the SAME .simultaneous_cis_impl() worker with the
# same alpha, so this is worker == worker -- it locks plumbing/alignment, not
# numerical correctness (which rests on test-simultaneous-cis.R's coverage test
# + Test 3's wider-than-pointwise check).
# ------------------------------------------------------------------------------
test_that("catt_df fit-time bounds are wired to the #192 cohort worker (row-aligned)", {
	sim <- make_ci_panel()
	fit <- fetwfeWithSimulatedData(sim)
	sci <- simultaneousCIs(fit, family = "cohort")
	expect_equal(fit$catt_df$ci_low, sci$ci$simultaneous_ci_low)
	expect_equal(fit$catt_df$ci_high, sci$ci$simultaneous_ci_high)
	expect_equal(nrow(sci$ci), nrow(fit$catt_df))
})

# ------------------------------------------------------------------------------
# Test 6: ROW-ALIGNMENT + WIRING LOCK (event-study analog of Test 5). Compared
# on finite-SE rows only -- the degenerate rows are NA in eventStudy()
# post-blocker-fix but 0 in the raw worker $ci (intentional; see Test 11).
# ------------------------------------------------------------------------------
test_that("eventStudy fit-time bounds are wired to the #192 event-study worker (row-aligned)", {
	sim <- make_ci_panel()
	fit <- fetwfeWithSimulatedData(sim)
	es <- eventStudy(fit)
	sci <- simultaneousCIs(fit, family = "event_study")
	fin <- is.finite(es$se)
	expect_equal(es$ci_low[fin], sci$ci$simultaneous_ci_low[fin])
	expect_equal(es$ci_high[fin], sci$ci$simultaneous_ci_high[fin])
})

# ------------------------------------------------------------------------------
# Test 7: calc_ses = FALSE degrades gracefully (no stop() at fit time). THE
# critical interaction-matrix cell. q = 1 forces calc_ses = FALSE.
# ------------------------------------------------------------------------------
test_that("ci_type = 'simultaneous' degrades gracefully when calc_ses = FALSE (q = 1)", {
	sim <- make_ci_panel()
	fit <- expect_no_error(fetwfeWithSimulatedData(sim, q = 1))
	expect_identical(fit$ci_type, "simultaneous")
	expect_false(isTRUE(fit$calc_ses))
	# Pointwise NA path is unchanged: bounds are NA (not 0, not widened).
	expect_true(all(is.na(fit$catt_df$ci_low)))
	expect_true(all(is.na(fit$catt_df$ci_high)))
	# eventStudy succeeds and returns NA CIs.
	es <- expect_no_error(eventStudy(fit))
	expect_true(all(is.na(es$ci_low)))
	expect_true(all(is.na(es$ci_high)))
})

# ------------------------------------------------------------------------------
# Test 8: K = 1 bypass via a DIRECT unit test of .apply_simultaneous_catt_band()
# on a hand-built 1-row catt_df. genCoefs(G = 1) errors (R >= 2 required) and
# genCoefs(T = 2) errors (T >= 3 required), so a K = 1 family is NOT
# simulator-constructible -- the direct unit test is the only viable route.
# The worker's K = 1 branch (R/simultaneous_cis.R: sum(nondeg) <= 1 ->
# crit = pointwise_crit) makes simultaneous == pointwise.
# ------------------------------------------------------------------------------
test_that("K = 1 cohort family: simultaneous bound equals pointwise (direct unit test)", {
	# Build a real 1-cohort etwfe-classed object is impossible via the
	# simulator; instead fit a normal etwfe and surgically reduce its cohort
	# family to a single cohort, then call the helper directly. The worker
	# reconstructs Sigma from the (single) selected treatment cell.
	sim <- make_ci_panel(R = 4, T = 6, eff_size = 2)
	fit <- etwfeWithSimulatedData(sim)
	# Sanity: the full helper widens for K >= 2.
	band_full <- fetwfe:::.apply_simultaneous_catt_band(
		fit,
		alpha = fit$alpha,
		has_valid_ses = fit$calc_ses
	)
	expect_false(is.null(band_full))

	# Direct K = 1 reach: a custom-family contrast of a single cohort routes
	# through the worker's K = 1 branch, returning crit = pointwise_crit.
	one_cohort <- matrix(0, nrow = 1L, ncol = length(fit$treat_inds))
	one_cohort[1, 1] <- 1
	sci1 <- fetwfe:::.simultaneous_cis_impl(
		x = fit,
		family = "custom",
		alpha = fit$alpha,
		contrasts = one_cohort,
		has_valid_ses = TRUE
	)
	expect_identical(sci1$K, 1L)
	# K = 1 -> simultaneous critical value == pointwise critical value.
	expect_equal(sci1$critical_value, sci1$pointwise_critical_value)
	# Hence the simultaneous and pointwise bounds coincide.
	expect_equal(sci1$ci$simultaneous_ci_low, sci1$ci$pointwise_ci_low)
	expect_equal(sci1$ci$simultaneous_ci_high, sci1$ci$pointwise_ci_high)
})

# ------------------------------------------------------------------------------
# Test 9: Validators accept both ci_type values; C9 catches a malformed
# ci_type; C10 catches a hand-narrowed simultaneous band.
# ------------------------------------------------------------------------------
test_that("validators accept both ci_type values; C9 / C10 catch violations", {
	sim <- make_ci_panel()
	fit_s <- fetwfeWithSimulatedData(sim)
	fit_p <- fetwfeWithSimulatedData(sim, ci_type = "pointwise")
	expect_silent(fetwfe:::.validate_fetwfe(fit_s))
	expect_silent(fetwfe:::.validate_fetwfe(fit_p))

	# C9: a bogus ci_type value.
	bad9 <- fit_s
	bad9$ci_type <- "bogus"
	expect_error(fetwfe:::.validate_fetwfe(bad9), "C9")

	# C10: a hand-narrowed simultaneous band (below pointwise width). Only
	# meaningful on a fit with >= 2 finite-SE cohorts.
	cd <- fit_s$catt_df
	fin <- is.finite(cd$se) & cd$se > 0
	skip_if_not(sum(fin) >= 2L, "fixture lacks >= 2 non-degenerate cohorts")
	bad10 <- fit_s
	cd10 <- bad10$catt_df
	# Narrow every interval to half a pointwise SE -- guaranteed < pointwise.
	cd10$ci_low <- cd10$estimate - 0.5 * cd10$se
	cd10$ci_high <- cd10$estimate + 0.5 * cd10$se
	bad10$catt_df <- cd10
	expect_error(fetwfe:::.validate_fetwfe(bad10), "C10")
})

# ------------------------------------------------------------------------------
# Test 10: Cross-estimator coverage (etwfe / betwfe / twfeCovs). Cohort family
# only for twfeCovs (no event-study surface). For betwfe include a
# ci_type x calc_ses = FALSE (q = 1) cell.
# ------------------------------------------------------------------------------
test_that("ci_type works across etwfe / betwfe / twfeCovs", {
	sim <- make_ci_panel()
	z <- stats::qnorm(1 - 0.05 / 2)

	for (est in c("etwfe", "betwfe", "twfeCovs")) {
		fn <- get(paste0(est, "WithSimulatedData"))
		fit_s <- fn(sim)
		fit_p <- fn(sim, ci_type = "pointwise")
		expect_identical(fit_s$ci_type, "simultaneous")
		expect_identical(fit_p$ci_type, "pointwise")
		# Pointwise reproduces the exact formula.
		cd_p <- fit_p$catt_df
		expect_equal(
			cd_p$ci_low,
			cd_p$estimate - z * cd_p$se,
			tolerance = 1e-12
		)
		# se identical across ci_type.
		expect_equal(fit_s$catt_df$se, fit_p$catt_df$se)
		# Alignment to the worker.
		sci <- simultaneousCIs(fit_s, family = "cohort")
		expect_equal(fit_s$catt_df$ci_low, sci$ci$simultaneous_ci_low)
	}

	# betwfe x calc_ses = FALSE (q = 1) cell: graceful degradation.
	fit_b1 <- expect_no_error(betwfeWithSimulatedData(sim, q = 1))
	expect_identical(fit_b1$ci_type, "simultaneous")
	expect_true(all(is.na(fit_b1$catt_df$ci_low)))
})

# ------------------------------------------------------------------------------
# Test 11: BLOCKER regression -- a degenerate event-time row stays NA, not 0,
# under the default simultaneous band. The worker returns 0/0 for the
# zero-variance effect; eventStudy() represents it as NA/NA; the fix preserves
# NA. WITHOUT the fix the worker's 0/0 would overwrite NA/NA.
# ------------------------------------------------------------------------------
test_that("degenerate event-time row stays NA (not 0) under the simultaneous default", {
	pdata <- make_scattered_degenerate_panel()
	fit <- etwfe(
		pdata = pdata,
		time_var = "year",
		unit_var = "unit",
		treatment = "treated",
		response = "y",
		verbose = FALSE
	)
	expect_identical(fit$ci_type, "simultaneous")
	es <- eventStudy(fit) # default simultaneous
	# There must be at least one degenerate (NA-se) event time for this test
	# to be non-vacuous.
	na_rows <- which(is.na(es$se))
	expect_gte(length(na_rows), 1L)
	# The blocker: those rows' bounds must be NA, NOT 0.
	expect_true(all(is.na(es$ci_low[na_rows])))
	expect_true(all(is.na(es$ci_high[na_rows])))
	# And finite-SE rows still carry (wider-than-pointwise) simultaneous bounds.
	fin <- is.finite(es$se) & es$se > 0
	expect_gte(sum(fin), 1L)
	expect_true(all(is.finite(es$ci_low[fin])))
})

# ------------------------------------------------------------------------------
# Test 12: OPTION B -- the broom tidiers reflect ci_type by passing through the
# object's ci_low/ci_high. Cross-surface consistency lock (tidy output vs the
# object's stored bounds -- independent surfaces, not a circular recompute).
# ------------------------------------------------------------------------------
test_that("broom tidiers reflect ci_type (Option B pass-through)", {
	skip_if_not_installed("broom")
	sim <- make_ci_panel()
	z <- stats::qnorm(1 - 0.05 / 2)

	fit_s <- fetwfeWithSimulatedData(sim)
	td_s <- broom::tidy(fit_s)
	cohort_rows <- td_s[td_s$term != "ATT", ]
	catt_sorted <- fit_s$catt_df[
		order(
			suppressWarnings(as.numeric(fit_s$catt_df$cohort)),
			fit_s$catt_df$cohort
		),
		,
		drop = FALSE
	]
	# (a) cohort rows pass through catt_df (simultaneous), same cohort order.
	expect_equal(cohort_rows$conf.low, catt_sorted$ci_low)
	expect_equal(cohort_rows$conf.high, catt_sorted$ci_high)
	# (b) overall-ATT row (row 1, K = 1) unchanged: scalar pointwise.
	expect_equal(td_s$conf.low[1], fit_s$att_hat - z * fit_s$att_se)
	expect_equal(td_s$conf.high[1], fit_s$att_hat + z * fit_s$att_se)
	# (c) tidy.eventStudy passes through eventStudy()$ci_low/high (finite rows).
	es_s <- eventStudy(fit_s)
	tes_s <- broom::tidy(es_s)
	fin <- is.finite(es_s$se)
	expect_equal(tes_s$conf.low[fin], es_s$ci_low[fin])
	expect_equal(tes_s$conf.high[fin], es_s$ci_high[fin])

	# (d) REVERT: ci_type = "pointwise" makes the cohort + event-study tidy
	# rows pointwise again.
	fit_p <- fetwfeWithSimulatedData(sim, ci_type = "pointwise")
	td_p <- broom::tidy(fit_p)
	cohort_rows_p <- td_p[td_p$term != "ATT", ]
	catt_p_sorted <- fit_p$catt_df[
		order(
			suppressWarnings(as.numeric(fit_p$catt_df$cohort)),
			fit_p$catt_df$cohort
		),
		,
		drop = FALSE
	]
	expect_equal(
		cohort_rows_p$conf.low,
		catt_p_sorted$estimate - z * catt_p_sorted$se,
		tolerance = 1e-12
	)
	es_p <- eventStudy(fit_p)
	tes_p <- broom::tidy(es_p)
	fin_p <- is.finite(es_p$se)
	expect_equal(
		tes_p$conf.low[fin_p],
		es_p$estimate[fin_p] - z * es_p$se[fin_p],
		tolerance = 1e-12
	)

	# (e) cross-surface: tidy.cohortStudy agrees with tidy cohort rows under
	# both ci_type values.
	cs_s <- cohortStudy(fit_s)
	tcs_s <- broom::tidy(cs_s)
	expect_equal(unname(tcs_s$conf.low), unname(cohort_rows$conf.low))
	cs_p <- cohortStudy(fit_p)
	tcs_p <- broom::tidy(cs_p)
	expect_equal(unname(tcs_p$conf.low), unname(cohort_rows_p$conf.low))
})

# ------------------------------------------------------------------------------
# Test 13: alpha interaction -- a non-default alpha flows into the band, and the
# conservative se_type does NOT emit chatter at fit time (suppressMessages in
# both fit-time helpers; the message is reachable only via a DIRECT call, not
# *WithSimulatedData() which forces the tight path via indep_counts).
# ------------------------------------------------------------------------------
test_that("ci_type honors a non-default alpha and suppresses the conservative message", {
	sim <- make_ci_panel()
	fit_90 <- fetwfeWithSimulatedData(sim, alpha = 0.1)
	expect_equal(fit_90$alpha, 0.1)
	# 90% simultaneous band: wider than 90% pointwise, narrower than the 95%
	# simultaneous band (more confidence -> wider).
	cd90 <- fit_90$catt_df
	fin <- is.finite(cd90$se) & cd90$se > 0
	skip_if_not(sum(fin) >= 2L)
	z90 <- stats::qnorm(1 - 0.1 / 2)
	pw90 <- 2 * z90 * cd90$se[fin]
	sim90 <- (cd90$ci_high - cd90$ci_low)[fin]
	expect_true(all(sim90 > pw90))

	fit_95 <- fetwfeWithSimulatedData(sim) # alpha = 0.05 default
	sim95 <- (fit_95$catt_df$ci_high - fit_95$catt_df$ci_low)[fin]
	expect_true(all(sim95 > sim90)) # 95% wider than 90%

	# Conservative se_type via a DIRECT call (not *WithSimulatedData, which
	# passes indep_counts and forces the tight path): no unprompted message().
	skip_if_not_installed("bacondecomp")
	# Build a small direct-call panel from the simulated object.
	pdata <- sim$pdata
	expect_silent(
		etwfe(
			pdata = pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			response = sim$response,
			covs = sim$covs,
			sig_eps_sq = sim$sig_eps_sq,
			sig_eps_c_sq = sim$sig_eps_c_sq,
			se_type = "conservative",
			verbose = FALSE
		)
	)
})

# ==============================================================================
# New-option interaction matrix (#197 post-review R1; PLAN.md Validation
# "New-option interaction matrix"). ci_type (default "simultaneous", plus a
# "pointwise" revert check) crossed with each fit-time option that can
# interact: q in {1, 2}, add_ridge, se_type = "cluster", indep_counts,
# allow_no_never_treated. The q = 0.5 default + alpha are already covered by
# Tests 1-13. Each cell asserts (per expect_simultaneous_matrix_cell):
# (1) the fit runs under the simultaneous default; (2) ci_type slot set +
# catt_df well-formed; (3) >= 2 finite-SE cohorts -> wider-or-equal vs
# pointwise; (4) a pointwise refit reverts. Behavioral locks (widen + revert),
# not well-formedness-only. Seeds pinned; the cluster cross-estimator cell
# (most refits) is skip_on_cran-gated.
# ==============================================================================

# ------------------------------------------------------------------------------
# Test 14: ci_type x q in {1, 2} (fetwfe + betwfe; lasso q = 1 and ridge q = 2).
# Both q >= 1 force calc_ses = FALSE (the oracle-property gate; SEs require
# q < 1 -- cf. test-fetwfe.R "q >= 1 returns NA SE"), so the simultaneous path
# degrades to NA bounds with NO stop() at fit time. NOTE: the plan's
# parenthetical "q = 2 -> all cohorts selected so the band spans all R cohorts"
# refers to ridge SELECTION; the BAND still degrades to NA because calc_ses is
# FALSE for any q >= 1. The substantive q-interaction lock is therefore the
# graceful-degradation NA path (the wider-band q-cell is q = 0.5, Tests 1/3).
# ------------------------------------------------------------------------------
test_that("ci_type x q in {1, 2}: q >= 1 (calc_ses FALSE) degrades to NA bounds (fetwfe/betwfe)", {
	sim <- make_ci_panel()
	for (qv in c(1, 2)) {
		for (est in c("fetwfe", "betwfe")) {
			fn <- get(paste0(est, "WithSimulatedData"))
			fit <- expect_no_error(fn(sim, q = qv))
			expect_identical(fit$ci_type, "simultaneous")
			expect_false(isTRUE(fit$calc_ses))
			# Pointwise NA path is unchanged: bounds NA (not 0, not widened).
			expect_true(all(is.na(fit$catt_df$ci_low)))
			expect_true(all(is.na(fit$catt_df$ci_high)))
			# Validator accepts the degraded simultaneous object.
			validator <- get(
				paste0(".validate_", est),
				envir = asNamespace("fetwfe")
			)
			expect_silent(validator(fit))
		}
	}
	# fetwfe q >= 1 has an event-study surface; it also degrades to NA.
	es <- expect_no_error(eventStudy(fetwfeWithSimulatedData(sim, q = 1)))
	expect_true(all(is.na(es$ci_low)))
})

# ------------------------------------------------------------------------------
# Test 15: ci_type x add_ridge = TRUE (fetwfe + etwfe + twfeCovs). add_ridge
# perturbs the coefficients pre-fit; it is orthogonal to the band computation
# (which reads se / Sigma off the fitted object). Smoke: the band still widens
# and reverts.
# ------------------------------------------------------------------------------
test_that("ci_type x add_ridge = TRUE: band widens and reverts (fetwfe/etwfe/twfeCovs)", {
	sim <- make_ci_panel()
	for (est in c("fetwfe", "etwfe", "twfeCovs")) {
		fn <- get(paste0(est, "WithSimulatedData"))
		fit_s <- expect_no_error(fn(sim, add_ridge = TRUE))
		fit_p <- fn(sim, add_ridge = TRUE, ci_type = "pointwise")
		expect_simultaneous_matrix_cell(fit_s, fit_p, label = est)
	}
})

# ------------------------------------------------------------------------------
# Test 16: ci_type x se_type = "cluster" -- THE genuine interaction. The
# cluster-robust sandwich changes Sigma_1, so the simultaneous band must
# consume the cluster covariance (not the default one). Run on an AR(1) panel
# so the sandwich bites AND >= 2 finite-SE cohorts survive (strict widening
# assertable). Cohort family for all four estimators; event-study widening
# checked for the estimators that expose it. skip_on_cran (cross-estimator
# refits).
# ------------------------------------------------------------------------------
test_that("ci_type x se_type = 'cluster': simultaneous band consumes the cluster cov (all estimators)", {
	skip_on_cran()
	mk <- make_ci_ar1_panel()
	z <- stats::qnorm(1 - 0.05 / 2)

	for (est in c("fetwfe", "etwfe", "betwfe", "twfeCovs")) {
		fn <- get(paste0(est, "WithSimulatedData"))
		fit_s <- expect_no_error(fn(mk, se_type = "cluster"))
		fit_p <- fn(mk, se_type = "cluster", ci_type = "pointwise")
		expect_identical(fit_s$se_type, "cluster")
		expect_simultaneous_matrix_cell(fit_s, fit_p, label = est)

		# The cluster panel keeps >= 2 finite-SE cohorts, so demand STRICT
		# widening here (the genuine-interaction substance, not >=).
		cd_s <- fit_s$catt_df
		cd_p <- fit_p$catt_df
		fin <- is.finite(cd_s$se) & cd_s$se > 0
		expect_gte(sum(fin), 2L)
		w_s <- (cd_s$ci_high - cd_s$ci_low)[fin]
		w_p <- (cd_p$ci_high - cd_p$ci_low)[fin]
		expect_true(all(w_s > w_p))

		# The band is wired to the cluster-aware worker (row-aligned).
		sci <- simultaneousCIs(fit_s, family = "cohort")
		expect_equal(cd_s$ci_low, sci$ci$simultaneous_ci_low)
	}

	# Event-study surface (fetwfe/etwfe/betwfe; twfeCovs has none) also widens
	# under the cluster covariance.
	for (est in c("fetwfe", "etwfe", "betwfe")) {
		fn <- get(paste0(est, "WithSimulatedData"))
		fit_s <- fn(mk, se_type = "cluster")
		fit_p <- fn(mk, se_type = "cluster", ci_type = "pointwise")
		es_s <- eventStudy(fit_s)
		es_p <- eventStudy(fit_p)
		es_fin <- is.finite(es_s$se) & es_s$se > 0
		if (sum(es_fin) >= 2L) {
			ew_s <- (es_s$ci_high - es_s$ci_low)[es_fin]
			ew_p <- (es_p$ci_high - es_p$ci_low)[es_fin]
			expect_true(all(ew_s > ew_p))
		}
	}
})

# ------------------------------------------------------------------------------
# Test 17: ci_type x indep_counts supplied. indep_counts makes the OVERALL-ATT
# two-sample exact; the COHORT / EVENT-STUDY families' Sigma is unaffected by
# it (it governs only the overall-ATT combination), so the simultaneous
# cohort/event-study band is well-formed and widens regardless. The
# *WithSimulatedData wrappers always pass the simulator's indep_counts (the
# tight path), so this cell locks that the band is correct WHEN indep_counts
# was consumed (indep_counts_used == TRUE).
# ------------------------------------------------------------------------------
test_that("ci_type x indep_counts: simultaneous band well-formed when indep_counts consumed", {
	sim <- make_ci_panel()
	expect_false(is.null(sim$indep_counts)) # fixture supplies them

	for (est in c("fetwfe", "etwfe", "betwfe", "twfeCovs")) {
		fn <- get(paste0(est, "WithSimulatedData"))
		fit_s <- expect_no_error(fn(sim))
		fit_p <- fn(sim, ci_type = "pointwise")
		# indep_counts was actually consumed (the two-sample overall-ATT path).
		expect_true(isTRUE(fit_s$indep_counts_used))
		expect_simultaneous_matrix_cell(fit_s, fit_p, label = est)
	}
})

# ------------------------------------------------------------------------------
# Test 18: ci_type x allow_no_never_treated = TRUE. A panel with NO
# never-treated units is auto-truncated (a panel-preprocessing toggle,
# orthogonal to the band). The truncated panel still produces a catt_df that
# flows through the finalizer. On this fixture selection zeroes every cohort
# (the .build_selected_out_result() degenerate fallback: estimate = se = 0,
# ci_low = ci_high = 0), so there are 0 finite-SE cohorts -- the wider-than-
# pointwise assertion is correctly gated off (no >= 2 finite-SE cohorts); the
# lock here is "runs + ci_type set + well-formed degenerate (0) bounds +
# pointwise revert + validates". Uses generate_bad_panel_data + a DIRECT
# fetwfe() call (the bad-panel helper returns a bare data.frame, and the
# truncation emits a one-time "auto-truncated" warning).
# ------------------------------------------------------------------------------
test_that("ci_type x allow_no_never_treated = TRUE: finalizer handles the auto-truncated degenerate fit", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)
	common <- list(
		pdata = df_bad,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		allow_no_never_treated = TRUE,
		verbose = FALSE
	)
	expect_warning(
		fit_s <- do.call(fetwfe, common),
		"auto-truncated"
	)
	expect_warning(
		fit_p <- do.call(fetwfe, c(common, list(ci_type = "pointwise"))),
		"auto-truncated"
	)
	# Cell assertions (gated: 0 finite-SE cohorts here, so no strict widening).
	expect_simultaneous_matrix_cell(
		fit_s,
		fit_p,
		label = "allow_no_never_treated"
	)
	# Explicitly confirm the degenerate-fallback convention on this fixture.
	expect_false(any(is.finite(fit_s$catt_df$se) & fit_s$catt_df$se > 0))
	expect_true(all(fit_s$catt_df$ci_low == 0))
	expect_true(all(fit_s$catt_df$ci_high == 0))
	# Validator accepts the simultaneous (degenerate) object.
	expect_silent(fetwfe:::.validate_fetwfe(fit_s))
})

# ------------------------------------------------------------------------------
# Test 19: REGRESSION (#204) -- twfeCovs() honors a non-default `alpha` when
# building the default simultaneous catt_df band. Before the fix, twfeCovs had
# no `alpha` slot and the entry-point hardcoded `.finalize_ci_type(out, alpha =
# 0.05)`, so the simultaneous band was ALWAYS the 0.05 band regardless of the
# requested level (and C10's `.alpha_of` returned the same wrong 0.05, so the
# band-width validator was blind). This test FAILS on the old code: it asserts
# (i) the new `alpha` slot is set; (ii) the catt_df bounds match
# `simultaneousCIs(fit, "cohort", alpha = 0.10)` (not the 0.05 band); (iii) the
# 0.10 band is strictly NARROWER than the 0.05 band on >= 2 finite-SE cohorts;
# (iv) `.validate_twfeCovs(fit)` is silent. Both the *WithSimulatedData wrapper
# and a DIRECT twfeCovs() call are exercised.
# ------------------------------------------------------------------------------
test_that("twfeCovs() honors a non-default alpha for the simultaneous band (#204)", {
	sim <- make_ci_panel()

	check_alpha_honored <- function(fit, label) {
		# (i) the new alpha slot carries the requested level.
		expect_equal(fit$alpha, 0.1, info = label)
		# (ii) catt_df bounds are the alpha = 0.10 simultaneous band, NOT 0.05.
		sci10 <- simultaneousCIs(fit, family = "cohort", alpha = 0.1)
		expect_equal(
			fit$catt_df$ci_low,
			sci10$ci$simultaneous_ci_low,
			info = label
		)
		expect_equal(
			fit$catt_df$ci_high,
			sci10$ci$simultaneous_ci_high,
			info = label
		)
		# (iii) the 0.10 band is strictly narrower than the 0.05 band (the bug:
		# it was wrongly EQUAL to the 0.05 band). Gated on >= 2 finite-SE cohorts.
		cd <- fit$catt_df
		fin <- is.finite(cd$se) & cd$se > 0
		skip_if_not(sum(fin) >= 2L)
		sci05 <- simultaneousCIs(fit, family = "cohort", alpha = 0.05)
		w10 <- (cd$ci_high - cd$ci_low)[fin]
		w05 <- (sci05$ci$simultaneous_ci_high -
			sci05$ci$simultaneous_ci_low)[fin]
		expect_true(all(w10 < w05), info = label)
		# (iv) the constructed object validates silently with alpha = 0.10.
		expect_silent(fetwfe:::.validate_twfeCovs(fit))
	}

	# Wrapper path.
	fit_w <- twfeCovsWithSimulatedData(sim, alpha = 0.1)
	check_alpha_honored(fit_w, "twfeCovsWithSimulatedData")

	# Direct twfeCovs() path (mirrors Test 13's direct-call construction).
	fit_d <- twfeCovs(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		alpha = 0.1,
		verbose = FALSE
	)
	check_alpha_honored(fit_d, "twfeCovs")
})
