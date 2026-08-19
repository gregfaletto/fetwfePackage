# Issue #433: a `p >= N*T` call that returns the fixed-`p` POST-SELECTION band
# rather than a uniformly-valid one now says so.
#
# One route out of four is valid in that regime -- a `fetwfe()` fit under
# `method = "bootstrap"`, the desparsified `debiased.highdim.joint.thm` band.
# Every other combination falls back to the post-selection selected-support
# band, which a #308 coverage study measures at ~0.09 against 0.95 nominal. The
# predicate at the shared step-7b site in `.simultaneous_cis_impl()` is exactly
# the complement of that one valid route:
#
#     p >= N * T_ && !(is_fetwfe && identical(method, "bootstrap"))
#                 && isTRUE(warn_highdim_postselection)
#
# so it fires for a non-`fetwfe` fit on EITHER method (this folds in the old
# bootstrap-only #308 warning) and for a `fetwfe` fit on the analytic default.
#
# THE ASSERTIONS COME IN TWO KINDS, and mixing them up is how a red-green
# measurement gets taken over the wrong set:
#   * "(red)" blocks fail on main and pass here -- they demonstrate the change.
#   * "(guardrail)" blocks pass on main too -- they constrain the implementation
#     so a too-wide predicate, a lost carve-out, or a moved number goes red.
#
# EVERY text anchor in this file passes `fixed = TRUE`. Two of the natural
# anchors here are broken as regexes, both measured: `grepl("use a fetwfe() fit",
# msg)` is FALSE (the parentheses are an empty capture group) and
# `grepl("[pointwise 95% CI]", <a simultaneous label>)` is TRUE (square brackets
# are a character class). Anchors are also unique to the message they target --
# never the shared word `unreliable`, which the #304 warning also contains
# (issue #458).
#
# COUNTS ARE COUNTED, NOT MATCHED. This suite runs testthat edition 2
# (`DESCRIPTION` carries no `Config/testthat/edition`), under which
# `expect_warning(expr, pattern)` PASSES when a second, unmatched warning is
# also signalled -- measured. So "exactly one warning" uses
# `capture_warnings()` + `expect_length()`, and "exactly once in the rendered
# output" uses `sum(grepl(..., fixed = TRUE))`, never `expect_match()`.

# ------------------------------------------------------------------------------
# Fixtures. Built LAZILY and memoized, not at file scope: every block below
# calls skip_on_cran() first, and a file-scope fixture would build four
# high-dimensional bridge fits on CRAN before the first skip fires.
# ------------------------------------------------------------------------------
.hpb433_cache <- new.env(parent = emptyenv())

.hpb433_build <- function(name) {
	switch(
		name,
		# p = 356 against N*T = 300, calc_ses = TRUE, three treatment
		# coefficients of which one survives the bridge -- high-dimensional and
		# NOT degenerate, which is what separates this from the #304 path.
		sim = simulateData(
			genCoefs(
				G = 3,
				T = 5,
				d = 20,
				density = 0.08,
				eff_size = 6,
				seed = 11
			),
			N = 60,
			sig_eps_sq = 0.5,
			sig_eps_c_sq = 0.5,
			seed = 1001
		),
		hd_betwfe = betwfeWithSimulatedData(.hpb433("sim")),
		hd_fetwfe = fetwfeWithSimulatedData(.hpb433("sim")),
		hd_betwfe_pointwise = betwfeWithSimulatedData(
			.hpb433("sim"),
			ci_type = "pointwise"
		),
		# p = 50 against N*T = 600.
		lo_sim = simulateData(
			genCoefs(
				G = 3,
				T = 5,
				d = 2,
				density = 0.5,
				eff_size = 2,
				seed = 7
			),
			N = 120,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 7
		),
		lo_betwfe = betwfeWithSimulatedData(.hpb433("lo_sim")),
		# The degenerate high-dimensional fixture, taken verbatim from
		# tests/testthat/test-degenerate-band-policy-429.R so the two files pin
		# the same regime: p = 311 against N*T = 180, every catt_hat exactly 0.
		degen_sim = simulateData(
			genCoefs(
				G = 4,
				T = 6,
				d = 12,
				density = 0.5,
				eff_size = 2,
				seed = 55
			),
			N = 30,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 55
		),
		degen_betwfe = betwfeWithSimulatedData(.hpb433("degen_sim"), q = 0.5),
		# The `gls = FALSE` shape from
		# tests/testthat/test-print-highdim-diagnostics-326.R: p = 274 against
		# N*T = 192, ci_type = "simultaneous", calc_ses = FALSE, and `catt_df`
		# bounds all NA because the fit-time band came back NULL.
		gls_false_fetwfe = .hpb433_gls_false_fit(),
		stop("unknown fixture: ", name)
	)
}

.hpb433 <- function(name) {
	if (!exists(name, envir = .hpb433_cache, inherits = FALSE)) {
		assign(name, .hpb433_build(name), envir = .hpb433_cache)
	}
	get(name, envir = .hpb433_cache, inherits = FALSE)
}

.hpb433_gls_false_fit <- function() {
	set.seed(5)
	d <- 10L
	N <- 24L
	T_ <- 8L
	adopt <- c(2L, 4L, 7L)
	covs <- matrix(stats::rnorm(N * d), N, d)
	cou <- c(
		rep(0L, N - 18L),
		rep(adopt[1], 6),
		rep(adopt[2], 6),
		rep(adopt[3], 6)
	)
	eff <- stats::setNames(c(0.5, 2, 3.5), as.character(adopt))
	panel <- do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cou[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T_,
				treat = as.integer(g > 0 & (1:T_) >= g)
			)
			for (j in 1:d) {
				df[[paste0("x", j)]] <- covs[i, j]
			}
			df$y <- 0.3 *
				(1:T_) /
				T_ +
				(if (g > 0) eff[[as.character(g)]] else 0) * df$treat +
				0.2 * covs[i, 1] +
				stats::rnorm(T_, 0, 0.4)
			df
		})
	)
	fetwfe(
		pdata = panel,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:10),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
}

# Warning-message anchors, resolved once. `no desparsified band` is carried by
# BOTH branches -- on the `fetwfe` branch it describes the analytic route, not
# the estimator -- and it is the phrase that discriminates this warning from the
# #304 degenerate one.
.HPB433_ANCHOR <- "no desparsified band"
.HPB433_REMEDY_OTHER <- "use a fetwfe() fit"
.HPB433_REMEDY_FETWFE <- "method = \"bootstrap\""
.HPB433_NOTICE <- "post-selection fallback"
.HPB433_CLASS <- "fetwfe_highdim_postselection_band"

# ------------------------------------------------------------------------------
# The fixture invariants come FIRST, so every block below runs on an asserted
# regime. Without this, a fixture that drifted out of `p >= N*T`, or that became
# degenerate, would make the silence assertions pass for the wrong reason and
# the warning assertions fail for an unrelated one.
# ------------------------------------------------------------------------------
test_that("the #433 fixtures really are high-dimensional and non-degenerate", {
	skip_on_cran()

	for (nm in c("hd_betwfe", "hd_fetwfe")) {
		fit <- .hpb433(nm)
		expect_true(fit$p >= fit$N * fit$T, info = nm) # measured: 356 vs 300
		expect_true(fit$calc_ses, info = nm)
		expect_identical(fit$ci_type, "simultaneous", info = nm)
		# NOT degenerate: the #304 early return would otherwise pre-empt the
		# site this file is about, and every warning assertion below would be
		# reading the wrong condition.
		expect_false(all(fit$catt_hats == 0), info = nm)
	}

	lo <- .hpb433("lo_betwfe")
	expect_true(lo$p < lo$N * lo$T) # measured: 50 vs 600
	expect_true(lo$calc_ses)
	expect_identical(lo$ci_type, "simultaneous")

	dg <- .hpb433("degen_betwfe")
	expect_true(dg$p >= dg$N * dg$T) # measured: 311 vs 180
	expect_true(all(dg$catt_hats == 0))
	expect_true(dg$calc_ses)
})

# ------------------------------------------------------------------------------
# 1. (red) A non-fetwfe high-dimensional fit warns on the ANALYTIC DEFAULT --
#    the route the issue was filed about, and the one that was silent.
# ------------------------------------------------------------------------------
test_that("betwfe p >= NT warns on the analytic default, every family (#433)", {
	skip_on_cran()

	fit <- .hpb433("hd_betwfe")

	for (fam in c("event_study", "cohort", "all_post_treatment")) {
		ws <- testthat::capture_warnings(simultaneousCIs(fit, family = fam))
		# Exactly one: a second condition here would mean the degenerate
		# carve-out or some other guard had started firing alongside.
		expect_length(ws, 1L)
		expect_true(grepl(.HPB433_ANCHOR, ws, fixed = TRUE), info = fam)
		# The remedy is the estimator-level one: no valid high-dimensional band
		# exists for a non-fetwfe fit on ANY route, so `method = "bootstrap"`
		# would be wrong advice here.
		expect_true(grepl(.HPB433_REMEDY_OTHER, ws, fixed = TRUE), info = fam)
		expect_false(grepl(.HPB433_REMEDY_FETWFE, ws, fixed = TRUE), info = fam)
	}

	# The predicate is family-agnostic, so the fourth (user-supplied) family is
	# covered too. One contrast over the per-cohort effect vector.
	ctr <- matrix(0, nrow = 1L, ncol = length(fit$treat_inds))
	ctr[1, 1] <- 1
	ws_custom <- testthat::capture_warnings(
		simultaneousCIs(fit, family = "custom", contrasts = ctr)
	)
	expect_length(ws_custom, 1L)
	expect_true(grepl(.HPB433_ANCHOR, ws_custom, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 2. (red) A fetwfe fit warns on the analytic default too -- the half of the
#    scope the issue text does not cover, folded in deliberately. Here the
#    remedy is actionable: the same fit under method = "bootstrap" is valid.
# ------------------------------------------------------------------------------
test_that("fetwfe p >= NT warns on the analytic default, naming bootstrap (#433)", {
	skip_on_cran()

	fit <- .hpb433("hd_fetwfe")
	ws <- testthat::capture_warnings(simultaneousCIs(
		fit,
		family = "event_study"
	))

	expect_length(ws, 1L)
	expect_true(grepl(.HPB433_ANCHOR, ws, fixed = TRUE))
	expect_true(grepl(.HPB433_REMEDY_FETWFE, ws, fixed = TRUE))
	# NOT the "use a fetwfe() fit" remedy: this IS a fetwfe fit, and telling its
	# user to switch estimator would be nonsense.
	expect_false(grepl(.HPB433_REMEDY_OTHER, ws, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 3. (guardrail) The one valid route stays silent. This is the assertion that
#    stops a too-wide predicate: drop the
#    `!(is_fetwfe && identical(method, "bootstrap"))` conjunct and this block
#    goes red while every other warning assertion in the file still passes.
# ------------------------------------------------------------------------------
test_that("fetwfe p >= NT under method = 'bootstrap' does NOT warn (#433)", {
	skip_on_cran()

	fit <- .hpb433("hd_fetwfe")
	ws <- testthat::capture_warnings(
		bo <- simultaneousCIs(
			fit,
			family = "event_study",
			method = "bootstrap",
			B = 200,
			seed = 1
		)
	)
	expect_length(ws, 0L)
	# ... and it really did take the desparsified path, rather than staying
	# quiet because it silently fell back somewhere.
	expect_identical(bo$regime, "high-dimensional")
})

# ------------------------------------------------------------------------------
# 4. (guardrail) The behavior the fold-in had to preserve: the old #308
#    bootstrap-only warning still fires for a non-fetwfe fit, now from the
#    shared site. `regime` pins that it is still the fixed-p path.
# ------------------------------------------------------------------------------
test_that("betwfe p >= NT under method = 'bootstrap' still warns (#308 fold-in)", {
	skip_on_cran()

	fit <- .hpb433("hd_betwfe")
	ws <- testthat::capture_warnings(
		bo <- simultaneousCIs(
			fit,
			family = "event_study",
			method = "bootstrap",
			B = 200,
			seed = 1
		)
	)
	expect_length(ws, 1L)
	expect_true(grepl(.HPB433_ANCHOR, ws, fixed = TRUE))
	expect_true(grepl(.HPB433_REMEDY_OTHER, ws, fixed = TRUE))
	expect_identical(bo$regime, "fixed-p")
})

# ------------------------------------------------------------------------------
# 5. (guardrail) Fit time stays silent. The existing coverage in
#    test-degenerate-band-policy-429.R uses a DEGENERATE fixture, which returns
#    before the site this file is about, so it cannot observe this.
#
#    THE ANTI-VACUITY HALF MATTERS MORE THAN THE SILENCE HALF.
#    `.fit_band_for_family()` degrades to NULL through a `tryCatch(error = )`
#    and a row-count check, and a fit that never reached the worker is
#    trivially warning-free. So assert positively that the band WAS built and
#    applied. Note what is deliberately NOT asserted: that the applied bounds
#    differ from the pointwise ones. Measured on this fixture, they are equal to
#    the last digit -- one cohort survives the bridge, so K = 1 and the worker
#    takes the bypass where the simultaneous critical value IS the pointwise
#    one. A width comparison here would be a guardrail that is blind by
#    construction.
# ------------------------------------------------------------------------------
test_that("a fit-time simultaneous band on a p >= NT fit is silent (#433)", {
	skip_on_cran()

	# A DELIBERATELY SEPARATE fit: expect_no_warning() has to wrap the call
	# itself, so the fit-time half of the policy costs one extra fit. Do not
	# "optimize" this into the memoized fixture -- that deletes the assertion.
	expect_no_warning(fit <- betwfeWithSimulatedData(.hpb433("sim")))

	expect_true(fit$p >= fit$N * fit$T)
	expect_identical(fit$ci_type, "simultaneous")
	expect_true(all(!is.na(fit$catt_df$ci_low)))
	expect_true(all(!is.na(fit$catt_df$ci_high)))

	# The fit-time band path runs to completion on this fixture and its output
	# is what is stored, so the silence above is silence from the flag rather
	# than from an early exit.
	band <- fetwfe:::.apply_simultaneous_catt_band(
		fit,
		alpha = fit$alpha,
		has_valid_ses = TRUE
	)
	expect_false(is.null(band))
	expect_length(band$ci_low, nrow(fit$catt_df))
	expect_equal(fit$catt_df$ci_low, band$ci_low)
	expect_equal(fit$catt_df$ci_high, band$ci_high)
})

# ------------------------------------------------------------------------------
# 6. (red) eventStudy() called directly warns. It is the last route that handed
#    a user a post-selection band programmatically while saying nothing.
#
#    THE CLASS IS ASSERTED, NOT ONLY THE TEXT. The three internal callers muffle
#    on the class, so an unclassed warning would silently break their
#    suppression while every text assertion in this file still passed.
# ------------------------------------------------------------------------------
test_that("eventStudy() on a p >= NT fit warns, as a classed condition (#433)", {
	skip_on_cran()

	for (nm in c("hd_betwfe", "hd_fetwfe")) {
		fit <- .hpb433(nm)

		ws <- testthat::capture_warnings(eventStudy(fit))
		expect_length(ws, 1L)
		expect_true(grepl(.HPB433_ANCHOR, ws, fixed = TRUE), info = nm)

		cond <- tryCatch(eventStudy(fit), warning = function(w) w)
		expect_s3_class(cond, .HPB433_CLASS)
		expect_s3_class(cond, "warning")

		# The muffle does not swallow the computation: the quiet path returns
		# the same data frame, band rows included.
		expect_equal(
			fetwfe:::.event_study_quiet(fit),
			suppressWarnings(eventStudy(fit)),
			info = nm
		)
	}
})

# ------------------------------------------------------------------------------
# 6b. (guardrail) The renderers muffle it. Asserting the silence ALONE would
#     pass if the warning were deleted outright, so the same block asserts that
#     eventStudy() -- which all three of them call internally -- is loud on the
#     very same fit.
# ------------------------------------------------------------------------------
test_that("print/summary/plot muffle the #433 condition while eventStudy warns", {
	skip_on_cran()

	fit <- .hpb433("hd_betwfe")

	# The pair. If this ever reads 0, the block below is vacuous.
	expect_length(testthat::capture_warnings(eventStudy(fit)), 1L)

	expect_length(testthat::capture_warnings(capture.output(print(fit))), 0L)
	expect_length(testthat::capture_warnings(summary(fit)), 0L)
	expect_length(
		testthat::capture_warnings(capture.output(print(summary(fit)))),
		0L
	)
	skip_if_not_installed("ggplot2")
	expect_length(
		testthat::capture_warnings(
			ggplot2::ggplot_build(plot(fit, type = "event_study"))
		),
		0L
	)
})

# ------------------------------------------------------------------------------
# 7. (guardrail) A low-dimensional fit warns on neither route. Drop the
#    `p >= N * T_` conjunct and this goes red.
# ------------------------------------------------------------------------------
test_that("a p < NT fit warns on neither route (#433)", {
	skip_on_cran()

	fit <- .hpb433("lo_betwfe")

	expect_length(
		testthat::capture_warnings(simultaneousCIs(
			fit,
			family = "event_study"
		)),
		0L
	)
	expect_length(
		testthat::capture_warnings(
			simultaneousCIs(
				fit,
				family = "event_study",
				method = "bootstrap",
				B = 200,
				seed = 1
			)
		),
		0L
	)
	expect_length(testthat::capture_warnings(eventStudy(fit)), 0L)
})

# ------------------------------------------------------------------------------
# 8. (guardrail) The degenerate carve-out. On an all-zero-support high-
#    dimensional fit EXACTLY ONE warning fires and it is the #304 one -- the new
#    site sits after that path's early return.
#
#    NOT written with expect_warning(): under edition 2 that passes with a
#    second warning present, so it would ship green under the exact mutation it
#    exists to catch (hoisting the new block above the degenerate early return).
# ------------------------------------------------------------------------------
test_that("a degenerate p >= NT fit keeps only its #304 warning (#433)", {
	skip_on_cran()

	fit <- .hpb433("degen_betwfe")

	for (fam in c("cohort", "event_study")) {
		ws <- testthat::capture_warnings(simultaneousCIs(fit, family = fam))
		expect_length(ws, 1L)
		expect_true(
			grepl("bridge selection is NOT ", ws, fixed = TRUE),
			info = fam
		)
		expect_false(grepl(.HPB433_ANCHOR, ws, fixed = TRUE), info = fam)
	}
})

# ------------------------------------------------------------------------------
# 9. (guardrail, numeric) NO NUMBER MOVES. The band is pinned from the base tree
#    at 5a69eeb, measured to 15 significant digits and asserted to 1e-8.
#
#    NOT `identical()` and NOT `tolerance = 0`: these are qmvnorm / Gram-inverse
#    quantities and the gate runs on six platform/BLAS combinations, four of
#    them OpenBLAS. `.workflow/PROFILE.md` section 12.8 records the exact
#    failure -- assertions that passed on Accelerate and Windows and failed all
#    four Linux jobs at 1-4 ULPs.
#
#    WHAT 1e-8 ACTUALLY ACCEPTS HERE, measured rather than reasoned about (the
#    edition-2 comparison switches between relative and absolute with the
#    magnitude of the expected value, so "1e-8" names no unit on its own):
#    perturbing ONE element passes at 5e-9 and fails at 1e-8, i.e. an effective
#    absolute band of ~1e-8 on values of order 1-5. `all.equal.numeric` defaults
#    to `countEQ = FALSE`, so the mean is taken over only the differing entries
#    and a single bad element is not diluted by the other three. That is ~7
#    orders of magnitude above the 1-4 ULP cross-BLAS noise section 12.8
#    records, and far below anything a real change produces: a 0.1% shift in one
#    bound fails, and the exact-zero row moving to 1e-6 fails.
#
#    The event_study family, not cohort: on this fixture the cohort band is
#    (3.31769, 0, 0) with critical_value = qnorm(0.975) via the K = 1 bypass, so
#    a cohort pin is very nearly vacuous. Here K = 4 with three live intervals
#    and a genuinely simultaneous critical value (2.3486 against a pointwise
#    1.9600), so the pin covers the critical-value machinery too.
# ------------------------------------------------------------------------------
test_that("the p >= NT band's numbers are unchanged by #433", {
	skip_on_cran()

	sc <- suppressWarnings(
		simultaneousCIs(.hpb433("hd_betwfe"), family = "event_study")
	)

	expect_identical(sc$K, 4L)
	expect_equal(
		sc$ci$estimate,
		c(0, 1.74819372973807, 2.59385952395926, 4.91099470563771),
		tolerance = 1e-8
	)
	expect_equal(
		sc$ci$simultaneous_ci_low,
		c(0, 0.835806940880722, 1.48120298660606, 4.33044766500044),
		tolerance = 1e-8
	)
	expect_equal(
		sc$ci$simultaneous_ci_high,
		c(0, 2.66058051859542, 3.70651606131246, 5.49154174627498),
		tolerance = 1e-8
	)
	expect_equal(sc$critical_value, 2.34857048574758, tolerance = 1e-8)
	# The band really is wider than pointwise here, so the pin above is on a
	# simultaneous quantity rather than on a disguised Wald interval.
	expect_gt(sc$critical_value, sc$pointwise_critical_value)
})

# ------------------------------------------------------------------------------
# 11. (red) print() and summary() render the caveat EXACTLY ONCE on a p >= NT
#     fit -- the interactive user's channel, who never sees the warning.
#
#     Counted, not matched: `expect_match()` on captured output passes on two
#     occurrences, and double emission is the natural bug when one helper is
#     wired into two renderers and one of them renders two previews.
# ------------------------------------------------------------------------------
test_that("print/summary render the #433 caveat exactly once (betwfe)", {
	skip_on_cran()

	fit <- .hpb433("hd_betwfe")

	out <- capture.output(print(fit))
	expect_equal(sum(grepl(.HPB433_NOTICE, out, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_OTHER, out, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_FETWFE, out, fixed = TRUE)), 0L)

	out_s <- capture.output(print(summary(fit)))
	expect_equal(sum(grepl(.HPB433_NOTICE, out_s, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_OTHER, out_s, fixed = TRUE)), 1L)

	# The field the summary object carries it on. It is inside the `keep`
	# whitelist, which silently drops anything not named in it.
	s <- summary(fit)
	expect_true("highdim_band_notice" %in% names(s))
	expect_type(s$highdim_band_notice, "character")
})

# ------------------------------------------------------------------------------
# 12. (red) The same for a fetwfe fit, with the other remedy.
# ------------------------------------------------------------------------------
test_that("print/summary render the #433 caveat exactly once (fetwfe)", {
	skip_on_cran()

	fit <- .hpb433("hd_fetwfe")

	out <- capture.output(print(fit))
	expect_equal(sum(grepl(.HPB433_NOTICE, out, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_FETWFE, out, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_OTHER, out, fixed = TRUE)), 0L)

	out_s <- capture.output(print(summary(fit)))
	expect_equal(sum(grepl(.HPB433_NOTICE, out_s, fixed = TRUE)), 1L)
	expect_equal(sum(grepl(.HPB433_REMEDY_FETWFE, out_s, fixed = TRUE)), 1L)
})

# ------------------------------------------------------------------------------
# 13. (guardrail) The notice does not over-fire: a low-dimensional fit and a
#     `ci_type = "pointwise"` high-dimensional fit render nothing. The second
#     one is the interesting half -- it is `p >= N*T` and differs from the
#     firing case only in `ci_type`, so it pins that conjunct on its own.
# ------------------------------------------------------------------------------
test_that("the #433 caveat does not render on p < NT or on ci_type = 'pointwise'", {
	skip_on_cran()

	lo <- .hpb433("lo_betwfe")
	expect_equal(
		sum(grepl(.HPB433_NOTICE, capture.output(print(lo)), fixed = TRUE)),
		0L
	)
	expect_equal(
		sum(grepl(
			.HPB433_NOTICE,
			capture.output(print(summary(lo))),
			fixed = TRUE
		)),
		0L
	)
	expect_null(summary(lo)$highdim_band_notice)

	pw <- .hpb433("hd_betwfe_pointwise")
	expect_true(pw$p >= pw$N * pw$T)
	expect_identical(pw$ci_type, "pointwise")
	expect_equal(
		sum(grepl(.HPB433_NOTICE, capture.output(print(pw)), fixed = TRUE)),
		0L
	)
	expect_equal(
		sum(grepl(
			.HPB433_NOTICE,
			capture.output(print(summary(pw))),
			fixed = TRUE
		)),
		0L
	)
})

# ------------------------------------------------------------------------------
# 14. (guardrail) `isTRUE(calc_ses)` is in the gate, and this is what pins it.
#     A plain `fetwfe(gls = FALSE)` fit can be p >= N*T with
#     ci_type = "simultaneous" and `catt_df` bounds all NA -- the fit-time band
#     came back NULL and `.finalize_ci_type()` left `ci_type` alone. Without the
#     conjunct the caveat asserts a full sentence about a band that does not
#     exist, and nothing else in the suite would say so.
# ------------------------------------------------------------------------------
test_that("the #433 caveat does not render when calc_ses is FALSE", {
	skip_on_cran()

	fit <- .hpb433("gls_false_fetwfe")

	# The shape this block exists for. If any of these drift, the assertion
	# below stops meaning what it says.
	expect_true(fit$p >= fit$N * fit$T) # measured: 274 vs 192
	expect_false(fit$calc_ses)
	expect_identical(fit$ci_type, "simultaneous")
	expect_true(all(is.na(fit$catt_df$ci_low)))

	expect_equal(
		sum(grepl(.HPB433_NOTICE, capture.output(print(fit)), fixed = TRUE)),
		0L
	)
	expect_equal(
		sum(grepl(
			.HPB433_NOTICE,
			capture.output(print(summary(fit))),
			fixed = TRUE
		)),
		0L
	)
})

# ------------------------------------------------------------------------------
# A DEVIATION FROM THE PLAN, PINNED RATHER THAN LEFT ACCIDENTAL.
#
# The ExecPlan's Milestone 4 says a degenerate p >= N*T fit "renders no notice,
# for the same reason the warning skips it". The reason does not transfer, and
# the signature the same plan prescribes -- six scalars, no view of the selected
# support -- cannot express the carve-out.
#
# The warning skips that path because the path EARLY-RETURNS before the warning
# site and carries its own #304 condition. Neither half holds for the rendered
# caveat: `print()` reads a stored `catt_df`, not a control-flow position; and at
# fit time #304 is a `message()` that `.fit_band_for_family()` suppresses, so a
# user who fits and prints a degenerate high-dimensional object would otherwise
# see all-zero bounds under a `[simultaneous 95% CI]` header with no caveat at
# all. Rendering there is accurate and its remedy is the right one.
#
# So the caveat DOES render on a degenerate high-dimensional fit. That is a
# decision the maintainer may want to revisit; this block makes revisiting it
# produce a deliberate red rather than silent drift.
# ------------------------------------------------------------------------------
test_that("the #433 caveat renders on a degenerate p >= NT fit (plan deviation)", {
	skip_on_cran()

	fit <- .hpb433("degen_betwfe")
	expect_true(all(fit$catt_hats == 0))

	expect_equal(
		sum(grepl(.HPB433_NOTICE, capture.output(print(fit)), fixed = TRUE)),
		1L
	)
	expect_equal(
		sum(grepl(
			.HPB433_NOTICE,
			capture.output(print(summary(fit))),
			fixed = TRUE
		)),
		1L
	)
})
