library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Issue #470: the BEHAVIOR of the matrix-valued variance floors, which is the
# half `test-cluster_floor.R` structurally cannot see.
#
# That guardrail is lexical: it checks that every cluster-sandwich quadratic
# form sits inside a floor call. It says nothing about whether the condition
# the floor raises ever reaches a user. On the internal routes it does not, by
# default: `.assemble_joint_cov_var1()`'s conditions are signalled two frames
# below `.fit_band_for_family()`, inside a `suppressMessages()` /
# `tryCatch(error = function(e) NULL)` region, and that helper's own
# `WHAT GETS THROUGH` list already records condition suppression more than one
# frame above a floor call as out of reach. #470 creates the first LIVE
# instance of that documented blind spot.
#
# THIS FILE IS THE SOLE GUARD ON THREE THINGS, and that is why it exists:
#
#   * `for (w in pending_floor) warning(w)` is pinned by EVERY `warn = 2`
#     block here -- the direct `eventStudy()` one and both renderer ones. A
#     muffle that is never re-raised cannot error under warnings-as-errors.
#     (Measured at a91efb0: deleting the loop turns exactly blocks 1, 4b and
#     4c red, 6 failures. Not one block -- do not trim 4b or 4c as redundant
#     with 1.)
#   * `if (!is.null(fatal)) stop(fatal)` is pinned by the fit-time `fetwfe()`
#     block below. Delete it and the catastrophic tier degrades to a NULL
#     band instead of propagating. (Measured: blocks 2, 4, 4c and 6 go red.)
#   * BOTH MODEL-BASED FLOORS -- `.floor_variance_diag()` at
#     `simultaneous_cis_impl/Sigma` and `/Sigma_1` -- are pinned by block 5
#     and by nothing else. `.floor_variance_diag` is deliberately outside
#     `.clf_floor_fns`, so the guardrail cannot see those sites at all.
#     (Measured: reverting both to their pre-#470 `pmax()` form leaves
#     `test-cluster_floor.R` at 20/20 blocks green, 101 passes, and reddens
#     block 5 alone.) Two of the three floors this PR adds live or die there.
#
# Delete either re-raise and THAT TIER of the #470 diagnostic goes silent on
# every internal route -- the warning tier with the loop, the catastrophic
# tier with the `stop()` -- while `test-cluster_floor.R` stays fully green
# (measured: 20/20 blocks, 101 passes, under both mutants). Nothing but this
# file sees it. `devtools::check()` does NOT stay clean under either mutant,
# because this file is in the suite -- that is the point; the blindness is the
# guardrail's, not the gate's.
#
# HOW EVERY MESSAGE ASSERTION HERE IS WRITTEN, without exception. Both rules
# are measured, and both are named entries in the skill's fixture
# anti-patterns:
#
#   1. `fixed = TRUE` on every one of them. The pinned message template
#      contains `(`, `)`, `.` and `-`. Measured: a parenthesised literal used
#      as a REGEX matches text with the parentheses stripped out entirely
#      (TRUE where the test wants FALSE), and an unbalanced one is a hard
#      `invalid regular expression ... Missing ')'` error rather than a quiet
#      failure. `expect_error()` / `expect_warning()` pass the pattern to
#      `grepl()` under edition 2, so both behaviours reach the assertion.
#   2. Every site label is ANCHORED with its closing delimiter --
#      `site 'simultaneous_cis_impl/Sigma'`, never the bare
#      `simultaneous_cis_impl/Sigma`. `fixed = TRUE` makes a pattern literal;
#      it does not make it DISCRIMINATING. Measured:
#      `simultaneous_cis_impl/Sigma` is a strict prefix of
#      `simultaneous_cis_impl/Sigma_1`, so the unanchored form matches the
#      WRONG case's message and injection B's two cases -- which differ only
#      in which model-based floor bound -- cannot be told apart. The
#      injection-B block discharges that explicitly, each case asserting
#      `FALSE` on the other's message.
#
# AND NEVER THE BARE PHRASE `Negative cluster-sandwich quadratic form`. The
# 2026 matrix warning deliberately reuses the 2024 scalar helper's exact
# subject, so `test-cluster_floor.R`'s smoke filter widens to the new site for
# free -- but the phrase cannot tell the two apart, and the pre-existing
# scalar sites sit on the SAME public routes, one frame earlier. Key on the
# anchored site label.
#
# WHY THE MOCK IS CALLER-GATED, and why that gate is not decoration.
# `.recompute_gram_and_sandwich()` is called twice per public call: once for
# `.event_study_fetwfe()` / `getCohortATTsFinal()`, which already carry the
# #139 SCALAR diagnostic, and once from inside `.simultaneous_cis_impl()`.
# An UNGATED mock fires the 2024 diagnostic first, so `expect_error()` passes
# on the unmodified tree for reasons that have nothing to do with #470 --
# measured, and it invalidated an earlier draft of two of these blocks. If a
# block here ever passes on a base tree, suspect the gate before the
# assertion, and read the message: a base-tree pass names
# `getCohortATTsFinal/cohort_te_se` or `event_study_fetwfe/var_1_e`, never
# `assemble_joint_cov_var1/Sigma_1`.
#
# `cohortTimeATTs()` and `cohortStudy()` also call
# `.recompute_gram_and_sandwich()`, but not via `.simultaneous_cis_impl()`, so
# the caller-aware gate leaves them untouched by construction. No block is
# needed for them; noted so the next reader does not re-derive it.
# ------------------------------------------------------------------------------

# The anchored site labels, written once. Every assertion below uses these.
.MFC470_SANDWICH <- "site 'assemble_joint_cov_var1/Sigma_1'"
.MFC470_SIGMA <- "site 'simultaneous_cis_impl/Sigma'"
.MFC470_SIGMA1 <- "site 'simultaneous_cis_impl/Sigma_1'"

.mfc470_orig <- fetwfe:::.recompute_gram_and_sandwich

# Negate `sandwich_full` ONLY when `gate_fn` is on the call stack. See the
# header: the gate is what keeps these blocks from passing on the base tree by
# firing the pre-existing scalar diagnostic. A call-index gate (`i == 2L`)
# also works but depends on call ordering staying what it is.
.mfc470_neg_sandwich <- function(gate_fn, scale) {
	function(...) {
		g <- .mfc470_orig(...)
		on_gate <- any(vapply(
			sys.calls(),
			function(cl) identical(deparse(cl[[1]])[1], gate_fn),
			logical(1)
		))
		if (on_gate && !is.null(g$sandwich_full)) {
			g$sandwich_full <- -scale * g$sandwich_full
		}
		g
	}
}

# Injection B: negate `gram_inv`, which is what the MODEL-BASED branch of
# `.assemble_joint_cov_var1()` reads. Injection A never reaches the two
# model-based floors, because that branch does not touch `sandwich_full` at
# all. A direct `simultaneousCIs()` call reaches the helper exactly once, so
# this injection needs no call gate.
.mfc470_neg_gram <- function(...) {
	g <- .mfc470_orig(...)
	g$gram_inv <- -g$gram_inv
	g
}

.mfc470_cache <- new.env(parent = emptyenv())

.mfc470_build <- function(name) {
	sim <- simulateData(
		genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 7),
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 7
	)
	switch(
		name,
		sim = sim,
		cluster = fetwfeWithSimulatedData(sim, se_type = "cluster"),
		default = fetwfeWithSimulatedData(sim, se_type = "default"),
		# `fetwfeWithSimulatedData()` can NEVER reach the `v1` site: it
		# forwards `sim$indep_counts`, so `indep_counts_used` is TRUE and
		# `is_indep` short-circuits `if (is_indep || !identical(se_type,
		# "conservative"))` into the tight branch. The conservative branch
		# needs `indep_counts = NA` and a hand-built call.
		conservative = fetwfe(
			pdata = sim$pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			response = sim$response,
			covs = sim$covs,
			indep_counts = NA,
			se_type = "conservative",
			sig_eps_sq = sim$sig_eps_sq,
			sig_eps_c_sq = sim$sig_eps_c_sq
		),
		stop("unknown fixture: ", name)
	)
}

.mfc470 <- function(name) {
	if (!exists(name, envir = .mfc470_cache, inherits = FALSE)) {
		assign(name, .mfc470_build(name), envir = .mfc470_cache)
	}
	get(name, envir = .mfc470_cache, inherits = FALSE)
}

.mfc470_skip <- function() {
	skip_if(
		packageVersion("testthat") < "3.2.0",
		"with_mocked_bindings() requires testthat >= 3.2.0"
	)
}

# ------------------------------------------------------------------------------
# 1. (red) THE LOAD-BEARING ONE. Under `options(warn = 2)`, an `eventStudy()`
#    call on a fit whose joint-covariance sandwich has been negated must NOT
#    return the pointwise band under a `[simultaneous 95% CI]` header. It must
#    error.
#
#    This block pins `for (w in pending_floor) warning(w)` in
#    `.fit_band_for_family()`. The muffle inside the `tryCatch()` beats the
#    `warn = 2` conversion, so without the re-raise below it nothing errors and
#    nothing warns.
#
#    DO NOT WRAP THE warn = 2 CALL IN suppressWarnings() -- nor in a muffling
#    withCallingHandlers(warning = ), which is the form a mocking test reaches
#    for and which defeats the conversion the same way. Measured during the
#    plan review: it reports "no error" where the bare call errors.
# ------------------------------------------------------------------------------
test_that("warn = 2 eventStudy errors rather than downgrading the band (#470)", {
	.mfc470_skip()
	fit <- .mfc470("cluster")

	# ANTI-VACUITY. If the simultaneous and pointwise bands coincided on this
	# fixture, an equality check could not tell a downgrade from a no-op.
	ref <- eventStudy(fit)
	pw <- eventStudy(fit, ci_type = "pointwise")
	expect_false(isTRUE(all.equal(ref$ci_low, pw$ci_low)))

	old <- options(warn = 2)
	on.exit(options(old), add = TRUE)

	err <- testthat::with_mocked_bindings(
		tryCatch(eventStudy(fit), error = function(e) e),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".simultaneous_cis_impl",
			1
		),
		.package = "fetwfe"
	)
	expect_s3_class(err, "error")
	expect_true(grepl(.MFC470_SANDWICH, conditionMessage(err), fixed = TRUE))
	# It is the #470 warning tier that was converted, not some other condition.
	expect_true(grepl(
		"clipped to 0 (indices",
		conditionMessage(err),
		fixed = TRUE
	))
})

# ------------------------------------------------------------------------------
# 2. (red) The catastrophic tier PROPAGATES out of a fit-time `fetwfe()` call
#    rather than degrading to a NULL band. This block pins
#    `if (!is.null(fatal)) stop(fatal)`.
#
#    `scale = 100` puts the Sigma_1 diagonal past -1 on this fixture.
#    Measured on this tree: most negative -1.68, three of three entries.
# ------------------------------------------------------------------------------
test_that("the catastrophic tier propagates out of a fetwfe() fit (#470)", {
	.mfc470_skip()
	sim <- .mfc470("sim")

	err <- testthat::with_mocked_bindings(
		tryCatch(
			fetwfeWithSimulatedData(sim, se_type = "cluster"),
			error = function(e) e
		),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".simultaneous_cis_impl",
			100
		),
		.package = "fetwfe"
	)
	expect_s3_class(err, "fetwfe_negative_variance_catastrophic")
	msg <- conditionMessage(err)
	expect_true(grepl(.MFC470_SANDWICH, msg, fixed = TRUE))
	expect_true(grepl("catastrophically negative", msg, fixed = TRUE))
	expect_true(grepl("entries below -1 (indices", msg, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 3. A direct `simultaneousCIs(family = "cohort")` call at `warn = 0` emits
#    EXACTLY ONE warning naming the site and its offending indices, and still
#    returns a band whose simultaneous interval has collapsed to a point at the
#    estimate.
#
#    NOT `family = "event_study"`: negating the sandwich breaks the
#    off-diagonals too, so the pre-existing invalid-correlation-matrix guard
#    fires there first on BOTH trees. A sweep of milder `sandwich - cI`
#    perturbations found no window in which the diagonal goes negative and the
#    correlation matrix stays usable, so there is no version of this block that
#    asserts "the band is still the simultaneous one" under injection A. That
#    assertion lives under injection B instead.
#
#    The muffling handler here is safe precisely because this block does NOT
#    run under `warn = 2` -- it is counting conditions, not testing the
#    conversion. Block 1 is the one where a muffle would hide the failure.
# ------------------------------------------------------------------------------
test_that("simultaneousCIs warns once and returns a collapsed band (#470)", {
	.mfc470_skip()
	fit <- .mfc470("cluster")

	# ANTI-VACUITY. On the clean tree the band is NOT collapsed, so the
	# collapse asserted below is caused by the injection.
	clean <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	expect_true(any(
		clean$ci$simultaneous_ci_high > clean$ci$simultaneous_ci_low
	))

	ws <- list()
	out <- testthat::with_mocked_bindings(
		withCallingHandlers(
			suppressMessages(simultaneousCIs(fit, family = "cohort")),
			warning = function(w) {
				ws[[length(ws) + 1L]] <<- w
				invokeRestart("muffleWarning")
			}
		),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".simultaneous_cis_impl",
			1
		),
		.package = "fetwfe"
	)

	# ONE aggregated condition, not one per offending diagonal entry.
	expect_length(ws, 1L)
	expect_s3_class(ws[[1]], "fetwfe_negative_variance_floored")
	msg <- conditionMessage(ws[[1]])
	expect_true(grepl(.MFC470_SANDWICH, msg, fixed = TRUE))
	expect_true(grepl("(indices 1, 2, 3;", msg, fixed = TRUE))

	# The band still comes back, with every interval collapsed to a point at
	# the estimate. (The returned frame has no `se` column, which is why the
	# criterion is phrased as interval collapse.)
	expect_true(is.data.frame(out$ci))
	expect_equal(out$ci$simultaneous_ci_low, out$ci$estimate)
	expect_equal(out$ci$simultaneous_ci_high, out$ci$estimate)
})

# ------------------------------------------------------------------------------
# 4. (red) THE RENDERER ROUTE. Same machinery, one string changed: the gate is
#    `.event_study_simultaneous_bounds` instead of `.simultaneous_cis_impl`, so
#    only the EVENT-STUDY family is perturbed.
#
#    This is the block behind the plan's `Decision Log` entry on staying loud,
#    and without it that entry is an assertion rather than a tested behavior.
#    The fit-time precompute runs `.fit_band_for_family(x, "cohort", ...)` and
#    the renderers run `.fit_band_for_family(x, "event_study", ...)`; the two
#    build DIFFERENT `Psi` matrices over the same `sandwich_full`, so a broken
#    PSD invariant can bind on one and not the other. Measured here: the fit
#    returns with ZERO warnings and `print()` / `summary()` then `stop()`.
#
#    That is a user-visible behavior change: a successfully-fitted object can
#    fail to print. It is deliberate -- `.fit_band_for_family()` returning NULL
#    falls through to the POINTWISE bounds under a `[simultaneous 95% CI]`
#    header, which is narrower and over-rejects, i.e. the #433 wrong-answer
#    shape in the function #433 rewrote to close it. The object stays
#    inspectable while a user diagnoses: `fit$catt_df`, `fit$att_hat` and
#    `eventStudy(fit, ci_type = "pointwise")` all still work.
#
#    `.event_study_quiet()`'s muffle is keyed to the
#    `fetwfe_highdim_postselection_band` class ONLY, which is what lets both
#    new conditions through to the renderers. DO NOT widen it to a blanket
#    muffle to make them quiet again: a class-blind muffle there re-opens
#    #470's own defect on three more doors.
# ------------------------------------------------------------------------------
test_that("print/summary stay loud on an event-study-only breakdown (#470)", {
	.mfc470_skip()
	sim <- .mfc470("sim")
	mock <- .mfc470_neg_sandwich(".event_study_simultaneous_bounds", 100)

	fit <- testthat::with_mocked_bindings(
		{
			ws <- list()
			f <- withCallingHandlers(
				fetwfeWithSimulatedData(sim, se_type = "cluster"),
				warning = function(w) {
					ws[[length(ws) + 1L]] <<- w
					invokeRestart("muffleWarning")
				}
			)
			# The FIT ITSELF succeeds, silently: the cohort family's contrasts
			# keep their diagonal non-negative under this perturbation.
			expect_length(ws, 0L)
			f
		},
		.recompute_gram_and_sandwich = mock,
		.package = "fetwfe"
	)

	testthat::with_mocked_bindings(
		{
			for (render in list(
				function() capture.output(print(fit)),
				function() capture.output(print(summary(fit)))
			)) {
				err <- tryCatch(render(), error = function(e) e)
				expect_s3_class(err, "fetwfe_negative_variance_catastrophic")
				expect_true(grepl(
					.MFC470_SANDWICH,
					conditionMessage(err),
					fixed = TRUE
				))
			}
			# ... while the object stays inspectable.
			expect_true(is.data.frame(fit$catt_df))
			expect_true(is.data.frame(eventStudy(fit, ci_type = "pointwise")))
		},
		.recompute_gram_and_sandwich = mock,
		.package = "fetwfe"
	)
})

# ------------------------------------------------------------------------------
# 4b. (red) The same three renderers error at the WARNING tier under
#     `options(warn = 2)` (`scale = 1`). Bare, per the no-`suppressWarnings()`
#     rule in block 1's comment.
# ------------------------------------------------------------------------------
test_that("warn = 2 makes the warning tier reach the renderers too (#470)", {
	.mfc470_skip()
	fit <- .mfc470("cluster")
	mock <- .mfc470_neg_sandwich(".event_study_simultaneous_bounds", 1)

	old <- options(warn = 2)
	on.exit(options(old), add = TRUE)

	testthat::with_mocked_bindings(
		{
			for (render in list(
				function() capture.output(print(fit)),
				function() capture.output(print(summary(fit)))
			)) {
				err <- tryCatch(render(), error = function(e) e)
				expect_s3_class(err, "error")
				expect_true(grepl(
					.MFC470_SANDWICH,
					conditionMessage(err),
					fixed = TRUE
				))
				expect_true(grepl(
					"clipped to 0 (indices",
					conditionMessage(err),
					fixed = TRUE
				))
			}
		},
		.recompute_gram_and_sandwich = mock,
		.package = "fetwfe"
	)
})

# ------------------------------------------------------------------------------
# 4c. `plot()` is the third renderer. Split into its own block, rather than
#     folded into 4 and 4b, so that a machine without `ggplot2` -- a Suggest --
#     skips only this and never one of the load-bearing blocks. `plot.fetwfe()`
#     checks `requireNamespace("ggplot2")` first and stops with its own message
#     if it is absent, so the assertion would otherwise be meaningless there.
# ------------------------------------------------------------------------------
test_that("plot() stays loud on an event-study-only breakdown (#470)", {
	.mfc470_skip()
	skip_if_not_installed("ggplot2")
	fit <- .mfc470("cluster")

	err <- testthat::with_mocked_bindings(
		tryCatch(print(plot(fit, type = "event")), error = function(e) e),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".event_study_simultaneous_bounds",
			100
		),
		.package = "fetwfe"
	)
	expect_s3_class(err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.MFC470_SANDWICH, conditionMessage(err), fixed = TRUE))

	old <- options(warn = 2)
	on.exit(options(old), add = TRUE)
	err2 <- testthat::with_mocked_bindings(
		tryCatch(print(plot(fit, type = "event")), error = function(e) e),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".event_study_simultaneous_bounds",
			1
		),
		.package = "fetwfe"
	)
	expect_s3_class(err2, "error")
	expect_true(grepl(.MFC470_SANDWICH, conditionMessage(err2), fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 5. (red) INJECTION B -- the two MODEL-BASED floors. Negate `gram_inv` rather
#    than `sandwich_full`; injection A never reaches these, because the
#    model-based branch of `.assemble_joint_cov_var1()` does not read the
#    sandwich at all. This is also where "a warning fires and a band is still
#    returned" lives, since injection A cannot support it (see block 3).
#
#    The two cases differ ONLY in which model-based floor bound, and their site
#    labels are in a strict prefix relation, so this block is the reason rule 2
#    in the header exists. Each case asserts the other's anchored literal is
#    FALSE on its own message -- the catalogue's own discharge for a
#    non-discriminating assertion.
# ------------------------------------------------------------------------------
test_that("the model-based floors warn and still return a band (#470)", {
	.mfc470_skip()

	grab <- function(fit) {
		# FORCE the fixture before the mock is installed. `grab(.mfc470("x"))`
		# passes a promise, and R would otherwise evaluate it at its first use
		# -- which is INSIDE `with_mocked_bindings()`, so the fit itself would
		# be built on a negated Gram inverse. Measured: eight warnings from a
		# garbage fit instead of the one this block is about.
		force(fit)
		ws <- list()
		out <- testthat::with_mocked_bindings(
			withCallingHandlers(
				suppressMessages(simultaneousCIs(fit, family = "cohort")),
				warning = function(w) {
					ws[[length(ws) + 1L]] <<- w
					invokeRestart("muffleWarning")
				}
			),
			.recompute_gram_and_sandwich = .mfc470_neg_gram,
			.package = "fetwfe"
		)
		expect_length(ws, 1L)
		expect_s3_class(ws[[1]], "fetwfe_negative_variance_floored")
		expect_true(is.data.frame(out$ci))
		conditionMessage(ws[[1]])
	}

	# se_type = "default" -> the combined `Sigma` diagonal.
	msg_d <- grab(.mfc470("default"))
	expect_true(grepl(.MFC470_SIGMA, msg_d, fixed = TRUE))
	# The subject is "variance", not the sandwich phrase: these diagonals are
	# model-based, so the #139 wording would be false about them.
	expect_true(grepl("Negative variance", msg_d, fixed = TRUE))

	# se_type = "conservative" with indep_counts = NA -> the `Sigma_1`
	# diagonal, in the Cauchy-Schwarz branch.
	fit_c <- .mfc470("conservative")
	expect_false(isTRUE(fit_c$indep_counts_used))
	msg_c <- grab(fit_c)
	expect_true(grepl(.MFC470_SIGMA1, msg_c, fixed = TRUE))

	# THE DISCRIMINATING DISCHARGE. `simultaneous_cis_impl/Sigma` is a strict
	# prefix of `simultaneous_cis_impl/Sigma_1`; anchored with the closing
	# quote, each literal is FALSE on the other case's message. Unanchored,
	# the first would match both and this block would be unable to fail.
	expect_false(grepl(.MFC470_SIGMA1, msg_d, fixed = TRUE))
	expect_false(grepl(.MFC470_SIGMA, msg_c, fixed = TRUE))
	# ... and the unanchored form really is the trap, not a hypothetical one.
	expect_true(grepl("simultaneous_cis_impl/Sigma", msg_c, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 6. (red) ONE CROSS-ESTIMATOR BLOCK. Everything above is `fetwfe`. `etwfe()`,
#    `betwfe()` and `twfeCovs()` reach the same `.finalize_ci_type()` ->
#    `.fit_band_for_family()` -> `.assemble_joint_cov_var1()` code, which is
#    estimator-agnostic, but `betwfe` / `twfeCovs` route through
#    `.event_study_etwfe_betwfe()` rather than `.event_study_fetwfe()`. This
#    checks the ROUTE, not the arithmetic.
# ------------------------------------------------------------------------------
test_that("the catastrophic tier propagates on a non-fetwfe fit too (#470)", {
	.mfc470_skip()
	sim <- .mfc470("sim")

	err <- testthat::with_mocked_bindings(
		tryCatch(
			betwfeWithSimulatedData(sim, se_type = "cluster"),
			error = function(e) e
		),
		.recompute_gram_and_sandwich = .mfc470_neg_sandwich(
			".simultaneous_cis_impl",
			100
		),
		.package = "fetwfe"
	)
	expect_s3_class(err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.MFC470_SANDWICH, conditionMessage(err), fixed = TRUE))
})
