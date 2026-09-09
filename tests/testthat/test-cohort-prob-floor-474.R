library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Issue #474: the BEHAVIOR of the COHORT-PROBABILITY variance floors.
#
# `att_var_2` / `Sigma_2` is the extra variance inherited from ESTIMATING,
# rather than knowing, what fraction of units belongs to each treated cohort.
# It is non-negative in exact arithmetic (`theta' J' Sigma_pi_hat J theta` with
# `Sigma_pi_hat` a multinomial covariance) and can come out negative in
# floating point. Before #474 three sites floored it with a BARE `max(., 0)` /
# `pmax(., 0)` and said nothing, and a fourth had no floor at all.
#
# THE FOUR SITES, and their labels, written once as constants below:
#
#   1. `.assemble_joint_cov_var2()`  (`R/variance_machinery.R`) -- the K x K
#      `Sigma_2` diagonal. THE ONLY ONE INSIDE `.fit_band_for_family()`'s
#      protected region, and the only one on the simultaneous-band path, so it
#      is the only one whose catastrophic message carries the band remedy.
#   2. `getSecondVarTermOLS()`       (`R/variance_machinery.R`) -- scalar,
#      OLS family (`etwfe` / `betwfe` / `twfeCovs`).
#   3. `getSecondVarTermDataApp()`   (`R/variance_machinery.R`) -- scalar,
#      FETWFE family.
#   4. `.event_study_var2_fetwfe()`  (`R/event_study.R`) -- scalar, per event
#      time. THE ONE THAT HAD NO FLOOR: a negative reached `sqrt()` and gave
#      `NaN`, or -- when `var_1(e)` kept the sum positive -- a silently
#      UNDERSTATED event-study standard error. Blocks 6 and 7 below are that
#      defect's regression tests, and they pin the VALUE, not just its
#      finiteness. Measured on the pre-#474 tree, one fixture shows both
#      shapes at once: `0.0501, 0.0501, NaN, 0.0804, 0.2128` where the
#      floored answer is `0.0954, 0.0954, 0.1088, 0.1238, 0.2128`.
#
# THIS FILE IS THE ENTIRE DETECTOR FOR THE CHANGE. Measured in the
# pre-implementation pass: with `.floor_cohort_prob_var()`'s body replaced by
# `pmax(v, 0)` -- all four sites silently reverting to their pre-#474 bare
# floor -- every other test file in the suite stays byte-identically green,
# `test-cluster_floor.R`'s A5c pin included, because A5c pins FORMALS and a
# body change does not move them. That is why the route coverage below is
# load-bearing rather than belt-and-braces: the family is deliberately outside
# `test-cluster_floor.R`'s recognized floor-function set `.clf_floor_fns`
# (adding it turns A9 red on four functions and A7b red on
# `.fit_band_for_family`; tracked as issue #478), so the lexical guardrail
# cannot see these sites at all.
#
# WHY EVERY TIER IS REACHED BY INJECTION. No naturally-occurring fixture
# produces a negative at any of these sites -- measured over ten fits,
# including the paper's `divorce` application, where this variance piece is
# exactly zero because every cohort fuses. So the battery injects, exactly as
# `test-matrix-floor-conditions-470.R` does.
#
# WHY THE INJECTION IS CALLER-GATED, AND WHY THE GATE DIFFERS PER SITE.
# Profile gotcha 12.8: mocking a helper that feeds two routes perturbs an
# earlier site one frame up and makes assertions pass on unmodified code.
# Sites 2, 3 and 4 build `Sigma_pi_hat` themselves via `.multinomial_cov()`;
# site 1 RECEIVES it as an argument from `.simultaneous_cis_impl()`. So one
# `.multinomial_cov()` mock serves all four, but the gate function is
# per-site: the site function itself for 2/3/4, and `.simultaneous_cis_impl` /
# `.event_study_simultaneous_bounds` for site 1's two live routes. Negating
# `Sigma_pi_hat` makes the quadratic form exactly `-scale x` its true value,
# so the tier reached is linear in `scale` -- which is why the two scales
# below are a measurement rather than a guess.
#
# HOW EVERY MESSAGE ASSERTION IS WRITTEN, inheriting
# `test-matrix-floor-conditions-470.R`'s two rules unchanged:
#
#   1. `fixed = TRUE` on every one. The message template contains `(`, `)`,
#      `.` and `-`; used as a regex a parenthesised literal matches text with
#      the parentheses stripped.
#   2. Every site label is ANCHORED with its closing quote --
#      `site 'getSecondVarTermOLS/att_var_2'`, never the bare label. There is
#      no prefix collision among the four new labels or against the six
#      existing ones, but `getSecondVarTermOLS/att_var_2` and
#      `getSecondVarTermDataApp/att_var_2` share the `att_var_2` SUFFIX, so an
#      assertion keyed on that alone cannot tell the two scalar sites apart --
#      and the coverage predicate here is per-site. The two blocks that cover
#      them assert the other's label is FALSE on their own message.
#
# THE REMEDY-CLAUSE CONTRAST is the assertion that keeps the false-remedy
# regression from returning silently. `.floor_psd_diag_core()`'s catastrophic
# tier hardcoded "refit with ci_type = \"pointwise\"" until #474. Measured:
# `ci_type = "pointwise"` does not avoid sites 2, 3 or 4 and does not even
# reduce their call count, so a user already on `pointwise` would be told to
# refit with `pointwise` and get the identical error. Every catastrophic-tier
# block below therefore asserts the clause PRESENT at site 1 and ABSENT at
# sites 2, 3 and 4.
#
# FOUR AXIS CELLS ARE OUTSIDE THIS BATTERY ON PURPOSE, each stated so the
# omission does not read as an oversight:
#
#   * `q >= 1` reaches NONE of the four sites -- `att_se` is `NA` by design
#     and the variance machinery does not run -- even though profile gotcha
#     12.6 names `q` as an interacting option.
#   * `simultaneousCIs(method = "bootstrap")` reaches none of the family
#     either. "The analytic path" is a scope boundary here, not an emphasis.
#   * `twfeCovs` has no event-study route at all (`print()` / `summary()` on a
#     `twfeCovs` fit reach none of the four sites), so a battery looping
#     uniformly over four estimators would write dead cells. It is covered at
#     fit time only.
#   * `se_type` across all three values, `gls`, `add_ridge` and `indep_counts`
#     leave the reached-site set IDENTICAL, so gotcha 12.6's interaction
#     matrix collapses to this line rather than needing a table. (`se_type`
#     does change the arithmetic DOWNSTREAM of site 4 -- the conservative
#     branch takes `2 * sqrt(var_1 * var_2)` -- which is what block 7 is
#     about.)
# ------------------------------------------------------------------------------

# The anchored site labels, written once. Every assertion below uses these.
.CPF474_SIGMA2 <- "site 'assemble_joint_cov_var2/Sigma_2'"
.CPF474_OLS <- "site 'getSecondVarTermOLS/att_var_2'"
.CPF474_DATAAPP <- "site 'getSecondVarTermDataApp/att_var_2'"
.CPF474_ES <- "site 'event_study_var2_fetwfe/var_2_e'"
# The remedy clause, and the family's own subject noun.
.CPF474_REMEDY <- "refit with ci_type = \"pointwise\""
.CPF474_SUBJECT <- "cohort-probability variance"
# The core sentence-cases `subject` in the CATASTROPHIC tier only (it opens
# that message), so the two tiers render the noun differently and an assertion
# keyed on the lowercase form silently fails at the error tier.
#
# The casing behavior ITSELF is already pinned, by two assertions in
# `test-cluster_floor.R` -- removing the core's `toupper()` step reddens those
# as well as the catastrophic-tier check here. What is new and unique to this
# file is the catastrophic-tier SUBJECT check on the RENDERER route, which is
# what a sibling-wrapper swap at site 1 drops. Measured both directions: the
# lowercase and capitalized assertions each catch a mutation the other misses.
.CPF474_SUBJECT_CAP <- "Cohort-probability variance"

# Measured on this fixture, not guessed. The quadratic form is exactly
# `-scale x` its true value, and the true values run 3.7e-3 .. 4.2e-2 across
# the four sites, so `scale = 1` lands every site in the warning band
# (-1, -1e-10) and `scale = 1e4` puts every site past -1 (measured most
# negative: -145 .. -423).
.CPF474_WARN_SCALE <- 1
.CPF474_ERR_SCALE <- 1e4
# Site 1's DIRECT `simultaneousCIs()` route needs a gentler warning-tier
# perturbation: `.assemble_joint_cov_var2()` floors only the DIAGONAL, so at
# `scale = 1` the untouched off-diagonals push `Sigma` past the pre-existing
# "Covariance matrix not positive semidefinite" guard and that error arrives
# first. At 1e-3 the diagonal is still -4.2e-05 -- comfortably inside the
# warning band -- and the off-diagonals stay usable. The catastrophic tier
# needs no such care: it `stop()`s inside `.assemble_joint_cov_var2()`, before
# `Sigma` is formed.
.CPF474_WARN_SCALE_DIRECT <- 1e-3

.cpf474_orig <- fetwfe:::.multinomial_cov

# Replace `Sigma_pi_hat` with `-scale * Sigma_pi_hat` ONLY when `gate_fn` is on
# the call stack. See the header for why the gate is not decoration.
.cpf474_neg <- function(gate_fn, scale) {
	force(gate_fn)
	force(scale)
	function(probs) {
		out <- .cpf474_orig(probs)
		if (.cpf474_on_gate(gate_fn)) {
			out <- -scale * out
		}
		out
	}
}

# The REFERENCE mock: `Sigma_pi_hat` is zeroed at `gate_fn`, so the quadratic
# form there is exactly `0` through the real code path. This is how blocks 7
# and 8 obtain `sqrt(var_1(e))` without recomputing it -- a separate mechanism
# from the value under test (a negative clipped to zero by the floor), rather
# than a round-trip through the same one.
.cpf474_zero <- function(gate_fn) {
	force(gate_fn)
	function(probs) {
		out <- .cpf474_orig(probs)
		if (.cpf474_on_gate(gate_fn)) {
			out <- 0 * out
		}
		out
	}
}

.cpf474_on_gate <- function(gate_fn) {
	any(vapply(
		sys.calls(),
		function(cl) identical(deparse(cl[[1]])[1], gate_fn),
		logical(1)
	))
}

.cpf474_cache <- new.env(parent = emptyenv())

.cpf474_build <- function(name) {
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
		fetwfe = fetwfeWithSimulatedData(sim, se_type = "cluster"),
		etwfe = etwfeWithSimulatedData(sim, se_type = "cluster"),
		betwfe = betwfeWithSimulatedData(sim, se_type = "cluster"),
		# `fetwfeWithSimulatedData()` forwards `sim$indep_counts`, so
		# `indep_counts_used` is TRUE and `.combine_att_variance()` takes its
		# TIGHT branch whatever `se_type` says. The Cauchy-Schwarz branch --
		# the one where a negative `att_var_2` reaches
		# `2 * sqrt(att_var_1 * att_var_2)` -- needs `indep_counts = NA` and a
		# hand-built call. (Same reasoning, same shape, as
		# `test-matrix-floor-conditions-470.R`'s `conservative` fixture.)
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

.cpf474 <- function(name) {
	if (!exists(name, envir = .cpf474_cache, inherits = FALSE)) {
		assign(name, .cpf474_build(name), envir = .cpf474_cache)
	}
	get(name, envir = .cpf474_cache, inherits = FALSE)
}

.cpf474_skip <- function() {
	skip_if(
		packageVersion("testthat") < "3.2.0",
		"with_mocked_bindings() requires testthat >= 3.2.0"
	)
}

# Run `expr` under a gated injection, collecting every warning and any error.
# Warnings are recorded by an OBSERVE-AND-MUFFLE handler, which is safe here
# because none of these blocks runs under `options(warn = 2)`; the `warn = 2`
# blocks call bare, per `test-matrix-floor-conditions-470.R` block 1's rule.
.cpf474_run <- function(mock, expr_fn) {
	ws <- list()
	err <- NULL
	testthat::with_mocked_bindings(
		withCallingHandlers(
			tryCatch(expr_fn(), error = function(e) err <<- e),
			warning = function(w) {
				ws[[length(ws) + 1L]] <<- w
				invokeRestart("muffleWarning")
			}
		),
		.multinomial_cov = mock,
		.package = "fetwfe"
	)
	list(ws = ws, err = err)
}

# The floor-family warnings only, so an unrelated warning cannot satisfy a
# count. Keyed on the CLASS, which is what the classing exists for.
.cpf474_floor_ws <- function(r) {
	Filter(function(w) inherits(w, "fetwfe_negative_variance_floored"), r$ws)
}

# ------------------------------------------------------------------------------
# 1. (red) SITE 2 -- `getSecondVarTermOLS()` at FIT TIME, both tiers, on all
#    three OLS-family estimators. This is the fit-time standard-error route,
#    and it is the cell that gains its FIRST fit-blocking variance diagnostic
#    here: `.compute_att_var1()`'s #139 `stop()` tier is inside
#    `if (identical(se_type, "cluster"))`, and #470's is on the simultaneous
#    band, so `se_type = "default"` with `ci_type = "pointwise"` had none.
# ------------------------------------------------------------------------------
test_that("site 2 fires both tiers at fit time on the OLS family (#474)", {
	.cpf474_skip()
	sim <- .cpf474("sim")

	fits <- list(
		etwfe = function() etwfeWithSimulatedData(sim, se_type = "cluster"),
		betwfe = function() betwfeWithSimulatedData(sim, se_type = "cluster"),
		twfeCovs = function() {
			twfeCovsWithSimulatedData(sim, se_type = "cluster")
		}
	)

	for (nm in names(fits)) {
		# --- warning tier ---
		r <- .cpf474_run(
			.cpf474_neg("getSecondVarTermOLS", .CPF474_WARN_SCALE),
			fits[[nm]]
		)
		fws <- .cpf474_floor_ws(r)
		expect_gt(length(fws), 0L)
		msg <- conditionMessage(fws[[1]])
		expect_s3_class(fws[[1]], "fetwfe_negative_variance_floored")
		expect_true(grepl(.CPF474_OLS, msg, fixed = TRUE))
		expect_true(grepl(.CPF474_SUBJECT, msg, fixed = TRUE))
		# The two scalar labels share the `att_var_2` suffix; anchored, each is
		# FALSE on the other's message.
		expect_false(grepl(.CPF474_DATAAPP, msg, fixed = TRUE))

		# --- catastrophic tier ---
		r <- .cpf474_run(
			.cpf474_neg("getSecondVarTermOLS", .CPF474_ERR_SCALE),
			fits[[nm]]
		)
		expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
		emsg <- conditionMessage(r$err)
		expect_true(grepl(.CPF474_OLS, emsg, fixed = TRUE))
		expect_true(grepl("catastrophically negative", emsg, fixed = TRUE))
		expect_true(grepl("entries below -1 (indices", emsg, fixed = TRUE))
		# NO REMEDY CLAUSE. This site is not on the band path, and
		# `ci_type = "pointwise"` neither avoids it nor reduces its call
		# count, so the pre-#474 hardcoded sentence would be false here.
		expect_false(grepl(.CPF474_REMEDY, emsg, fixed = TRUE))
	}
})

# ------------------------------------------------------------------------------
# 2. (red) SITE 3 -- `getSecondVarTermDataApp()` at FIT TIME on `fetwfe()`,
#    both tiers. The FETWFE-family counterpart of block 1; measured, a
#    `fetwfe()` fit reaches THIS site and never `getSecondVarTermOLS()`, and
#    the OLS-family fits reach the other and never this one.
# ------------------------------------------------------------------------------
test_that("site 3 fires both tiers at fit time on fetwfe (#474)", {
	.cpf474_skip()
	sim <- .cpf474("sim")
	fit_it <- function() fetwfeWithSimulatedData(sim, se_type = "cluster")

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermDataApp", .CPF474_WARN_SCALE),
		fit_it
	)
	fws <- .cpf474_floor_ws(r)
	expect_gt(length(fws), 0L)
	msg <- conditionMessage(fws[[1]])
	expect_s3_class(fws[[1]], "fetwfe_negative_variance_floored")
	expect_true(grepl(.CPF474_DATAAPP, msg, fixed = TRUE))
	expect_true(grepl(.CPF474_SUBJECT, msg, fixed = TRUE))
	expect_false(grepl(.CPF474_OLS, msg, fixed = TRUE))

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermDataApp", .CPF474_ERR_SCALE),
		fit_it
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	emsg <- conditionMessage(r$err)
	expect_true(grepl(.CPF474_DATAAPP, emsg, fixed = TRUE))
	expect_true(grepl("catastrophically negative", emsg, fixed = TRUE))
	expect_false(grepl(.CPF474_OLS, emsg, fixed = TRUE))
	expect_false(grepl(.CPF474_REMEDY, emsg, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 2b. (red) SITE 3 at fit time again, this time at `se_type = "default"` with
#     `ci_type = "pointwise"` -- the cell the `NEWS.md` breaking-change bullet
#     singles out as the one that had NO fit-blocking variance diagnostic
#     before this change. `.compute_att_var1()`'s #139 `stop()` tier sits
#     inside `if (identical(se_type, "cluster"))` and #470's sits on the
#     simultaneous band, so neither reaches here.
#
#     It exists because that distinction is a user-facing claim in
#     `### Breaking changes` and every other fit-time block in this file fits
#     at `se_type = "cluster"`, which left the claim resting on prose. An
#     assertion is cheap; a NEWS bullet nobody can falsify is not.
# ------------------------------------------------------------------------------

test_that("site 3 fires at fit time on se_type = default, ci_type = pointwise (#474)", {
	.cpf474_skip()
	sim <- .cpf474("sim")
	fit_it <- function() {
		fetwfeWithSimulatedData(
			sim,
			se_type = "default",
			ci_type = "pointwise"
		)
	}

	# Control: this cell fits clean on unperturbed data, so the assertions
	# below are about the injection rather than about the cell being broken.
	expect_no_error(fit_it())

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermDataApp", .CPF474_WARN_SCALE),
		fit_it
	)
	fws <- .cpf474_floor_ws(r)
	expect_gt(length(fws), 0L)
	expect_s3_class(fws[[1]], "fetwfe_negative_variance_floored")
	expect_true(grepl(
		.CPF474_DATAAPP,
		conditionMessage(fws[[1]]),
		fixed = TRUE
	))

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermDataApp", .CPF474_ERR_SCALE),
		fit_it
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	emsg <- conditionMessage(r$err)
	expect_true(grepl(.CPF474_DATAAPP, emsg, fixed = TRUE))
	# `ci_type = "pointwise"` is already set, so a remedy clause telling the
	# user to set it would be the false-remedy defect this PR exists to avoid.
	expect_false(grepl(.CPF474_REMEDY, emsg, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 3. (red) SITE 2 via `eventStudy()` -- the OLS-family event-study route,
#    through `.event_study_etwfe_betwfe()`. A different route to the same
#    site, which is what the per-route half of the coverage predicate is
#    about.
# ------------------------------------------------------------------------------
test_that("site 2 fires both tiers via eventStudy on an etwfe fit (#474)", {
	.cpf474_skip()
	fit <- .cpf474("etwfe")

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermOLS", .CPF474_WARN_SCALE),
		function() eventStudy(fit)
	)
	fws <- .cpf474_floor_ws(r)
	expect_gt(length(fws), 0L)
	expect_s3_class(fws[[1]], "fetwfe_negative_variance_floored")
	expect_true(grepl(.CPF474_OLS, conditionMessage(fws[[1]]), fixed = TRUE))

	r <- .cpf474_run(
		.cpf474_neg("getSecondVarTermOLS", .CPF474_ERR_SCALE),
		function() eventStudy(fit)
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.CPF474_OLS, conditionMessage(r$err), fixed = TRUE))
	expect_false(grepl(.CPF474_REMEDY, conditionMessage(r$err), fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 4. (red) SITE 4 -- `.event_study_var2_fetwfe()` via `eventStudy()` on a
#    `fetwfe` fit, both tiers, plus the AGGREGATION BOUNDARY.
#
#    `.floor_psd_diag_core()`'s header advertises "at most ONE aggregated
#    condition per CALL" as the thing that stops K warnings burying the
#    signal. That contract holds -- and this site is called once per event
#    time inside `.event_study_fetwfe()`'s loop, so the aggregation does not
#    span the loop and a user sees one warning per event time that pools more
#    than one cohort. Asserted as the derived predicate (one warning per event
#    time with `n_cohorts > 1`) rather than as a hardcoded count, so the
#    fixture and the expectation cannot drift apart.
#
#    THE `|V_e| > 1` PRECONDITION. The site returns a hard structural `0` when
#    the valid cohort set has one member, so an injected negative never
#    reaches the quadratic form at a singleton event time; without this guard
#    a fixture whose every event time was a singleton would satisfy the block
#    vacuously. The repo's own idiom for the guard is
#    `expect_gt(length(V_e), 1L)` in `test-degenerate-jacobian-225.R`.
# ------------------------------------------------------------------------------
test_that("site 4 fires both tiers via eventStudy on a fetwfe fit (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	# ANTI-VACUITY: at least one event time pools more than one cohort, so the
	# early structural `return(0)` is not the only path exercised.
	clean <- eventStudy(fit)
	n_pooled <- sum(clean$n_cohorts > 1L)
	expect_gt(n_pooled, 0L)

	r <- .cpf474_run(
		.cpf474_neg(".event_study_var2_fetwfe", .CPF474_WARN_SCALE),
		function() eventStudy(fit)
	)
	fws <- .cpf474_floor_ws(r)
	# ONE aggregated condition per CALL, and one call per pooled event time.
	expect_length(fws, n_pooled)
	msg <- conditionMessage(fws[[1]])
	expect_s3_class(fws[[1]], "fetwfe_negative_variance_floored")
	expect_true(grepl(.CPF474_ES, msg, fixed = TRUE))
	expect_true(grepl(.CPF474_SUBJECT, msg, fixed = TRUE))

	r <- .cpf474_run(
		.cpf474_neg(".event_study_var2_fetwfe", .CPF474_ERR_SCALE),
		function() eventStudy(fit)
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.CPF474_ES, conditionMessage(r$err), fixed = TRUE))
	expect_false(grepl(.CPF474_REMEDY, conditionMessage(r$err), fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 5. (red) SITE 4 reached through the RENDERERS. `print()` / `summary()` /
#    `plot()` call `eventStudy()` internally, so they reach the scalar site as
#    well as `Sigma_2` -- which is where the new warning becomes repeated
#    output on an ordinary `print()`.
#
#    BOTH TIERS, and the warning tier under `options(warn = 2)` specifically,
#    because that is the only condition under which a muffle is visible. The
#    renderers wrap their `eventStudy()` call in `.event_study_quiet()`, whose
#    muffle is keyed to `fetwfe_highdim_postselection_band` ONLY. Were it ever
#    widened to a blanket warning muffle, site 4's WARNING would be swallowed
#    here while still arriving on the direct `eventStudy()` route of block 4 --
#    and a catastrophic-tier-only assertion could not notice, since an `error`
#    is not muffled by a warning muffle. The `warn = 2` half is what closes
#    that, and it is measured: replacing `.event_study_quiet()`'s class-keyed
#    handler with a blanket `warning = function(w) invokeRestart(...)` reddens
#    this block, block 8 and block 8b, while block 4's DIRECT `eventStudy()`
#    route stays green -- which is the discriminating half, since the muffle
#    only sits on the renderer route. Without the `warn = 2` assertions this
#    block stayed green under that mutation. Called BARE:
#    `suppressWarnings()`, and a muffling `withCallingHandlers()`, both defeat
#    the conversion.
# ------------------------------------------------------------------------------
test_that("site 4 reaches print/summary through eventStudy (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	for (render in list(
		function() capture.output(print(fit)),
		function() capture.output(print(summary(fit)))
	)) {
		r <- .cpf474_run(
			.cpf474_neg(".event_study_var2_fetwfe", .CPF474_ERR_SCALE),
			render
		)
		expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
		expect_true(grepl(.CPF474_ES, conditionMessage(r$err), fixed = TRUE))
		expect_false(grepl(
			.CPF474_REMEDY,
			conditionMessage(r$err),
			fixed = TRUE
		))
	}

	old <- options(warn = 2)
	on.exit(options(old), add = TRUE)
	for (render in list(
		function() capture.output(print(fit)),
		function() capture.output(print(summary(fit)))
	)) {
		err <- testthat::with_mocked_bindings(
			tryCatch(render(), error = function(e) e),
			.multinomial_cov = .cpf474_neg(
				".event_study_var2_fetwfe",
				.CPF474_WARN_SCALE
			),
			.package = "fetwfe"
		)
		expect_s3_class(err, "error")
		expect_true(grepl(.CPF474_ES, conditionMessage(err), fixed = TRUE))
		expect_true(grepl(
			"clipped to 0 (indices",
			conditionMessage(err),
			fixed = TRUE
		))
	}
})

# ------------------------------------------------------------------------------
# 6. (red) SITE 4's LIVE DEFECT, default branch. Before #474 this site had no
#    floor: a negative `var_2(e)` flowed into
#    `sqrt(.combine_att_variance(var_1_e, var_2_e, ...))`, giving `NaN` when
#    the sum went negative and a SILENTLY UNDERSTATED standard error when
#    `var_1(e)` kept it positive.
#
#    THE VALUE IS PINNED, not just its finiteness: `is.finite(se)` passes on
#    an SE that is fifty percent wrong. With `var_2(e)` floored to zero the
#    event-time SE is exactly `sqrt(var_1(e))`, which the reference mock
#    computes through the real code path by zeroing `Sigma_pi_hat` at the same
#    gate. That is a separate mechanism from the one under test, not a round
#    trip through it.
# ------------------------------------------------------------------------------
test_that("a negative var_2(e) is floored, not NaN'd, on the default branch (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	ref <- testthat::with_mocked_bindings(
		eventStudy(fit),
		.multinomial_cov = .cpf474_zero(".event_study_var2_fetwfe"),
		.package = "fetwfe"
	)
	# ANTI-VACUITY: `var_2(e)` is genuinely non-zero on this fixture, so the
	# reference is NOT the unperturbed answer and the comparison below can
	# tell a floor from a no-op.
	clean <- eventStudy(fit)
	expect_false(isTRUE(all.equal(clean$se, ref$se)))
	# ... and `var_1(e) > 0` everywhere, so `sqrt(var_1(e))` is a real number
	# to hit rather than a degenerate zero.
	expect_true(all(ref$se > 0))

	got <- testthat::with_mocked_bindings(
		suppressWarnings(eventStudy(fit)),
		.multinomial_cov = .cpf474_neg(
			".event_study_var2_fetwfe",
			.CPF474_WARN_SCALE
		),
		.package = "fetwfe"
	)
	expect_false(any(is.nan(got$se)))
	expect_equal(got$se, ref$se)
})

# ------------------------------------------------------------------------------
# 7. (red) SITE 4's LIVE DEFECT, CONSERVATIVE branch -- the one
#    `getSecondVarTermOLS()`'s own #127 comment named as the reason its floor
#    exists, and the one place in the family where the hazard was identified
#    and never closed. `.combine_att_variance()` returns
#    `var_1 + var_2 + 2 * sqrt(var_1 * var_2)` there, so ANY negative
#    `var_2(e)` gives `NaN` regardless of magnitude -- provided
#    `var_1(e) > 0`, which is asserted rather than assumed: at `var_1(e)`
#    exactly zero the product is `-0`, `sqrt()` returns `0`, and there is no
#    `NaN` for the floor to prevent.
# ------------------------------------------------------------------------------
test_that("a negative var_2(e) does not NaN the conservative branch (#474)", {
	.cpf474_skip()
	fit <- .cpf474("conservative")
	# The Cauchy-Schwarz branch is only taken when the fit did NOT use
	# independent counts.
	expect_false(isTRUE(fit$indep_counts_used))

	ref <- testthat::with_mocked_bindings(
		eventStudy(fit),
		.multinomial_cov = .cpf474_zero(".event_study_var2_fetwfe"),
		.package = "fetwfe"
	)
	clean <- eventStudy(fit)
	expect_false(isTRUE(all.equal(clean$se, ref$se)))
	# `att_var_1 > 0`: without this the assertion below is vacuous.
	expect_true(all(ref$se > 0))

	got <- testthat::with_mocked_bindings(
		suppressWarnings(eventStudy(fit)),
		.multinomial_cov = .cpf474_neg(
			".event_study_var2_fetwfe",
			.CPF474_WARN_SCALE
		),
		.package = "fetwfe"
	)
	expect_false(any(is.nan(got$se)))
	expect_equal(got$se, ref$se)
})

# ------------------------------------------------------------------------------
# 8. (red) SITE 1 -- `.assemble_joint_cov_var2()` on the RENDERER BAND route,
#    `print()` / `summary()` -> `eventStudy()` ->
#    `.event_study_simultaneous_bounds()` -> `.fit_band_for_family()`. This is
#    the route the classing exists for: the site sits INSIDE
#    `.fit_band_for_family()`'s `tryCatch(error = function(e) NULL)` region,
#    so an unclassed condition here is swallowed and the band degrades to the
#    pointwise one under a `[simultaneous 95% CI]` header -- the #433/#470
#    wrong-answer shape. This block, 8b and 8c are the ONLY guards on
#    `.fit_band_for_family()`'s two re-raises for THIS family, and the two
#    red sets are NOT the same set -- do not trim either as redundant.
#    Measured, deleting each re-raise against this tree:
#      * `for (w in pending_floor) warning(w)` -- blocks 8 and 8b go red;
#        8c stays green, because it drives only the catastrophic tier.
#      * `if (!is.null(fatal)) stop(fatal)` -- blocks 8 and 8c go red;
#        8b stays green, because it drives only the warning tier.
#    Block 9's direct `simultaneousCIs()` route is OUTSIDE the protected
#    region and guards neither; it stays green under both deletions.
#    (Until #474, `test-matrix-floor-conditions-470.R` was the SOLE guard on
#    both re-raises. It still reddens on both mutants, and so now does this
#    file, independently.)
#
#    THE GATE IS `.event_study_simultaneous_bounds`, per profile gotcha 12.8:
#    the fit-time band and the rendered band are different families over the
#    same sandwich, so gating on the fit-time path would measure the wrong
#    route.
# ------------------------------------------------------------------------------
test_that("site 1 survives the renderer band route in both tiers (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	for (render in list(
		function() capture.output(print(fit)),
		function() capture.output(print(summary(fit)))
	)) {
		# --- warning tier: captured inside the protected region, muffled,
		#     and RE-RAISED below it. Without the re-raise it never arrives.
		r <- .cpf474_run(
			.cpf474_neg(
				".event_study_simultaneous_bounds",
				.CPF474_WARN_SCALE
			),
			render
		)
		fws <- .cpf474_floor_ws(r)
		expect_gt(length(fws), 0L)
		hit <- Filter(
			function(w) {
				grepl(.CPF474_SIGMA2, conditionMessage(w), fixed = TRUE)
			},
			fws
		)
		expect_gt(length(hit), 0L)
		expect_s3_class(hit[[1]], "fetwfe_negative_variance_floored")
		# SUBJECT NOUN, on this route specifically. Class, site label and
		# remedy are all preserved by `.floor_variance_diag()`, which inherits
		# the same remedy from the core's default -- so without this line a
		# sibling-wrapper swap at site 1 is invisible HERE and caught only by
		# the direct-`simultaneousCIs()` block, i.e. on the route OUTSIDE the
		# protected region rather than the one the classing exists for. The
		# #482 review measured that asymmetry; the noun is what discriminates.
		expect_true(grepl(
			.CPF474_SUBJECT,
			conditionMessage(hit[[1]]),
			fixed = TRUE
		))

		# --- catastrophic tier: PROPAGATES rather than degrading to a NULL
		#     band, and carries the band remedy, which is TRUE here.
		r <- .cpf474_run(
			.cpf474_neg(".event_study_simultaneous_bounds", .CPF474_ERR_SCALE),
			render
		)
		expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
		emsg <- conditionMessage(r$err)
		expect_true(grepl(.CPF474_SIGMA2, emsg, fixed = TRUE))
		expect_true(grepl(.CPF474_SUBJECT_CAP, emsg, fixed = TRUE))
		# THE REMEDY IS PRESENT HERE and absent at the other three sites.
		# This half and its three negations are the pair that keeps the
		# false-remedy regression from returning silently.
		expect_true(grepl(.CPF474_REMEDY, emsg, fixed = TRUE))
	}
})

# ------------------------------------------------------------------------------
# 8b. (red) `warn = 2` on the renderer band route. The muffle inside
#     `.fit_band_for_family()`'s `tryCatch()` beats the warnings-as-errors
#     conversion, so without `for (w in pending_floor) warning(w)` below it
#     nothing errors and nothing warns -- the band silently degrades. Called
#     BARE: `suppressWarnings()`, and a muffling `withCallingHandlers()`,
#     both defeat the conversion this block is about.
# ------------------------------------------------------------------------------
test_that("warn = 2 keeps site 1 loud on the renderer route (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	old <- options(warn = 2)
	on.exit(options(old), add = TRUE)

	err <- testthat::with_mocked_bindings(
		tryCatch(capture.output(print(fit)), error = function(e) e),
		.multinomial_cov = .cpf474_neg(
			".event_study_simultaneous_bounds",
			.CPF474_WARN_SCALE
		),
		.package = "fetwfe"
	)
	expect_s3_class(err, "error")
	expect_true(grepl(.CPF474_SIGMA2, conditionMessage(err), fixed = TRUE))
	expect_true(grepl(
		"clipped to 0 (indices",
		conditionMessage(err),
		fixed = TRUE
	))
})

# ------------------------------------------------------------------------------
# 8c. `plot()` is the third renderer, split out so a machine without
#     `ggplot2` -- a Suggest -- skips only this and never a load-bearing
#     block. Same split, same reason, as
#     `test-matrix-floor-conditions-470.R`'s block 4c.
# ------------------------------------------------------------------------------
test_that("plot() stays loud on the cohort-probability breakdown (#474)", {
	.cpf474_skip()
	skip_if_not_installed("ggplot2")
	fit <- .cpf474("fetwfe")

	r <- .cpf474_run(
		.cpf474_neg(".event_study_simultaneous_bounds", .CPF474_ERR_SCALE),
		function() print(plot(fit, type = "event"))
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.CPF474_SIGMA2, conditionMessage(r$err), fixed = TRUE))
	expect_true(grepl(.CPF474_REMEDY, conditionMessage(r$err), fixed = TRUE))

	r <- .cpf474_run(
		.cpf474_neg(".event_study_var2_fetwfe", .CPF474_ERR_SCALE),
		function() print(plot(fit, type = "event"))
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	expect_true(grepl(.CPF474_ES, conditionMessage(r$err), fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 9. (red) SITE 1's second live route: a DIRECT `simultaneousCIs()` call on a
#    family other than `"cohort"`. Outside the protected region, so the
#    condition reaches the caller directly rather than through the deferral.
#    The band still comes back at the warning tier.
# ------------------------------------------------------------------------------
test_that("site 1 fires on a direct simultaneousCIs call (#474)", {
	.cpf474_skip()
	fit <- .cpf474("fetwfe")

	r <- .cpf474_run(
		.cpf474_neg(".simultaneous_cis_impl", .CPF474_WARN_SCALE_DIRECT),
		function() {
			suppressMessages(simultaneousCIs(fit, family = "event_study"))
		}
	)
	expect_null(r$err)
	fws <- .cpf474_floor_ws(r)
	expect_length(fws, 1L)
	msg <- conditionMessage(fws[[1]])
	expect_true(grepl(.CPF474_SIGMA2, msg, fixed = TRUE))
	expect_true(grepl(.CPF474_SUBJECT, msg, fixed = TRUE))

	r <- .cpf474_run(
		.cpf474_neg(".simultaneous_cis_impl", .CPF474_ERR_SCALE),
		function() {
			suppressMessages(simultaneousCIs(fit, family = "event_study"))
		}
	)
	expect_s3_class(r$err, "fetwfe_negative_variance_catastrophic")
	emsg <- conditionMessage(r$err)
	expect_true(grepl(.CPF474_SIGMA2, emsg, fixed = TRUE))
	expect_true(grepl(.CPF474_REMEDY, emsg, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# 10. SITE 1's THIRD route is STRUCTURALLY DEAD, and this block is what makes
#     that a measured fact rather than an assumption.
#
#     The fit-time band -- `.finalize_ci_type()` ->
#     `.apply_simultaneous_catt_band()` -> `.fit_band_for_family(x, "cohort",
#     ...)` -- is hardcoded to the COHORT family (profile gotcha 12.8), and
#     every cohort-family `J_k` is identically zero, so `diag(Sigma_2)` is
#     exactly zero there and no perturbation of `Sigma_pi_hat` can make the
#     floor bind. Measured: a `-1e6 * I` cohort covariance raises nothing on
#     that route while firing on both live ones.
#
#     THE ASSERTION READS THE UNFLOORED `J_list`, NOT THE FLOORED DIAGONAL.
#     Written as `diag(.assemble_joint_cov_var2(...))`, the floor this issue
#     adds would sit between the observation and the check: a tree whose
#     cohort Jacobians stopped being identically zero AND went negative would
#     read back as exactly zero, so the assertion would be blind by
#     construction to the very change it is meant to notice, erased by the
#     transform #474 introduced. The floored-diagonal form stands below as
#     CORROBORATION, never as the assertion.
#
#     An `expect_no_condition()` on the fit-time route is not the alternative:
#     it is a vacuous fixture that passes on a tree where this site's floor
#     has been deleted outright.
# ------------------------------------------------------------------------------
test_that("the cohort family's Sigma_2 Jacobians are identically zero (#474)", {
	d_inv <- matrix(seq_len(20 * 4) / 7, nrow = 20, ncol = 4)
	args <- list(
		K = 3L,
		G = 3L,
		T = 6L,
		cohort_offsets_int = c(2L, 3L, 4L),
		first_inds = c(1L, 7L, 13L),
		cohort_probs_overall = c(0.3, 0.3, 0.4),
		d_inv_treat_sel = d_inv
	)

	j_cohort <- do.call(
		fetwfe:::.build_j_list_for_family,
		c(list(family = "cohort"), args)
	)
	expect_length(j_cohort, args$K)
	expect_true(all(vapply(
		j_cohort,
		function(m) all(m == 0),
		logical(1)
	)))

	# ANTI-VACUITY. The builder is capable of returning a non-zero Jacobian on
	# these very inputs, so "all zero" above is a property of the COHORT
	# family and not of the fixture.
	j_es <- do.call(
		fetwfe:::.build_j_list_for_family,
		c(list(family = "event_study"), utils::modifyList(args, list(K = 5L)))
	)
	expect_true(any(vapply(j_es, function(m) any(m != 0), logical(1))))

	# CORROBORATION ONLY (see the block comment): the assembled diagonal is
	# exactly zero, which follows from the Jacobians above and is therefore
	# not independent evidence.
	sigma_pi <- fetwfe:::.multinomial_cov(args$cohort_probs_overall)
	diag_2 <- diag(fetwfe:::.assemble_joint_cov_var2(
		J_list = j_cohort,
		theta_sel = seq_len(ncol(d_inv)) / 3,
		Sigma_pi_hat = sigma_pi,
		N = 200L,
		T = args$T
	))
	expect_true(all(diag_2 == 0))
})
