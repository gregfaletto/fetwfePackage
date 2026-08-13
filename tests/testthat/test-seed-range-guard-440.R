library(testthat)
library(fetwfe)

# ==============================================================================
# Out-of-range seed guard at every door (#440).
#
# `set.seed()` accepts only values representable as a 32-bit signed integer.
# Anything larger in magnitude is coerced to NA, so base R emits a coercion
# WARNING and then an opaque ERROR that never names the argument or the limit.
# This file pins the package's own message at all four guarded sites, plus the
# `match.arg()` deletion in the four *WithSimulatedData() wrappers.
#
# TWO ASSERTION RULES, both learned the hard way:
#
#   1. Every error assertion pins the literal "integer.max" with fixed = TRUE.
#      A looser pattern such as "seed" is green BEFORE the change too, because
#      base R's own message -- "supplied seed is not a valid integer" --
#      contains the word "seed". It does NOT contain "integer.max".
#
#   2. Every "no coercion warning" assertion wraps the call in
#      try(..., silent = TRUE). Without the wrapper the call errors in both
#      states, expect_no_warning() re-raises that error out of the expectation,
#      and the test reports an ERROR rather than a pass -- red before AND after.
#
# Groups B, C, and G are unchanged-behavior anchors and pass in both states;
# Groups A, D, E, F, H, and I flip.
# ==============================================================================

# ------------------------------------------------------------------------------
# File-scope fixture. Built once here, above Group A, because Groups G and I
# both need it: Group G's fits live inside test_that() blocks and do not escape
# them, so "reuse Group G's fit" is not executable.
# ------------------------------------------------------------------------------
coefs <- genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 1, seed = 1)
sim <- simulateData(
	coefs,
	N = 40,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 2
)
fit <- suppressWarnings(suppressMessages(fetwfeWithSimulatedData(sim)))

# ==============================================================================
# Group A -- .apply_seed() rejects out-of-range seeds with the package message
# ==============================================================================

test_that("Group A: .apply_seed() rejects out-of-range seeds", {
	expect_error(fetwfe:::.apply_seed(3e9), "integer.max", fixed = TRUE)
	expect_error(fetwfe:::.apply_seed(-3e9), "integer.max", fixed = TRUE)
	# R's integer MINIMUM is -2147483647, not -2147483648, so the negative side
	# overflows one value earlier than a two's-complement reader expects. That
	# is why the predicate uses abs() rather than a bare comparison.
	expect_error(
		fetwfe:::.apply_seed(-2147483648),
		"integer.max",
		fixed = TRUE
	)
	expect_error(fetwfe:::.apply_seed(Inf), "integer.max", fixed = TRUE)
	expect_error(fetwfe:::.apply_seed(-Inf), "integer.max", fixed = TRUE)
})

test_that("Group A: the base-R coercion warning no longer fires", {
	# try() is load-bearing here -- see rule 2 in the header.
	expect_no_warning(try(fetwfe:::.apply_seed(3e9), silent = TRUE))
})

# ==============================================================================
# Group B -- the boundary is strict (`>`), and that is load-bearing
# ==============================================================================

test_that("Group B: exactly .Machine$integer.max is still accepted", {
	# The strict `>` boundary is load-bearing: getBetaCV() and
	# .fit_q1_nuisance() (R/bridge_selection.R) and .cv_lambda_node()
	# (R/riesz_lasso.R) all clip a computed default seed to exactly
	# .Machine$integer.max and feed it through .apply_seed(). A `>=` boundary
	# would turn ordinary large-N*T fits into errors.
	expect_no_error(fetwfe:::.apply_seed(.Machine$integer.max))
	expect_no_warning(fetwfe:::.apply_seed(.Machine$integer.max))
	expect_no_error(fetwfe:::.apply_seed(-.Machine$integer.max))
	expect_no_warning(fetwfe:::.apply_seed(-.Machine$integer.max))
})

test_that("Group B: the boundary seed is actually applied, not skipped", {
	# Deliberately NOT expect_equal(.with_preserved_rng(int.max, 42), 42):
	# that asserts only that a literal survives a pass-through, and its only
	# failure mode is "the call errored", so a silently-skipped seed would
	# still pass. Comparing draws proves set.seed() ran.
	got <- fetwfe:::.with_preserved_rng(.Machine$integer.max, runif(1))
	set.seed(.Machine$integer.max)
	expect_identical(got, runif(1))
})

# ==============================================================================
# Group C -- the guard's placement did not break the ambient-RNG contracts
# ==============================================================================

test_that("Group C: NULL / NA / NaN seeds still mean 'ambient RNG'", {
	# A guard hoisted above the is.na() early return would make these throw
	# "missing value where TRUE/FALSE needed", because abs(NA) >
	# .Machine$integer.max is NA and `if (NA)` errors.
	expect_no_error(fetwfe:::.apply_seed(NULL))
	expect_no_error(fetwfe:::.apply_seed(NA))
	expect_no_error(fetwfe:::.apply_seed(NA_real_))
	expect_no_error(fetwfe:::.apply_seed(NA_integer_))
	# NaN belongs in this group, not Group A: is.na(NaN) is TRUE, so NaN is
	# caught by the is.na() early return and stays a supported ambient value.
	# The roxygen must not claim the function "rejects non-finite seeds".
	expect_no_error(fetwfe:::.apply_seed(NaN))
})

test_that("Group C: an in-range seed still reaches set.seed()", {
	fetwfe:::.apply_seed(1)
	got <- runif(1)
	set.seed(1)
	expect_identical(got, runif(1))
})

# ==============================================================================
# Group D -- .validate_boot_args() is the bootstrap path's front door
# ==============================================================================

test_that("Group D: .validate_boot_args() range-checks `seed`", {
	expect_error(
		fetwfe:::.validate_boot_args(100, 3e9, "simultaneousCIs"),
		"integer.max",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.validate_boot_args(100, 3e9, "simultaneousCIs"),
		"simultaneousCIs()",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.validate_boot_args(100, -3e9, "simultaneousCIs"),
		"integer.max",
		fixed = TRUE
	)
	# Inf slips through the pre-existing branches entirely: Inf != round(Inf)
	# is FALSE and is.na(Inf) is FALSE. The !is.finite() term catches it.
	expect_error(
		fetwfe:::.validate_boot_args(100, Inf, "simultaneousCIs"),
		"integer.max",
		fixed = TRUE
	)
	# Both flavours of the message prefix are witnessed.
	expect_error(
		fetwfe:::.validate_boot_args(100, 3e9, "debiasedATT"),
		"debiasedATT()",
		fixed = TRUE
	)
})

test_that("Group D: .validate_boot_args() range-checks `B`", {
	# Today this returns NA_integer_ after a coercion warning, which surfaces
	# several frames later as "vector size cannot be NA".
	expect_error(
		fetwfe:::.validate_boot_args(3e9, NULL, "simultaneousCIs"),
		"integer.max",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.validate_boot_args(3e9, NULL, "simultaneousCIs"),
		"simultaneousCIs()",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.validate_boot_args(Inf, NULL, "simultaneousCIs"),
		"integer.max",
		fixed = TRUE
	)
	# `B = 3e9` is the one path here where the coercion warning genuinely fires
	# today. Deliberately NOT asserted for `seed = Inf`, which is already
	# silent on main and would therefore be green in both states.
	expect_no_warning(
		try(
			fetwfe:::.validate_boot_args(3e9, NULL, "simultaneousCIs"),
			silent = TRUE
		)
	)
})

test_that("Group D: unchanged-behavior anchors", {
	expect_identical(
		fetwfe:::.validate_boot_args(100, .Machine$integer.max, "x"),
		100L
	)
	expect_no_warning(
		fetwfe:::.validate_boot_args(100, .Machine$integer.max, "x")
	)
	expect_identical(fetwfe:::.validate_boot_args(100, NULL, "x"), 100L)
	# This PR is a RANGE fix only. The non-integral contract is unchanged:
	# .validate_boot_args() still rejects 1.5 while .apply_seed() still
	# truncates it silently.
	expect_error(
		fetwfe:::.validate_boot_args(100, 1.5, "x"),
		"must be NULL or a single integer",
		fixed = TRUE
	)
	expect_no_error(fetwfe:::.apply_seed(1.5))
	expect_error(
		fetwfe:::.validate_boot_args(0, NULL, "x"),
		"must be a single positive integer",
		fixed = TRUE
	)
})

# ==============================================================================
# Group E -- cv_seed fails at input validation, for both estimators
# ==============================================================================

test_that("Group E: fetwfe() rejects an out-of-range cv_seed at entry", {
	df <- generate_panel_data()
	expect_error(
		fetwfe(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			cv_seed = 3e9
		),
		"integer.max",
		fixed = TRUE
	)
	expect_error(
		fetwfe(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			cv_seed = 3e9
		),
		"Invalid inputs:",
		fixed = TRUE
	)
})

test_that("Group E: betwfe() rejects an out-of-range cv_seed at entry", {
	# Both estimators route through the same checkFetwfeInputs(); the
	# shared-validator claim deserves its own witness.
	df <- generate_panel_data()
	expect_error(
		betwfe(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			cv_seed = 3e9
		),
		"integer.max",
		fixed = TRUE
	)
})

test_that("Group E: the cv_seed violation joins the collect-all list", {
	# Two bad arguments, one message. This is what proves the new check landed
	# in the `violations` collector rather than as a bare stop().
	df <- generate_panel_data()
	msg <- tryCatch(
		fetwfe(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			cv_seed = 3e9,
			cv_folds = 1L
		),
		error = conditionMessage
	)
	expect_true(grepl("integer.max", msg, fixed = TRUE))
	expect_true(grepl("cv_folds must be an integer >= 2", msg, fixed = TRUE))
})

test_that("Group E: cv_seed = .Machine$integer.max is not rejected", {
	# Positive control, called directly rather than through a fit: it costs
	# ~0 s and tests exactly the claim, whereas a full fit's non-rejection
	# could be confounded by any other failure.
	df <- generate_panel_data()
	expect_no_error(
		fetwfe:::checkFetwfeInputs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2"),
			cv_seed = .Machine$integer.max
		)
	)
})

# ==============================================================================
# Group F -- the end-to-end data-generation seed doors
# ==============================================================================

test_that("Group F: simulateData() and simulateDataCore() reject 3e9", {
	expect_error(
		simulateData(
			coefs,
			N = 40,
			sig_eps_sq = 0.5,
			sig_eps_c_sq = 0.5,
			seed = 3e9
		),
		"integer.max",
		fixed = TRUE
	)
	expect_no_warning(
		try(
			simulateData(
				coefs,
				N = 40,
				sig_eps_sq = 0.5,
				sig_eps_c_sq = 0.5,
				seed = 3e9
			),
			silent = TRUE
		)
	)
	# simulateDataCore() is the function that literally calls .apply_seed(),
	# and it is exported in its own right -- so it is a door, not a detail.
	expect_error(
		simulateDataCore(
			N = 40,
			T = 4,
			G = 3,
			d = 2,
			sig_eps_sq = 0.5,
			sig_eps_c_sq = 0.5,
			beta = coefs$beta,
			seed = 3e9
		),
		"integer.max",
		fixed = TRUE
	)
})

test_that("Group F: genCoefs() and genCoefsCore() reject 3e9", {
	expect_error(
		genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 1, seed = 3e9),
		"integer.max",
		fixed = TRUE
	)
	expect_no_warning(
		try(
			genCoefs(
				G = 3,
				T = 4,
				d = 2,
				density = 0.5,
				eff_size = 1,
				seed = 3e9
			),
			silent = TRUE
		)
	)
	expect_error(
		genCoefsCore(
			G = 3,
			T = 4,
			d = 2,
			density = 0.5,
			eff_size = 1,
			seed = 3e9
		),
		"integer.max",
		fixed = TRUE
	)
})

# ==============================================================================
# Group G -- deleting the wrappers' match.arg() preserves partial matching
# ==============================================================================

test_that("Group G: fetwfeWithSimulatedData() still partial-matches", {
	res <- suppressWarnings(suppressMessages(fetwfeWithSimulatedData(
		sim,
		se_type = "conserv",
		ci_type = "point",
		fusion_structure = "event"
	)))
	expect_identical(res$se_type, "conservative")
	expect_identical(res$ci_type, "pointwise")
	expect_identical(res$fusion_structure, "event_study")
})

test_that("Group G: betwfeWithSimulatedData() still partial-matches", {
	res <- suppressWarnings(suppressMessages(betwfeWithSimulatedData(
		sim,
		se_type = "conserv",
		ci_type = "point"
	)))
	expect_identical(res$se_type, "conservative")
	expect_identical(res$ci_type, "pointwise")
})

test_that("Group G: etwfeWithSimulatedData() still partial-matches", {
	res <- suppressWarnings(suppressMessages(etwfeWithSimulatedData(
		sim,
		se_type = "conserv",
		ci_type = "point"
	)))
	expect_identical(res$se_type, "conservative")
	expect_identical(res$ci_type, "pointwise")
})

test_that("Group G: twfeCovsWithSimulatedData() still partial-matches", {
	res <- suppressWarnings(suppressMessages(twfeCovsWithSimulatedData(
		sim,
		se_type = "conserv",
		ci_type = "point"
	)))
	expect_identical(res$se_type, "conservative")
	expect_identical(res$ci_type, "pointwise")
})

test_that("Group G: a bad se_type still errors, from the core's match.arg()", {
	# Cheap: the core's match.arg() fires before any fitting.
	expect_error(
		fetwfeWithSimulatedData(sim, se_type = "bogus"),
		"should be one of"
	)
	expect_error(
		betwfeWithSimulatedData(sim, se_type = "bogus"),
		"should be one of"
	)
	expect_error(
		etwfeWithSimulatedData(sim, se_type = "bogus"),
		"should be one of"
	)
	expect_error(
		twfeCovsWithSimulatedData(sim, se_type = "bogus"),
		"should be one of"
	)
})

# ==============================================================================
# Group H -- the documented error-ordering change
# ==============================================================================

test_that("Group H: a bad simulated_obj now outranks a bad se_type", {
	# Before the match.arg() deletion the wrapper's own match.arg() fired first
	# and reported the bad se_type; now .unpack_simulated_obj() fires first.
	# A real, disclosed behavior change.
	expect_error(
		fetwfeWithSimulatedData(list(a = 1), se_type = "bogus"),
		"must be an object of class 'FETWFE_simulated'",
		fixed = TRUE
	)
})

# ==============================================================================
# Group I -- the public bootstrap doors
# ==============================================================================

test_that("Group I: simultaneousCIs() reports the range error itself", {
	expect_error(
		simultaneousCIs(fit, method = "bootstrap", seed = 3e9),
		"integer.max",
		fixed = TRUE
	)
	# Today this is `vector size cannot be NA`, thrown several frames deeper
	# inside the multiplier draw.
	expect_error(
		simultaneousCIs(fit, method = "bootstrap", B = 3e9),
		"integer.max",
		fixed = TRUE
	)
})

test_that("Group I: debiasedATT() reports the range error itself", {
	expect_error(
		debiasedATT(fit, method = "bootstrap", B = 3e9),
		"integer.max",
		fixed = TRUE
	)
})
