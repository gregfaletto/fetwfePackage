library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for the cluster-floor family -- `.floor_cluster_quad()` (issue #139,
# version 1.11.2) and its vectorized siblings `.floor_cluster_quad_diag()` /
# `.floor_variance_diag()` (issue #470) -- and the guardrail keeping every
# cluster-sandwich quadratic-form site routed through one of them (issue #463).
#
# The family layers a two-tier diagnostic on top of the pre-existing
# `max(q, 0)` / `pmax(diag(.), 0)` floor at each cluster-sandwich
# quadratic-form site. The forms are PSD in exact arithmetic, so any negative
# value is either FP-noise (silently clipped to 0) or a bug signal
# (warn / stop).
#
# The guardrail below states a PREDICATE, not a count and not a file list --
# both of those went stale here once already. It walks the AST of every
# function in the package NAMESPACE and asserts that every cluster-sandwich
# quadratic form (a `%*%` / `crossprod()` / `tcrossprod()` node whose operands
# reach the cluster-robust sandwich) sits lexically inside a call to one of
# the RECOGNIZED FLOOR FUNCTIONS, `.clf_floor_fns` below. Reading the
# namespace instead of `R/*.R` as text is what makes it run under
# `R CMD check`, where the tests execute against the *installed* package and
# there is no `R/` directory at all -- the whole block used to skip there,
# which is every CI job and every CRAN machine (#463).
#
# One block here is not about the floor sites at all. A11 / A12 walk the
# namespace CALL GRAPH rather than the AST, to pin the reachability claim
# `R/cluster_floor.R`'s header rests on: that no `.floor_cluster_quad()` call
# site sits inside `.fit_band_for_family()`'s protected region, which is what
# lets the scalar helper's conditions stay unclassed. They live here because
# that claim is the scalar floor's, and because the caller set they read is the
# one A10 pins.
#
# On the per-site labels asserted below: `getTeResultsOLS/att_var_1` and
# `getTeResults2/att_var_1` are two LABELS at one SITE. #344 merged those two
# functions' floors into the shared `.compute_att_var1()`, and each caller
# passes its own label through the `label` formal -- which is also why the
# label assertion has no power at that site, and why the floor-call inventory
# (A4) exists.
#
# THE TWO RECOGNIZED FLOOR FUNCTIONS ARE NOT INTERCHANGEABLE, which is why A10
# pins each site to one of them by name. `.floor_cluster_quad()` takes a SCALAR
# quadratic form; `.floor_cluster_quad_diag()` takes the K x K MATRIX and floors
# its diagonal. Handing either the other's argument is a silent revert rather
# than a type error, because each simply declines to act on what it does not
# recognize -- and until #476 the matrix helper declined by returning its
# argument untouched, so a scalar site swapped to it lost both the #139
# diagnostic and the `max(q, 0)` floor underneath it with the suite green. The
# helper now `stop()`s on a non-matrix, and A10 pins the pairing lexically; the
# defect needed both halves, since either alone still leaves a green tree in
# some spelling.
#
# THE RULE THESE ASSERTIONS ENFORCE, stated once, positively: a cluster-sandwich
# quadratic form and its floor are written TOGETHER -- one expression, one
# function, the form appearing literally as the floor call's argument. Not
# because that spelling is inherently better, but because it is the only
# arrangement a lexical check can verify. Split them across a statement, a
# helper, or a `...` forwarder and no static check can tell a correct
# refactor from a silent revert: a helper that returns the form for its caller
# to floor, and a decoy floor call in front of a helper that floors nothing,
# are the SAME SHAPE. Measured, on this tree, in both directions.
#
# So the code below rejects both, and the convention is the price of the
# guarantee. If a refactor needs to break it, that is a decision for a human,
# which is what a red test is for.
#
# WHAT GOES RED ON CORRECT CODE. All measured, none hypothetical:
#   * a form hoisted into its own statement and floored on the next line;
#   * a form computed by a callee and floored by the caller;
#   * a floor call reached through `...` forwarding;
#   * a new floored site, or new sandwich-touching code, added as a NEW
#     FUNCTION -- because the expected sets are exact rather than subsets,
#     which is what stops a floor added tomorrow being removed invisibly the
#     day after. Added INSIDE a function already in those sets, it is
#     invisible: they pin function names, not call counts;
#   * a SECOND floor site added inside a function already in the sets, calling
#     the OTHER member of `.clf_floor_fns` -- A4 is blind to it (function names,
#     not call counts) but A10 pairs owner with callee, so a legitimately
#     matrix-valued new floor inside, say, `getCohortATTsFinal()` fires A10
#     alone. Extend A10's literal, deliberately, rather than relaxing it;
#   * a wrapper between the floor call and its form that is not in
#     `.clf_transparent` -- extend that vector rather than deleting the
#     assertion.
#
# WHAT GETS THROUGH. A lexical guardrail has a boundary, and this list is the
# boundary. It is written out because the block this one replaced advertised
# coverage it did not have -- "the five sites", and an `any(grepl())` check
# sold as catching duplicate labels, which it could not -- and that is the
# failure #463 IS. Every entry below was measured green on a mutant tree.
#   * An ADDITIONAL unfloored form inside an already-listed function, computed
#     by a helper whose formal is named `S` rather than `sandwich_full`. The
#     function's own correct floor call satisfies everything else.
#   * A decoy floor call handed a REAL form on a degenerate input, or one
#     whose result is discarded. The checks inspect a floor call's argument,
#     never its result or its runtime value.
#   * An alias built by anything other than symbol-to-symbol assignment inside
#     an already-listed function: a submatrix slice, a list element,
#     `assign()`.
#   * Every floor site of `.floor_variance_diag()`. It is deliberately absent
#     from `.clf_floor_fns` below, so no assertion here sees its call sites at
#     all (A5b, a text grep over every body, is the one exception). Its
#     arguments are MODEL-BASED covariance diagonals rather than sandwich
#     forms, so adding it would make A9 -- "every floor call is handed a real
#     quadratic form" -- report `.simultaneous_cis_impl` as a decoy. Measured:
#     adding it fires A9 naming exactly that function. Refining A9 by owning
#     function does not escape it either, since `.simultaneous_cis_impl`
#     mentions `sandwich_full` and is in A6's inventory (#470).
#   * `abs()` or any sanitizer applied to an OPERAND inside the product --
#     the path from floor to form is checked, the form's operands are not.
#   * Condition suppression more than ONE frame above a floor call. One frame
#     is checked; a call graph would be needed for more. #470 created the
#     first LIVE instance of this blind spot: `.assemble_joint_cov_var1()`'s
#     conditions are captured and muffled two frames up, inside
#     `.fit_band_for_family()`, and what makes that correct rather than a
#     silent #139 revert is the pair of re-raises below its `tryCatch()`.
#     Delete either and the whole #470 diagnostic goes silent on every
#     internal route while A1-A9 stay green. `test-matrix-floor-conditions-470.R`
#     is the only thing guarding them.
#
# None of these is a regression. The block this replaces was blind to all of
# them AND to whole files, and it ran on no automated machine. The single
# thing it did better -- reading formal defaults, which it got for free by
# scanning raw text -- is closed: `.ns_code_exprs()` scans defaults alongside
# bodies, so a bare `max(... sandwich_full ..., 0)` in a default is now
# rejected everywhere.
#
# Do not add a claim to either list without a mutant behind it. Three earlier
# drafts of this header each stated a coverage claim one case too wide, and
# each was caught only by someone running the mutation rather than reading the
# sentence.
# ------------------------------------------------------------------------------

# --- Unit tests on the helper -------------------------------------------------

test_that(".floor_cluster_quad returns positive values unchanged", {
	expect_equal(fetwfe:::.floor_cluster_quad(0.5, "test"), 0.5)
	expect_equal(fetwfe:::.floor_cluster_quad(1e-12, "test"), 1e-12)
	expect_equal(fetwfe:::.floor_cluster_quad(100, "test"), 100)
})

test_that(".floor_cluster_quad returns 0 unchanged", {
	expect_equal(fetwfe:::.floor_cluster_quad(0, "test"), 0)
})

test_that(".floor_cluster_quad silently floors FP-noise negatives to 0", {
	# Magnitudes inside the FP-noise band ([-1e-10, 0]) clip without warning.
	expect_silent(out1 <- fetwfe:::.floor_cluster_quad(-1e-15, "test"))
	expect_equal(out1, 0)
	expect_silent(out2 <- fetwfe:::.floor_cluster_quad(-1e-11, "test"))
	expect_equal(out2, 0)
	# Exactly at the warning threshold is still silent (strict `<`).
	expect_silent(out3 <- fetwfe:::.floor_cluster_quad(-1e-10, "test"))
	expect_equal(out3, 0)
})

test_that(".floor_cluster_quad warns on warning-range negatives", {
	expect_warning(
		out <- fetwfe:::.floor_cluster_quad(-1e-8, "test_site"),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out, 0)
	# Site name in the message.
	expect_warning(
		fetwfe:::.floor_cluster_quad(-1e-8, "test_site"),
		"test_site"
	)
	# Boundary: just-inside the error threshold still warns, doesn't stop.
	expect_warning(
		out_boundary <- fetwfe:::.floor_cluster_quad(-0.999, "test"),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out_boundary, 0)
})

test_that(".floor_cluster_quad errors on catastrophic negatives", {
	expect_error(
		fetwfe:::.floor_cluster_quad(-2, "test_site"),
		"catastrophically negative"
	)
	expect_error(
		fetwfe:::.floor_cluster_quad(-1.5, "another_site"),
		"another_site"
	)
	expect_error(
		fetwfe:::.floor_cluster_quad(-1000, "test"),
		"file an issue"
	)
})

test_that(".floor_cluster_quad pass-through on NA / multi-element / non-numeric", {
	# NA passes through.
	expect_silent(out_na <- fetwfe:::.floor_cluster_quad(NA_real_, "test"))
	expect_true(is.na(out_na))
	# Multi-element vector passes through unchanged.
	expect_silent(out_vec <- fetwfe:::.floor_cluster_quad(c(-2, 1), "test"))
	expect_equal(out_vec, c(-2, 1))
	# Non-numeric (character) passes through.
	expect_silent(
		out_chr <- fetwfe:::.floor_cluster_quad("not_numeric", "test")
	)
	expect_equal(out_chr, "not_numeric")
	# Length-0 numeric passes through.
	expect_silent(out_e <- fetwfe:::.floor_cluster_quad(numeric(0), "test"))
	expect_equal(out_e, numeric(0))
})

test_that(".floor_cluster_quad respects custom thresholds", {
	# Tighter warn threshold -> a value that would normally pass silently warns.
	expect_warning(
		out <- fetwfe:::.floor_cluster_quad(
			-1e-12,
			"test",
			warn_threshold = -1e-13
		),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out, 0)
	# Looser error threshold -> a value that would normally warn errors.
	expect_error(
		fetwfe:::.floor_cluster_quad(
			-0.5,
			"test",
			err_threshold = -0.1
		),
		"catastrophically negative"
	)
})

# --- Unit tests on the vectorized family (#470) -------------------------------
# `.floor_cluster_quad_diag()` (matrix in / matrix out) and
# `.floor_variance_diag()` (vector in / vector out) over the shared
# `.floor_psd_diag_core()`. They mirror the scalar block above, plus the three
# properties the scalar helper has no analogue for: ONE aggregated condition
# per call however many entries offend, untouched off-diagonals, and a classed
# condition on both tiers.
#
# Every message assertion here is `fixed = TRUE`. The rendered messages carry
# `(`, `)`, `.` and `-`, and a parenthesised literal used as a REGEX matches
# text with the parentheses stripped out entirely, while an unbalanced one is a
# hard `invalid regular expression` error. `expect_error()` / `expect_warning()`
# pass the pattern to `grepl()` under edition 2, so both reach the assertion.

# Collect every condition of one class raised by `expr`, without muffling
# anything else, and return them. An observe-only handler on purpose: a
# muffling one would stand in for the behavior under test.
.clf_collect <- function(expr, class) {
	got <- list()
	withCallingHandlers(
		force(expr),
		condition = function(cnd) {
			if (inherits(cnd, class)) {
				got[[length(got) + 1L]] <<- cnd
				if (inherits(cnd, "warning")) {
					invokeRestart("muffleWarning")
				}
			}
		}
	)
	got
}

test_that(".floor_variance_diag leaves non-negative diagonals alone", {
	expect_silent(out <- fetwfe:::.floor_variance_diag(c(0.5, 1e-12, 100), "t"))
	expect_equal(out, c(0.5, 1e-12, 100))
	# Exact zero is unchanged, and so is the boundary of the FP-noise band.
	expect_silent(out0 <- fetwfe:::.floor_variance_diag(c(0, 0, 0), "t"))
	expect_equal(out0, c(0, 0, 0))
})

test_that(".floor_variance_diag silently floors FP-noise negatives", {
	expect_silent(
		out <- fetwfe:::.floor_variance_diag(c(1, -1e-15, -1e-11, -1e-10), "t")
	)
	# Strict `<` at the threshold, matching the scalar helper.
	expect_equal(out, c(1, 0, 0, 0))
})

test_that(".floor_variance_diag warns once and names every offending index", {
	ws <- .clf_collect(
		out <- fetwfe:::.floor_variance_diag(
			c(1, -1e-8, 2, -3e-8, 3),
			"simultaneous_cis_impl/Sigma"
		),
		"fetwfe_negative_variance_floored"
	)
	# EXACTLY ONE condition, not one per offending entry.
	expect_length(ws, 1L)
	msg <- conditionMessage(ws[[1]])
	expect_true(grepl("(indices 2, 4;", msg, fixed = TRUE))
	expect_true(grepl("2 of 5 entries clipped to 0", msg, fixed = TRUE))
	expect_true(grepl("most negative -3e-08", msg, fixed = TRUE))
	expect_true(grepl("site 'simultaneous_cis_impl/Sigma'", msg, fixed = TRUE))
	# The subject is "variance", NOT the sandwich phrase: these diagonals are
	# model-based, so the #139 wording would be a false statement about them.
	expect_true(grepl("Negative variance", msg, fixed = TRUE))
	expect_false(
		grepl("cluster-sandwich quadratic form", msg, fixed = TRUE)
	)
	expect_equal(out, c(1, 0, 2, 0, 3))
})

test_that(".floor_variance_diag errors on catastrophic negatives", {
	e <- tryCatch(
		fetwfe:::.floor_variance_diag(c(1, -2.4, -1e-8), "some_site"),
		error = function(e) e
	)
	expect_s3_class(e, "fetwfe_negative_variance_catastrophic")
	msg <- conditionMessage(e)
	# The error tier reports the ERROR-range indices only, and its subject is
	# sentence-cased because it opens the message.
	expect_true(grepl("Variance on the covariance diagonal", msg, fixed = TRUE))
	expect_true(grepl("1 of 3 entries below -1 (indices", msg, fixed = TRUE))
	expect_true(grepl("(indices 2;", msg, fixed = TRUE))
	expect_true(grepl("site 'some_site'", msg, fixed = TRUE))
	expect_true(grepl("file an issue", msg, fixed = TRUE))
})

test_that(".floor_variance_diag propagates NA and passes non-numerics through", {
	# `pmax()` propagates NA; `which(v < threshold)` drops it, so NA is
	# neither diagnosed nor floored and needs no special case.
	expect_silent(out <- fetwfe:::.floor_variance_diag(c(NA, -1e-15, 2), "t"))
	expect_equal(out, c(NA, 0, 2))
	# A warning-range negative alongside an NA still warns, and the NA does
	# not enter the index list or the "most negative" figure.
	ws <- .clf_collect(
		fetwfe:::.floor_variance_diag(c(NA, -1e-8), "t"),
		"fetwfe_negative_variance_floored"
	)
	expect_length(ws, 1L)
	expect_true(
		grepl("(indices 2;", conditionMessage(ws[[1]]), fixed = TRUE)
	)
	# Non-numeric passes through unchanged.
	expect_silent(out_chr <- fetwfe:::.floor_variance_diag("nope", "t"))
	expect_equal(out_chr, "nope")
	expect_silent(out_e <- fetwfe:::.floor_variance_diag(numeric(0), "t"))
	expect_equal(out_e, numeric(0))
})

test_that(".floor_cluster_quad_diag floors the diagonal and nothing else", {
	M <- matrix(c(-1e-8, -7, 3, 2), 2, 2)
	ws <- .clf_collect(
		out <- fetwfe:::.floor_cluster_quad_diag(
			M,
			"assemble_joint_cov_var1/Sigma_1"
		),
		"fetwfe_negative_variance_floored"
	)
	expect_length(ws, 1L)
	# OFF-DIAGONALS SURVIVE UNTOUCHED, negative ones included: they can
	# legitimately take either sign.
	expect_equal(out[1, 2], 3)
	expect_equal(out[2, 1], -7)
	expect_equal(diag(out), c(0, 2))
	# The subject IS the scalar helper's exact phrase here, on purpose -- the
	# integration smoke test at the bottom of this file filters on it, so its
	# coverage widens to this site for free.
	expect_true(grepl(
		"Negative cluster-sandwich quadratic form",
		conditionMessage(ws[[1]]),
		fixed = TRUE
	))
	expect_true(grepl(
		"site 'assemble_joint_cov_var1/Sigma_1'",
		conditionMessage(ws[[1]]),
		fixed = TRUE
	))
})

test_that(".floor_cluster_quad_diag passes through what it should", {
	# Positives and exact zeros are untouched, silently.
	P <- matrix(c(1, 0.5, 0.5, 0), 2, 2)
	expect_silent(outP <- fetwfe:::.floor_cluster_quad_diag(P, "t"))
	expect_equal(outP, P)
	# FP-noise negatives clip silently.
	expect_silent(
		outN <- fetwfe:::.floor_cluster_quad_diag(
			matrix(c(-1e-15, 1, 1, -1e-11), 2, 2),
			"t"
		)
	)
	expect_equal(diag(outN), c(0, 0))
	# Catastrophic entries error, classed.
	e <- tryCatch(
		fetwfe:::.floor_cluster_quad_diag(
			matrix(c(-2.14, 0, 0, 1), 2, 2),
			"assemble_joint_cov_var1/Sigma_1"
		),
		error = function(e) e
	)
	expect_s3_class(e, "fetwfe_negative_variance_catastrophic")
	# Sentence-cased at the error tier, so a lowercase `fixed = TRUE` literal
	# does NOT match it -- which is what keeps the smoke-test filter above
	# matching the warning tier only.
	expect_true(grepl(
		"Cluster-sandwich quadratic form on the covariance diagonal",
		conditionMessage(e),
		fixed = TRUE
	))
	expect_false(grepl(
		"Negative cluster-sandwich quadratic form",
		conditionMessage(e),
		fixed = TRUE
	))
	# A character matrix passes through byte-identically (the numeric test
	# lives in the core, so `diag()` still reads and rewrites the same
	# values).
	chr <- matrix(letters[1:4], 2, 2)
	expect_silent(out_chr <- fetwfe:::.floor_cluster_quad_diag(chr, "t"))
	expect_identical(out_chr, chr)
})

test_that(".floor_cluster_quad_diag REJECTS a non-matrix argument", {
	# A non-two-dimensional argument is a PROGRAMMING ERROR, not an input to
	# tolerate: `diag(5)` BUILDS a 5 x 5 identity rather than reading a
	# diagonal, and the only caller always hands this function a matrix. The
	# silent pass-through this `stop()` replaced is the source half of the
	# defect A10 below pins -- a one-token swap at the SCALAR
	# `.compute_att_var1()` site reverted both the #139 diagnostic and the
	# #84-item-9 `max(q, 0)` floor with the whole suite green (#476).
	expect_error(
		fetwfe:::.floor_cluster_quad_diag(5, "t"),
		"requires a two-dimensional matrix",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.floor_cluster_quad_diag(-3.7, "probe_att_var1"),
		"site 'probe_att_var1'",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.floor_cluster_quad_diag(c(1, 2, 3), "t"),
		"numeric of length 3",
		fixed = TRUE
	)
	# Deliberately UNCLASSED: this is a bug report about the caller, not a
	# variance diagnostic, so `.fit_band_for_family()` must neither capture it
	# nor re-raise it below its `tryCatch()` the way it does the two classed
	# tiers -- it degrades to NULL there like any other ordinary error.
	e <- tryCatch(
		fetwfe:::.floor_cluster_quad_diag(5, "t"),
		error = function(e) e
	)
	expect_identical(class(e), c("simpleError", "error", "condition"))
})

test_that("both tiers of the vectorized family carry their condition class", {
	# The classes are the interface `.fit_band_for_family()` keys on, and the
	# two wrappers SHARE them on purpose so one handler pair covers the family.
	w_mat <- .clf_collect(
		fetwfe:::.floor_cluster_quad_diag(matrix(c(-1e-8, 0, 0, 1), 2, 2), "t"),
		"fetwfe_negative_variance_floored"
	)[[1]]
	w_vec <- .clf_collect(
		fetwfe:::.floor_variance_diag(-1e-8, "t"),
		"fetwfe_negative_variance_floored"
	)[[1]]
	for (w in list(w_mat, w_vec)) {
		expect_identical(
			class(w),
			c("fetwfe_negative_variance_floored", "warning", "condition")
		)
	}
	e_mat <- tryCatch(
		fetwfe:::.floor_cluster_quad_diag(matrix(c(-2, 0, 0, 1), 2, 2), "t"),
		error = function(e) e
	)
	e_vec <- tryCatch(
		fetwfe:::.floor_variance_diag(-2, "t"),
		error = function(e) e
	)
	for (e in list(e_mat, e_vec)) {
		expect_identical(
			class(e),
			c(
				"fetwfe_negative_variance_catastrophic",
				"error",
				"condition"
			)
		)
	}
})

test_that(".floor_psd_diag_core caps the rendered index list", {
	# A large `K` must not produce a multi-kilobyte condition message; #431
	# shipped a half-megabyte one from a caller that captured the condition.
	ws <- .clf_collect(
		fetwfe:::.floor_psd_diag_core(rep(-1e-8, 25), "t", "variance"),
		"fetwfe_negative_variance_floored"
	)
	msg <- conditionMessage(ws[[1]])
	expect_true(grepl("25 of 25 entries", msg, fixed = TRUE))
	expect_true(grepl(
		"indices 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...",
		msg,
		fixed = TRUE
	))
	expect_false(grepl("11", msg, fixed = TRUE))
	expect_lt(nchar(msg), 400L)
})

test_that(".floor_psd_diag_core respects custom thresholds", {
	# The core is the only member of the family that names a threshold, and
	# A5b below pins that. A tighter warn threshold makes a silent value warn.
	ws <- .clf_collect(
		fetwfe:::.floor_psd_diag_core(
			-1e-12,
			"t",
			"variance",
			warn_threshold = -1e-13
		),
		"fetwfe_negative_variance_floored"
	)
	expect_length(ws, 1L)
	# A looser error threshold makes a warning-range value error.
	expect_s3_class(
		tryCatch(
			fetwfe:::.floor_psd_diag_core(
				-0.5,
				"t",
				"variance",
				err_threshold = -0.1
			),
			error = function(e) e
		),
		"fetwfe_negative_variance_catastrophic"
	)
})

# --- Coverage / regression test ----------------------------------------------
# Assert every cluster-sandwich quadratic-form site in the package routes
# through one of the recognized floor functions, catching a future revert to
# the bare `max(., 0)` / `pmax(diag(.), 0)` floor. See the file header for the
# predicate and the known false-positive classes.
#
# The package-agnostic AST primitives (`.ns_*`) come from
# helper-namespace-inspect.R, which testthat sources before this file under
# both `devtools::test()` and `test_check()` inside `R CMD check`. Everything
# below is specific to the cluster floor and deliberately stays here.

# The matrix operators a quadratic form can be spelled with.
.clf_quad_ops <- c("%*%", "crossprod", "tcrossprod")

# The floor functions this guardrail RECOGNIZES: a cluster-sandwich quadratic
# form counts as floored when it sits lexically inside a call to one of these.
# `.floor_cluster_quad()` is the scalar helper (#139); `.floor_cluster_quad_diag()`
# is its matrix-valued sibling, which floors a covariance diagonal (#470).
#
# `.floor_variance_diag()` is DELIBERATELY ABSENT -- see WHAT GETS THROUGH in
# the header. Its arguments are model-based diagonals, not sandwich forms.
#
# Membership here says a call counts as a floor; it does NOT say the two are
# substitutable at a given site. A10 pins the owner-to-callee pairing, because
# they are not.
.clf_floor_fns <- c(".floor_cluster_quad", ".floor_cluster_quad_diag")

# Calls that swallow the #139 diagnostic if one wraps a floor call. Wrapping
# the call in `suppressWarnings()` returns the site to pre-#139 behavior while
# leaving every other assertion satisfied, and "this warning is noisy, wrap it"
# is a one-line edit somebody makes without thinking about #139 at all.
.clf_suppressors <- c(
	"suppressWarnings",
	"suppressMessages",
	"try",
	"tryCatch",
	"withCallingHandlers"
)

# Calls that may sit between a `.floor_cluster_quad()` call and the quadratic
# form it floors. Anything else changes the value the floor is asked to judge.
# `abs()` is the case that motivated this: `.floor_cluster_quad(abs(form), ...)`
# leaves the call, the label, the arity and both inventories intact while making
# the diagnostic unable to fire at all -- and `max(abs(q), 0)` is `abs(q)`, so a
# genuinely negative quadratic form comes back POSITIVE rather than floored to
# zero. That inflates a standard error instead of zeroing it, which is a wrong
# answer rather than a missing warning, and it is the most natural way in:
# wrapping in `abs()` is what somebody reaches for to silence a
# negative-variance complaint.
.clf_transparent <- c("as.numeric", "drop", "(", "as.vector", "unname")

# Is any ancestor node a call to one of `fn`?
.clf_any_ancestor <- function(ancestors, fn) {
	for (a in ancestors) {
		if (.ns_is_call_to(a, fn)) {
			return(TRUE)
		}
	}
	FALSE
}

# The name of the function a call invokes, with any `::` / `:::` qualifier
# stripped, or `NA_character_` when the head is not a symbol.
#
# The unwrapping MIRRORS `.ns_is_call_to()`, which handles a namespace-qualified
# head deliberately. Two callers depend on it. `.clf_floor_first_formal()` below
# would return NULL on `fetwfe:::.floor_cluster_quad_diag(...)` -- loud rather
# than blinding (A8 and A9 would go red), but avoidable. A10 compares CALLEE
# NAMES against a literal set, so there the qualified and bare spellings must
# resolve to the SAME string or a correct site reads as the wrong one.
.clf_callee_name <- function(cl) {
	fn_head <- cl[[1]]
	if (
		is.call(fn_head) &&
			length(fn_head) == 3L &&
			is.symbol(fn_head[[1]]) &&
			as.character(fn_head[[1]]) %in% c("::", ":::")
	) {
		fn_head <- fn_head[[3]]
	}
	if (!is.symbol(fn_head)) {
		return(NA_character_)
	}
	as.character(fn_head)
}

# The name of a called floor function's FIRST formal, resolved through the
# package namespace, or NULL when the head cannot be resolved to a function.
#
# Derived rather than hardcoded because the family's members do not agree on
# it: `.floor_cluster_quad()`'s is `q` and `.floor_cluster_quad_diag()`'s is
# `M`. A hardcoded `"q"` made `.floor_cluster_quad_diag(M = form, site = "x")`
# -- behaviour-identical to the positional spelling -- read as having no
# argument at all, turning A8 red on correct code (#470).
.clf_floor_first_formal <- function(floor_call) {
	nm <- .clf_callee_name(floor_call)
	if (is.na(nm)) {
		return(NULL)
	}
	ns <- asNamespace("fetwfe")
	if (!exists(nm, envir = ns, inherits = FALSE)) {
		return(NULL)
	}
	f <- get(nm, envir = ns, inherits = FALSE)
	if (!is.function(f)) {
		return(NULL)
	}
	fmls <- names(formals(f))
	if (length(fmls) == 0L) {
		return(NULL)
	}
	fmls[[1]]
}

# The expression a floor call is actually asked to floor, matched the way R
# matches it: the first formal named EXPLICITLY wins wherever it sits in source
# order, otherwise the first positional argument. Reading `call[[2]]` instead
# would reject `.floor_cluster_quad(site = "x", q = form)`, which is
# behaviour-identical.
.clf_floor_q_arg <- function(floor_call) {
	args <- as.list(floor_call)[-1]
	if (length(args) == 0L) {
		return(NULL)
	}
	nms <- names(args)
	if (is.null(nms)) {
		nms <- rep("", length(args))
	}
	first_formal <- .clf_floor_first_formal(floor_call)
	if (!is.null(first_formal) && first_formal %in% nms) {
		return(args[[which(nms == first_formal)[1]]])
	}
	positional <- which(nms == "")
	if (length(positional) == 0L) {
		return(NULL)
	}
	args[[positional[1]]]
}

# Index of the INNERMOST recognized floor call among `ancestors`
# (outermost-first), or NA if the form is not floored at all.
.clf_floor_ancestor_index <- function(ancestors) {
	hit <- NA_integer_
	for (i in seq_along(ancestors)) {
		if (.ns_is_call_to(ancestors[[i]], .clf_floor_fns)) {
			hit <- i
		}
	}
	hit
}

# Does the quadratic form reach its floor call UNALTERED? Two conditions, both
# needed: every node between the floor call and the form is a transparent
# wrapper, and the chain hangs off the floor call's FIRST argument (element 2)
# rather than off the label or a threshold. Returns TRUE for an unfloored form
# so that this predicate reports only the sanitizer defect -- A1 already owns
# "not floored at all".
.clf_reaches_floor_cleanly <- function(ancestors, node, floor_idx) {
	if (is.na(floor_idx)) {
		return(TRUE)
	}
	floor_call <- ancestors[[floor_idx]]
	between <- if (floor_idx < length(ancestors)) {
		ancestors[(floor_idx + 1L):length(ancestors)]
	} else {
		list()
	}
	for (a in between) {
		if (!.ns_is_call_to(a, .clf_transparent)) {
			return(FALSE)
		}
	}
	chain_head <- if (length(between) > 0L) between[[1]] else node
	q_arg <- .clf_floor_q_arg(floor_call)
	!is.null(q_arg) && identical(q_arg, chain_head)
}

# The symbols inside one function body that hold the cluster-robust sandwich:
# `sandwich_full` plus anything assigned directly from it, to a fixed point so
# a chain `a <- sandwich_full; b <- a` is covered. Without this, a revert that
# assigns the matrix to a local first is invisible -- no site is detected, so
# the universal "every detected site is floored" is satisfied by finding
# nothing.
#
# The set is MONOTONE: a symbol later reassigned to something else is never
# removed, and an assignment inside a branch that never runs still counts. So
# it can over-report a site and can never under-report one. That direction is
# deliberate -- a false positive puts a human in the loop, a false negative is
# a silently unfloored standard error. It is not dataflow analysis and does
# not try to be.
.clf_alias_symbols <- function(code_exprs) {
	aliases <- "sandwich_full"
	repeat {
		found <- aliases
		for (code in code_exprs) {
			.ns_walk_ast(code, function(node, ancestors) {
				if (
					.ns_is_call_to(node, c("<-", "=", "<<-")) &&
						length(node) == 3L &&
						is.symbol(node[[2]]) &&
						is.symbol(node[[3]]) &&
						as.character(node[[3]]) %in% found
				) {
					found <<- union(found, as.character(node[[2]]))
				}
			})
		}
		if (setequal(found, aliases)) {
			break
		}
		aliases <- found
	}
	aliases
}

# One pass over a namespace's functions, returning
#   $sites        one entry per OUTERMOST cluster-sandwich quadratic form: the
#                 owning function, whether a recognized floor call
#                 (`.clf_floor_fns`) is among its AST ancestors, and its
#                 deparsed text.
#   $floor_calls  one entry per recognized floor call: the owning function, the
#                 floor function it CALLS, its argument count, and whether a
#                 condition-suppressing call is among ITS ancestors.
#
# "Outermost" drops the inner node of a chain: `t(a) %*% B %*% a` parses as
# `(t(a) %*% B) %*% a`, so without it one form would be reported twice.
.clf_scan <- function(fns) {
	sites <- list()
	floor_calls <- list()
	for (nm in names(fns)) {
		code_exprs <- .ns_code_exprs(fns[[nm]])
		if (length(code_exprs) == 0L) {
			next
		}
		aliases <- .clf_alias_symbols(code_exprs)
		quad_nodes <- list()
		call_nodes <- list()
		for (code in code_exprs) {
			.ns_walk_ast(code, function(node, ancestors) {
				if (
					.ns_is_call_to(node, .clf_quad_ops) &&
						.ns_subtree_has_symbol(node, aliases)
				) {
					quad_nodes[[length(quad_nodes) + 1L]] <<- list(
						node = node,
						ancestors = ancestors
					)
				}
				if (.ns_is_call_to(node, .clf_floor_fns)) {
					call_nodes[[length(call_nodes) + 1L]] <<- list(
						node = node,
						ancestors = ancestors
					)
				}
			})
		}
		for (quad in quad_nodes) {
			n_anc <- length(quad$ancestors)
			parent <- if (n_anc > 0L) quad$ancestors[[n_anc]] else NULL
			if (
				!is.null(parent) &&
					.ns_is_call_to(parent, .clf_quad_ops) &&
					.ns_subtree_has_symbol(parent, aliases)
			) {
				next
			}
			floor_idx <- .clf_floor_ancestor_index(quad$ancestors)
			sites[[length(sites) + 1L]] <- list(
				fn = nm,
				floored = !is.na(floor_idx),
				judged = .clf_reaches_floor_cleanly(
					quad$ancestors,
					quad$node,
					floor_idx
				),
				text = paste(deparse(quad$node), collapse = " ")
			)
		}
		for (cl in call_nodes) {
			q_arg <- .clf_floor_q_arg(cl$node)
			carries_form <- FALSE
			if (!is.null(q_arg)) {
				.ns_walk_ast(q_arg, function(node, ancestors) {
					if (
						.ns_is_call_to(node, .clf_quad_ops) &&
							.ns_subtree_has_symbol(node, aliases)
					) {
						carries_form <<- TRUE
					}
				})
			}
			floor_calls[[length(floor_calls) + 1L]] <- list(
				fn = nm,
				callee = .clf_callee_name(cl$node),
				nargs = length(cl$node) - 1L,
				carries_form = carries_form,
				suppressed = .clf_any_ancestor(
					cl$ancestors,
					.clf_suppressors
				)
			)
		}
	}
	list(sites = sites, floor_calls = floor_calls)
}

test_that("every cluster-sandwich floor routes through .floor_cluster_quad", {
	bodies <- .ns_deparsed_code("fetwfe")
	fns <- .ns_functions("fetwfe")
	scanned <- .clf_scan(fns)

	# --- A1: the structural universal -------------------------------------
	# Every collected cluster-sandwich quadratic form is floored in place.
	# EXACT set equality against an empty vector, not containment.
	#
	# A1's non-vacuity rests on A4, not on this set: A4 pins the floor-call
	# inventory as a LITERAL, so a `.clf_scan()` that silently collected
	# nothing leaves `floor_call_fns` empty and A4 fires. (Measured: emptying
	# the walk fires A4 alone -- A1 and A6 stay green.) A6 is an independent
	# text inventory over the deparsed bodies and covers a DIFFERENT evasion,
	# a renamed or aliased `sandwich_full`; it does not observe the walk at
	# all.
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	expected_unfloored <- character(0)
	unfloored <- sort(unique(vapply(
		Filter(function(s) !s$floored, scanned$sites),
		`[[`,
		character(1),
		"fn"
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(unfloored, expected_unfloored)

	# --- A8: the floored form reaches its floor UNALTERED ------------------
	# A1 asks whether a form is floored; it cannot ask whether the floor is
	# able to judge what it was handed. `.floor_cluster_quad(abs(form), ...)`
	# satisfies A1, A4 and every inventory while guaranteeing the diagnostic
	# never fires -- and it returns `abs(q)`, so a negative quadratic form is
	# inflated rather than floored. Only `as.numeric()`, `drop()` and parens
	# may sit between the floor call and the form, and the chain must hang
	# off the floor call's first argument.
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	sanitized <- sort(unique(vapply(
		Filter(function(site) !site$judged, scanned$sites),
		`[[`,
		character(1),
		"fn"
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(sanitized, character(0))

	# --- A9: every floor call is handed a real form (no decoys) ------------
	# A8 asks whether a form reaches its floor unaltered. A9 asks the converse,
	# and it is the converse that closes a COMPLETE revert: A4 and A2 are
	# satisfied by a `.floor_cluster_quad()` call on ANY argument, a constant
	# or an already-floored value included. So moving the arithmetic into a
	# helper with a generic formal -- which A1 and A6 cannot see -- while
	# leaving a decoy floor call behind keeps every other assertion green
	# while the #139 diagnostic can never fire again. Measured: that mutation
	# passed every other assertion and the whole suite before this one
	# existed. What A9 verifies is that the floor call is handed a form at
	# all -- not that the form is the right one, nor that the floored value is
	# used. A decoy on a degenerate input, or one whose result is discarded,
	# is in WHAT GETS THROUGH in the header.
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	decoys <- sort(unique(vapply(
		Filter(function(fc) !fc$carries_form, scanned$floor_calls),
		`[[`,
		character(1),
		"fn"
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(decoys, character(0))

	# --- A2: per-site label coverage, exactly once ------------------------
	# Anchored with the surrounding double quotes: `deparse()` renders string
	# literals including their quotes, so an unanchored
	# "cohort_time_atts/var_1" would also match a future
	# "cohort_time_atts/var_12". Counted over OCCURRENCES with `gregexpr()`
	# rather than over elements containing a match -- `sum(grepl(...))` scores
	# two calls sharing one deparsed line as one and passes, which is exactly
	# the copy-paste case this assertion advertises catching. `fixed = TRUE`
	# throughout, stated rather than relied upon.
	expected_labels <- c(
		"getTeResultsOLS/att_var_1",
		"getTeResults2/att_var_1",
		"getCohortATTsFinal/cohort_te_se",
		"event_study_etwfe_betwfe/var_1_e",
		"event_study_fetwfe/var_1_e",
		"cohort_time_atts/var_1",
		"assemble_joint_cov_var1/Sigma_1"
	)
	all_body_text <- paste(unlist(bodies), collapse = " ")
	label_counts <- vapply(
		expected_labels,
		function(lbl) {
			hits <- gregexpr(
				paste0('"', lbl, '"'),
				all_body_text,
				fixed = TRUE
			)[[1]]
			if (hits[[1]] == -1L) 0L else length(hits)
		},
		integer(1)
	)
	# `sprintf()` rather than `paste0()`: `paste0("x", character(0))` recycles
	# the zero-length argument to `""` and returns `"x"`, so an empty offender
	# list would render as one bogus entry and fail on a clean tree.
	# `sprintf()` returns `character(0)` for zero-length inputs.
	off <- label_counts != 1L
	wrong_counts <- sprintf(
		"%s (%d occurrences)",
		expected_labels[off],
		label_counts[off]
	)
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(wrong_counts, character(0))

	# --- A3: the negative-direction guard ---------------------------------
	# No bare `max(... sandwich_full ..., 0)` may remain. Scanned PER BODY,
	# not over a concatenation of all of them -- that is a change of meaning,
	# not a port of the old file-text scan. The regex has a bounded reach in
	# characters, so joining every body makes adjacency an artifact of
	# alphabetical name order and a `max(` in one function can pair with a
	# `sandwich_full` in the next, which both fakes detections and can go red
	# on a correct tree. Per-body scanning also lets the failure name the
	# offending function.
	bare_max_re <- "max\\(.{0,300}sandwich_full.{0,200}?,\\s*0\\)"
	# Scanned per EXPRESSION, not per function: joining a body to its formal
	# defaults with a space recreates, within one function, the same
	# adjacency artifact that joining whole bodies created across functions.
	bare_max_fns <- sort(names(fns)[vapply(
		fns,
		function(f) {
			any(vapply(
				.ns_code_exprs(f),
				function(e) {
					txt <- paste(
						deparse(e, control = c("keepInteger", "keepNA")),
						collapse = " "
					)
					grepl(bare_max_re, txt, perl = TRUE)
				},
				logical(1)
			))
		},
		logical(1)
	)])
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(bare_max_fns, character(0))

	# --- A4: the floor-call inventory -------------------------------------
	# The set of namespace functions that call `.floor_cluster_quad()`, as an
	# EXACT set against a LITERAL expected vector.
	#
	# The expectation is written out here on purpose. An expectation
	# re-derived from the namespace is `setequal(x, x)` and passes every
	# mutation, including a bare revert of the floor at any site. What went
	# stale in #463 was the hand-maintained FILE LIST that produced the
	# OBSERVATION -- and the observation is what is namespace-derived below.
	#
	# This is the absolute pin: it does not care how the arithmetic is
	# spelled, so it survives the aliasing and `crossprod` refactors that
	# defeat a shape-based predicate, and it is what makes A1 non-vacuous --
	# without it, a site that stops being DETECTED satisfies A1 by absence.
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	expected_floor_call_fns <- c(
		".assemble_joint_cov_var1",
		".compute_att_var1",
		".event_study_etwfe_betwfe",
		".event_study_fetwfe",
		"cohortTimeATTs",
		"getCohortATTsFinal"
	)
	floor_call_fns <- sort(unique(vapply(
		scanned$floor_calls,
		`[[`,
		character(1),
		"fn"
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(floor_call_fns, expected_floor_call_fns)

	# --- A10: WHICH floor function each site calls ------------------------
	# A4 pins the set of functions that contain a floor call. It does NOT pin
	# which member of `.clf_floor_fns` each one calls, and the two are not
	# interchangeable: `.floor_cluster_quad()` takes a SCALAR quadratic form,
	# `.floor_cluster_quad_diag()` a MATRIX whose diagonal it floors. Handing
	# either the other's argument is a silent revert of the floor, not a type
	# error, so A4 and every other assertion here stay green through it.
	#
	# Measured (#476): swapping the one token `.floor_cluster_quad(` to
	# `.floor_cluster_quad_diag(` at `.compute_att_var1()` -- arguments
	# untouched, since both sites write the form literally inside the call --
	# reverted the #139 diagnostic AND the underlying #84-item-9 `max(q, 0)`
	# floor, returned a NEGATIVE `att_var_1` and an understated overall-ATT
	# standard error, and left the whole 5440-assertion suite BYTE-IDENTICALLY
	# green. The dimensionless scalar hit `.floor_cluster_quad_diag()`'s old
	# `if (length(dim(M)) != 2L) return(M)` pass-through and came straight back.
	#
	# Note the direction the guardrail had it backwards. The SAFE edit --
	# swapping to `.floor_variance_diag()`, which really does floor a scalar --
	# was REJECTED (A1 reports the form as unfloored and A4 loses
	# `.compute_att_var1`, because that helper is deliberately outside
	# `.clf_floor_fns`), while the DESTRUCTIVE one was accepted. A10 and
	# `.floor_cluster_quad_diag()`'s own `stop()` are the two halves of the fix;
	# on the destructive swap A10 now fires, and so do the fifteen-odd runtime
	# assertions that reach a cluster-SE fit.
	#
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	expected_floor_callees <- c(
		".assemble_joint_cov_var1 -> .floor_cluster_quad_diag",
		".compute_att_var1 -> .floor_cluster_quad",
		".event_study_etwfe_betwfe -> .floor_cluster_quad",
		".event_study_fetwfe -> .floor_cluster_quad",
		"cohortTimeATTs -> .floor_cluster_quad",
		"getCohortATTsFinal -> .floor_cluster_quad"
	)
	floor_callees <- sort(unique(vapply(
		scanned$floor_calls,
		function(fc) paste0(fc$fn, " -> ", fc$callee),
		character(1)
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(floor_callees, expected_floor_callees)

	# --- A6: the sandwich-mention inventory -------------------------------
	# The set of namespace functions whose deparsed body mentions
	# `sandwich_full`, as an exact set against a literal vector. Its
	# predicate is exactly that -- the literal string, anywhere in the body --
	# which is what lets it survive the aliasing and `crossprod` refactors
	# that defeat a shape-based check: a NEW function that reaches the matrix
	# through a list element, a submatrix slice or `assign()`, or that
	# receives it as a formal still SPELLED `sandwich_full`, lands in this set
	# even though A1 cannot see the arithmetic. It also closes the coupling in
	# the SHRINK direction: an already-listed function that stops spelling the
	# symbol drops out and fails.
	#
	# What it does NOT do is the grow direction inside existing code -- see
	# KNOWN BLIND SPOTS in the file header. (A `::`-qualified operator is
	# caught by A1, not here; see `.ns_is_call_to()`.)
	expected_sandwich_fns <- c(
		".assemble_cluster_robust_sandwich",
		".assemble_joint_cov_var1",
		".compute_att_var1",
		".event_study_etwfe_betwfe",
		".event_study_fetwfe",
		".ols_estimator_core",
		".recompute_gram_and_sandwich",
		".simultaneous_cis_impl",
		"betwfe_core",
		"cohortTimeATTs",
		"fetwfe_core",
		"getCohortATTsFinal",
		"getTeResults2",
		"getTeResultsOLS"
	)
	sandwich_fns <- sort(names(bodies)[vapply(
		bodies,
		function(x) any(grepl("sandwich_full", x, fixed = TRUE)),
		logical(1)
	)])
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(sandwich_fns, expected_sandwich_fns)

	# --- A7: no condition suppression around a floor call -----------------
	suppressed_fns <- sort(unique(vapply(
		Filter(function(fc) fc$suppressed, scanned$floor_calls),
		`[[`,
		character(1),
		"fn"
	)))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(suppressed_fns, character(0))

	# --- A7b: ... and none at the CALLER, one frame up --------------------
	# A7 above scans the ancestors of the floor call within its own body, so
	# `suppressWarnings(.compute_att_var1(...))` at the call site passes it
	# while suppressing exactly the same diagnostic. That is the MORE likely
	# edit of the two -- wrapping a whole internal call is easier to write
	# than wrapping one inner expression. Measured: without this assertion
	# that mutation leaves the suite fully green.
	wrapped_callers <- character(0)
	for (nm in names(fns)) {
		for (code in .ns_code_exprs(fns[[nm]])) {
			.ns_walk_ast(code, function(node, ancestors) {
				if (
					.ns_is_call_to(node, .clf_suppressors) &&
						.ns_subtree_has_symbol(node, expected_floor_call_fns)
				) {
					wrapped_callers <<- union(wrapped_callers, nm)
				}
			})
		}
	}
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(sort(wrapped_callers), character(0))
})

test_that("no call site neuters .floor_cluster_quad's diagnostic contract", {
	bodies <- .ns_deparsed_code("fetwfe")
	scanned <- .clf_scan(.ns_functions("fetwfe"))

	# --- A5a: arity -- every floor call passes exactly two arguments ------
	# This is the load-bearing one: a per-call-site threshold override spelled
	# POSITIONALLY -- `.floor_cluster_quad(form, "label", -Inf, -Inf)` --
	# mentions no parameter name, so the scan below cannot see it, and it
	# returns the site to a silent `max(q, 0)` with the whole suite green.
	# The arity rule rejects both spellings and does not depend on the
	# parameters keeping their current names.
	bad_arity <- vapply(
		Filter(function(fc) fc$nargs != 2L, scanned$floor_calls),
		function(fc) paste0(fc$fn, " (", fc$nargs, " args)"),
		character(1)
	)
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(bad_arity, character(0))

	# --- A5b: no per-call-site threshold override, spelled by NAME --------
	# No namespace function other than the two helpers that legitimately CARRY
	# the thresholds may mention either: a NAMED per-call-site override is the
	# other spelling of the same edit A5a catches positionally. Exactly those
	# two and no wider -- `.floor_cluster_quad_diag()` and
	# `.floor_variance_diag()` must NOT name a threshold, which is what keeps
	# A5b's real target closed. Being an exact set, it also PARTLY bounds the
	# deliberate duplication between the scalar helper and
	# `.floor_psd_diag_core()`: a third copy of the tiering logic goes red here
	# automatically IF it names either threshold. One that writes the constants
	# inline (`if (any(v < -1))`) names neither and is invisible to the whole
	# guardrail -- measured, all green -- since it has no call site for A1/A4/A9
	# to collect and is not in A5c's pin set. Do not read this assertion as
	# bounding the duplication outright; it bounds one of the two spellings,
	# and the inline one is the likelier accident.
	#
	# This is also the only assertion in the file that covers
	# `.floor_variance_diag()`'s call sites at all -- A5a iterates the floor
	# calls `.clf_scan()` collected via `.clf_floor_fns`, from which it is
	# deliberately absent. See WHAT GETS THROUGH in the header.
	threshold_fns <- sort(names(bodies)[vapply(
		bodies,
		function(x) {
			any(grepl("err_threshold", x, fixed = TRUE)) ||
				any(grepl("warn_threshold", x, fixed = TRUE))
		},
		logical(1)
	)])
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(
		threshold_fns,
		c(".floor_cluster_quad", ".floor_psd_diag_core")
	)

	# --- A5c: the helpers' own formals, pinned by name --------------------
	# Pinned by NAME, not by value. The unit tests at the top of this file
	# already catch every way of neutering the DEFAULTS, so a value pin would
	# add nothing; the name pin is what catches a new opt-out formal (say
	# `diagnose = TRUE`, passed `FALSE` at one site), which leaves the
	# defaults untouched and so keeps those unit tests green.
	#
	# EVERY function in the floor family is pinned, stated as that predicate
	# rather than as a count. `.floor_psd_diag_core()` carries the thresholds,
	# so it is the one a new opt-out formal would most usefully subvert, and
	# it also carries two pieces of message machinery -- the index cap and the
	# error tier's sentence-casing -- that live INLINE there on purpose. If
	# either is ever factored into a fourth helper, that helper joins this
	# set.
	expect_identical(
		paste(names(formals(fetwfe:::.floor_cluster_quad)), collapse = ", "),
		"q, site, err_threshold, warn_threshold"
	)
	expect_identical(
		paste(names(formals(fetwfe:::.floor_psd_diag_core)), collapse = ", "),
		"v, site, subject, err_threshold, warn_threshold"
	)
	expect_identical(
		paste(
			names(formals(fetwfe:::.floor_cluster_quad_diag)),
			collapse = ", "
		),
		"M, site"
	)
	expect_identical(
		paste(names(formals(fetwfe:::.floor_variance_diag)), collapse = ", "),
		"v, site"
	)
})

# --- A11: the reachability claim the SCALAR helper's unclassed conditions rest
# on ---------------------------------------------------------------------------
#
# `R/cluster_floor.R`'s header justifies leaving `.floor_cluster_quad()`'s
# `warning()` / `stop()` UNCLASSED on a reachability argument: every one of its
# call sites is outside `.fit_band_for_family()`'s protected region, so there is
# nothing for a class to key on. That claim is load-bearing. If it stops
# holding, a scalar `stop()` raised inside the region is swallowed by
# `error = function(e) NULL` and the band degrades to the POINTWISE one under a
# `[simultaneous 95% CI]` header -- a narrower interval, so it over-rejects.
# That is #470's own defect, at the sites #139 was written for.
#
# An untested claim about reachability rots; this one pins it. (Same reasoning,
# and the same sentence, as `test-fit-time-singular-gram-degrade-400.R`.)
#
# The protected region is exactly the transitive closure of
# `.simultaneous_cis_impl()`: that is the single call inside the `tryCatch()`.

# Every namespace function `f` calls, by CALL HEAD, `::` / `:::` unwrapped.
# Deliberately NOT "every symbol or string naming a namespace function": that
# over-approximation is a measured false-positive generator here, because
# `.simultaneous_cis_impl()` mentions the CLASS NAME `"fetwfe"`, which is also
# the name of the top-level estimator, and the resulting phantom edge makes
# nearly half the package look reachable. A12 below is what covers the opposite
# risk -- a real edge this lexical rule cannot see.
.clf_call_edges <- function(fns) {
	nms <- names(fns)
	lapply(fns, function(f) {
		out <- character(0)
		for (code in .ns_code_exprs(f)) {
			.ns_walk_ast(code, function(node, ancestors) {
				if (is.call(node)) {
					nm <- .clf_callee_name(node)
					if (!is.na(nm) && nm %in% nms) {
						out <<- c(out, nm)
					}
				}
			})
		}
		unique(out)
	})
}

# Transitive closure of `seed` over `edges`, INCLUDING the seed itself.
.clf_reachable_from <- function(edges, seed) {
	seen <- seed
	frontier <- seed
	while (length(frontier) > 0L) {
		nxt <- setdiff(
			unique(unlist(edges[intersect(frontier, names(edges))])),
			seen
		)
		seen <- c(seen, nxt)
		frontier <- nxt
	}
	sort(seen)
}

test_that("no .floor_cluster_quad() site is reachable from the protected region", {
	fns <- .ns_functions("fetwfe")
	edges <- .clf_call_edges(fns)
	reachable <- .clf_reachable_from(edges, ".simultaneous_cis_impl")

	# The scalar helper's callers, derived from the same walk `.clf_scan()`
	# uses. A10 above is what pins this set as a literal; here it is an
	# OBSERVATION, and the claim is the intersection below.
	scalar_callers <- sort(unique(vapply(
		Filter(
			function(fc) identical(fc$callee, ".floor_cluster_quad"),
			.clf_scan(fns)$floor_calls
		),
		`[[`,
		character(1),
		"fn"
	)))

	# Non-vacuity, in both directions. Without these an empty intersection is
	# equally consistent with a broken walk or an empty caller set.
	expect_gt(length(scalar_callers), 0L)
	expect_gt(length(reachable), 1L)
	# The POSITIVE CONTROL: the walk really does reach floor-calling code from
	# the seed. `.assemble_joint_cov_var1()` calls the MATRIX floor and IS
	# inside the protected region -- which is exactly why the vectorized
	# family's conditions are classed and the scalar family's are not.
	expect_true(".assemble_joint_cov_var1" %in% reachable)
	expect_true(".floor_cluster_quad_diag" %in% reachable)

	# --- A11: the claim itself --------------------------------------------
	# If this goes red, do NOT delete it and do NOT relax it to containment:
	# class the scalar helper's two conditions and teach
	# `.fit_band_for_family()` to capture and re-raise them, exactly as it does
	# the vectorized family's, then update `R/cluster_floor.R`'s header.
	expect_setequal(intersect(reachable, scalar_callers), character(0))

	# --- A12: nothing in the region dispatches out of the lexical graph ----
	# A11's walk follows call heads, so a call assembled at RUNTIME is invisible
	# to it -- and this package really does dispatch that way elsewhere
	# (`.call_te()` builds a call from a character `te_fn_name` and `eval()`s
	# it, which is how `fetwfe_core()` reaches `getTeResults2()` and hence
	# `.compute_att_var1()`). One such construct inside the protected region
	# would make A11's empty intersection unsound. Exactly one reachable
	# function uses one, and its use cannot name a namespace function:
	# `.build_propensity_if()` calls `do.call(rbind, blocks)`. An EXACT literal,
	# so a second one puts a human in the loop.
	dyn_verbs <- c(
		"eval",
		"evalq",
		"do.call",
		"match.fun",
		"get",
		"get0",
		"mget",
		"Recall",
		"getFromNamespace",
		"getExportedValue"
	)
	dynamic_dispatchers <- sort(Filter(
		function(nm) {
			hit <- FALSE
			for (code in .ns_code_exprs(fns[[nm]])) {
				.ns_walk_ast(code, function(node, ancestors) {
					if (.ns_is_call_to(node, dyn_verbs)) {
						hit <<- TRUE
					}
				})
			}
			hit
		},
		reachable
	))
	# The expected set is a LITERAL on purpose. Never regenerate it from a
	# failure's `Needs:` / `Absent:` output -- that makes it `setequal(x, x)`.
	expect_setequal(dynamic_dispatchers, ".build_propensity_if")
})

# --- Integration smoke test --------------------------------------------------
# Fit a small cluster-SE model on well-conditioned simulated data and
# verify no warning/error from `.floor_cluster_quad` fires.

test_that("cluster-SE fit on well-conditioned data does not trigger the diagnostic", {
	# Small but well-conditioned: same recipe as several existing tests.
	set.seed(2026)
	coefs <- genCoefs(
		G = 3,
		T = 5,
		density = 0.5,
		eff_size = 1,
		d = 1,
		seed = 2026
	)
	sim <- simulateData(
		coefs,
		N = 80,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 2026
	)
	# Run the fit; assert no warning at all surfaces with the message
	# pattern from the floor family. (Other unrelated warnings are not the
	# concern of this test; we filter by substring.)
	#
	# The filter phrase is rendered by BOTH the scalar helper (#139) and
	# `.floor_cluster_quad_diag()`'s warning tier (#470), so this test's
	# coverage widened to the matrix site for free. It cannot tell the two
	# apart, though -- which is why the failure text names the SITE from the
	# message rather than a hardcoded helper name, and why no assertion in
	# `test-matrix-floor-conditions-470.R` keys on this phrase alone.
	withCallingHandlers(
		res <- fetwfeWithSimulatedData(sim, q = 0.5, se_type = "cluster"),
		warning = function(w) {
			msg <- conditionMessage(w)
			if (grepl("cluster-sandwich quadratic form", msg, fixed = TRUE)) {
				stop(
					"Unexpected cluster-floor warning on well-conditioned data: ",
					msg
				)
			}
			invokeRestart("muffleWarning")
		}
	)
	expect_true(!is.null(res$att_se))
	expect_true(is.finite(res$att_se))
})
