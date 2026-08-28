library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for `.floor_cluster_quad()` (issue #139, version 1.11.2), and the
# guardrail keeping every cluster-sandwich quadratic-form site routed through
# it (issue #463).
#
# `.floor_cluster_quad()` layers a two-tier diagnostic on top of the
# pre-existing `max(q, 0)` floor at each cluster-sandwich quadratic-form
# site. The forms are PSD in exact arithmetic, so any negative value is
# either FP-noise (silently clipped to 0) or a bug signal (warn / stop).
#
# The guardrail below states a PREDICATE, not a count and not a file list --
# both of those went stale here once already. It walks the AST of every
# function in the package NAMESPACE and asserts that every cluster-sandwich
# quadratic form (a `%*%` / `crossprod()` / `tcrossprod()` node whose operands
# reach the cluster-robust sandwich) sits lexically inside a
# `.floor_cluster_quad()` call. Reading the namespace instead of `R/*.R` as
# text is what makes it run under `R CMD check`, where the tests execute
# against the *installed* package and there is no `R/` directory at all -- the
# whole block used to skip there, which is every CI job and every CRAN
# machine (#463).
#
# ONE DOCUMENTED EXCEPTION, asserted rather than tolerated (it is part of an
# exact set, so this site silently becoming floored fails too):
#   `.assemble_joint_cov_var1()` in R/variance_machinery.R builds a K x K
#   covariance block `t(Psi_full) %*% sandwich_full %*% Psi_full` and floors
#   only its diagonal, with a bare `pmax(diag(.), 0)` and no #139 diagnostic.
#   It cannot route through `.floor_cluster_quad()`, which is scalar-only by
#   contract: any input of length other than one is passed through unchanged.
#   Extending the diagnostic to the matrix-valued floors is tracked as #470.
#   The exemption is function-granular, so a future scalar form added inside
#   that same function would inherit it silently.
#
# On the per-site labels asserted below: `getTeResultsOLS/att_var_1` and
# `getTeResults2/att_var_1` are two LABELS at one SITE. #344 merged those two
# functions' floors into the shared `.compute_att_var1()`, and each caller
# passes its own label through the `label` formal -- which is also why the
# label assertion has no power at that site, and why the floor-call inventory
# (A4) exists.
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
#   * A scalar form inside `.assemble_joint_cov_var1()`, which inherits that
#     function's exemption. The exemption is function-granular.
#   * `abs()` or any sanitizer applied to an OPERAND inside the product --
#     the path from floor to form is checked, the form's operands are not.
#   * Condition suppression more than ONE frame above a floor call. One frame
#     is checked; a call graph would be needed for more.
#
# None of these is a regression. The block this replaces was blind to all of
# them AND to whole files, and it ran on no automated machine. The single
# thing it did better -- reading formal defaults, which it got for free by
# scanning raw text -- is closed: `.ns_code_exprs()` scans defaults alongside
# bodies, so a bare `max(... sandwich_full ..., 0)` in a default is now
# rejected wherever the function-granular exemption above does not apply.
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

# --- Coverage / regression test ----------------------------------------------
# Assert every cluster-sandwich quadratic-form site in the package routes
# through `.floor_cluster_quad`, catching a future revert to the bare
# `max(., 0)` floor. See the file header for the predicate, the one documented
# exception, and the known false-positive classes.
#
# The package-agnostic AST primitives (`.ns_*`) come from
# helper-namespace-inspect.R, which testthat sources before this file under
# both `devtools::test()` and `test_check()` inside `R CMD check`. Everything
# below is specific to the cluster floor and deliberately stays here.

# The matrix operators a quadratic form can be spelled with.
.clf_quad_ops <- c("%*%", "crossprod", "tcrossprod")

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

# The expression `.floor_cluster_quad()` is actually asked to floor, matched the
# way R matches it: a named `q =` wins wherever it sits in source order,
# otherwise the first positional argument. Reading `call[[2]]` instead would
# reject `.floor_cluster_quad(site = "x", q = form)`, which is
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
	if ("q" %in% nms) {
		return(args[[which(nms == "q")[1]]])
	}
	positional <- which(nms == "")
	if (length(positional) == 0L) {
		return(NULL)
	}
	args[[positional[1]]]
}

# Index of the INNERMOST `.floor_cluster_quad()` call among `ancestors`
# (outermost-first), or NA if the form is not floored at all.
.clf_floor_ancestor_index <- function(ancestors) {
	hit <- NA_integer_
	for (i in seq_along(ancestors)) {
		if (.ns_is_call_to(ancestors[[i]], ".floor_cluster_quad")) {
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
#                 owning function, whether a `.floor_cluster_quad()` call is
#                 among its AST ancestors, and its deparsed text.
#   $floor_calls  one entry per `.floor_cluster_quad()` call: the owning
#                 function, its argument count, and whether a
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
				if (.ns_is_call_to(node, ".floor_cluster_quad")) {
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
	# Every collected cluster-sandwich quadratic form is floored in place,
	# except the one documented matrix-valued site. EXACT set equality, not
	# containment: a listed site that has BECOME floored fails too, which is
	# what stops the walk collapsing to nothing and satisfying the universal
	# by finding no sites at all.
	# A failure here may be correct code -- see WHAT GOES RED in the header.
	expected_unfloored <- c(
		# K x K covariance block; floors its diagonal with `pmax(diag(.), 0)`
		# because `.floor_cluster_quad()` is scalar-only. See the header.
		".assemble_joint_cov_var1"
	)
	unfloored <- sort(unique(vapply(
		Filter(function(s) !s$floored, scanned$sites),
		`[[`,
		character(1),
		"fn"
	)))
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
		"cohort_time_atts/var_1"
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
	expect_setequal(floor_call_fns, expected_floor_call_fns)

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
	expect_setequal(sandwich_fns, expected_sandwich_fns)

	# --- A7: no condition suppression around a floor call -----------------
	suppressed_fns <- sort(unique(vapply(
		Filter(function(fc) fc$suppressed, scanned$floor_calls),
		`[[`,
		character(1),
		"fn"
	)))
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
	expect_setequal(sort(wrapped_callers), character(0))
})

test_that("no call site neuters .floor_cluster_quad's diagnostic contract", {
	bodies <- .ns_deparsed_code("fetwfe")
	scanned <- .clf_scan(.ns_functions("fetwfe"))

	# Arity. Every call passes exactly two arguments. This is the
	# load-bearing one: a per-call-site threshold override spelled
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
	expect_setequal(bad_arity, character(0))

	# No namespace function other than the helper itself may mention either
	# threshold: a NAMED per-call-site override is the other spelling of the
	# same edit.
	threshold_fns <- sort(names(bodies)[vapply(
		bodies,
		function(x) {
			any(grepl("err_threshold", x, fixed = TRUE)) ||
				any(grepl("warn_threshold", x, fixed = TRUE))
		},
		logical(1)
	)])
	expect_setequal(threshold_fns, ".floor_cluster_quad")

	# The helper's formals pinned by NAME, not by value. The unit tests above
	# already catch every way of neutering the DEFAULTS, so a value pin would
	# add nothing; the name pin is what catches a new opt-out formal (say
	# `diagnose = TRUE`, passed `FALSE` at one site), which leaves the
	# defaults untouched and so keeps those unit tests green.
	expect_identical(
		paste(names(formals(fetwfe:::.floor_cluster_quad)), collapse = ", "),
		"q, site, err_threshold, warn_threshold"
	)
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
	# pattern from `.floor_cluster_quad()`. (Other unrelated warnings are
	# not the concern of this test; we filter by substring.)
	withCallingHandlers(
		res <- fetwfeWithSimulatedData(sim, q = 0.5, se_type = "cluster"),
		warning = function(w) {
			msg <- conditionMessage(w)
			if (grepl("cluster-sandwich quadratic form", msg, fixed = TRUE)) {
				stop(
					"Unexpected .floor_cluster_quad warning on well-conditioned data: ",
					msg
				)
			}
			invokeRestart("muffleWarning")
		}
	)
	expect_true(!is.null(res$att_se))
	expect_true(is.finite(res$att_se))
})
