library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Direct unit tests for the shared namespace-inspection primitives defined in
# helper-namespace-inspect.R (issue #463).
#
# WHY THIS FILE EXISTS. Those primitives are the machinery every namespace-
# walking guardrail in this package stands on, and until this file they were
# exercised only *through* one such guardrail, in test-cluster_floor.R. That is
# the wrong direction of coverage: a bug that made a primitive UNDER-report --
# a walker that stops recursing, a formals scan that silently skips defaults, a
# `::`-qualified head that stops matching -- makes the guardrail find fewer
# sites, and a guardrail that finds fewer sites PASSES. The failure would look
# exactly like the defect #463 exists to fix: a green test that never ran the
# check it advertises.
#
# So every assertion below is written in the under-report direction. Node counts
# are pinned rather than bounded, "must find" is preferred to "must not find",
# and each empty-symbol case asserts that the walk continued PAST the empty
# symbol rather than merely that it did not error.
#
# The fixtures are defined here, with known contents, rather than taken from the
# `fetwfe` namespace, so no expectation below can drift when the package
# changes. The two package-level primitives (`.ns_functions()` and
# `.ns_deparsed_code()`) take a package name and cannot be pointed at a fixture;
# they are asserted on drift-free PROPERTIES of the real namespace --
# sortedness, function-ness, body-and-defaults coverage -- never on a list of
# names or a count.
#
# `str2lang()` rather than `quote()` wherever the expression is a function
# definition or contains an empty symbol: it parses with `keep.source = FALSE`,
# so the shape does not depend on the option in force when this file was parsed.
# ------------------------------------------------------------------------------

# Every node `.ns_walk_ast()` visited, in visit order, with its ancestors.
.nsi_visited <- function(expr) {
	seen <- list()
	.ns_walk_ast(expr, function(node, ancestors) {
		seen[[length(seen) + 1L]] <<- list(node = node, ancestors = ancestors)
	})
	seen
}

# Was `target` among the visited nodes?
.nsi_saw <- function(visited, target) {
	any(vapply(visited, function(v) identical(v$node, target), logical(1)))
}

# --- .ns_walk_ast -------------------------------------------------------------

test_that(".ns_walk_ast visits every node of a known expression", {
	# `a + b * c` parses as `+`(a, `*`(b, c)): the outer call, its head `+`,
	# `a`, the inner call, its head `*`, `b`, `c`. Seven nodes, counted rather
	# than bounded -- a walker that stopped recursing into nested calls would
	# report four and still satisfy any "visits at least the root" assertion.
	visited <- .nsi_visited(quote(a + b * c))
	expect_identical(length(visited), 7L)
	for (node in list(
		quote(a + b * c),
		as.symbol("+"),
		quote(a),
		quote(b * c),
		as.symbol("*"),
		quote(b),
		quote(c)
	)) {
		expect_true(.nsi_saw(visited, node))
	}

	# The root itself is visited, with no ancestors, even when it is a leaf.
	leaf <- .nsi_visited(quote(x))
	expect_identical(length(leaf), 1L)
	expect_identical(leaf[[1]]$node, quote(x))
	expect_identical(leaf[[1]]$ancestors, list())
})

test_that(".ns_walk_ast passes ancestors outermost-first", {
	# The whole ancestor-order contract: `.clf_floor_ancestor_index()` in
	# test-cluster_floor.R takes the LAST matching ancestor as the innermost
	# floor call, so a reversed list would silently pick the wrong one.
	visited <- .nsi_visited(quote(a + b * c))
	deepest <- Filter(function(v) identical(v$node, quote(c)), visited)
	expect_identical(length(deepest), 1L)
	anc <- deepest[[1]]$ancestors
	expect_identical(length(anc), 2L)
	expect_identical(anc[[1]], quote(a + b * c))
	expect_identical(anc[[2]], quote(b * c))
	# `ancestors` excludes the node itself.
	expect_false(identical(anc[[length(anc)]], quote(c)))
})

test_that(".ns_walk_ast does not choke on R's empty symbol", {
	# R represents an absent index and a formal with no default with the same
	# EMPTY SYMBOL, and binding one to a local turns the local into a missing
	# argument -- so the walker tests and recurses inline. Each case below
	# asserts the walk CONTINUED past the empty symbol, not just that it
	# survived: a `break` where the helper writes `next` would leave every
	# later sibling unvisited and every guardrail built on this quietly blind.

	# `x[, 1]`: the empty index sits BETWEEN `x` and `1`.
	idx <- .nsi_visited(str2lang("x[, 1]"))
	expect_identical(length(idx), 4L)
	expect_true(.nsi_saw(idx, quote(x)))
	expect_true(.nsi_saw(idx, 1))

	# A formals pairlist: `a` has no default, `b` does, and `b`'s comes second.
	fmls <- .nsi_visited(formals(function(a, b = 2) NULL))
	expect_identical(length(fmls), 2L)
	expect_true(.nsi_saw(fmls, 2))
	# All-empty formals: the pairlist itself, and nothing else.
	expect_identical(length(.nsi_visited(formals(function(a, b) NULL))), 1L)

	# A nested function definition inside a body. A `function` call always has
	# four elements -- head, formals pairlist, body, srcref slot -- and the
	# srcref slot is a single leaf whether it holds NULL (source not kept) or a
	# srcref object, so this count is stable under either `keep.source`.
	# Measured both ways.
	nested <- .nsi_visited(str2lang("function(a, b = 2) a + b"))
	expect_identical(length(nested), 9L)
	expect_true(.nsi_saw(nested, quote(a + b)))
	expect_true(.nsi_saw(nested, 2))

	# An `expression` object, the third container the walker recurses into.
	expect_identical(
		length(.nsi_visited(parse(text = "a; b", keep.source = FALSE))),
		3L
	)
})

# --- .ns_code_exprs -----------------------------------------------------------

test_that(".ns_code_exprs returns the body and every code-carrying default", {
	f <- function(x, y = 1 + 2, z) {
		x + 1
	}
	exprs <- .ns_code_exprs(f)
	# Body first, then defaults in formals order; non-defaulted formals
	# contribute nothing, so the length is exactly two.
	expect_identical(length(exprs), 2L)
	expect_identical(exprs[[1]], body(f))
	expect_identical(exprs[[2]], quote(1 + 2))
})

test_that(".ns_code_exprs sees code that body() alone cannot", {
	# The regression this primitive exists to close, stated as a fixture: a
	# guardrail that scanned raw file text saw formal defaults for free, and a
	# `body()`-only walk does not. `S` appears ONLY in the default here.
	g <- function(x, S, v = max(as.numeric(t(x) %*% S %*% x), 0)) {
		v
	}
	expect_false(.ns_subtree_has_symbol(body(g), "S"))
	expect_true(any(vapply(
		.ns_code_exprs(g),
		function(e) .ns_subtree_has_symbol(e, "S"),
		logical(1)
	)))
})

test_that(".ns_code_exprs handles a NULL body and a primitive", {
	# A NULL body contributes nothing and must not error.
	expect_identical(length(.ns_code_exprs(function(a) NULL)), 0L)

	# A primitive has neither body nor formals.
	expect_identical(length(.ns_code_exprs(sum)), 0L)
	expect_identical(length(.ns_code_exprs(length)), 0L)

	# A default of literal `NULL` is DROPPED, and that is a property of the
	# language rather than a choice: `out[[length(out) + 1L]] <- NULL` deletes
	# rather than appends. Pinned so nobody assumes such a default is scanned.
	# It costs no coverage -- `NULL` is a leaf with no subexpression to find --
	# but a reader who assumed otherwise would mis-read every scan built on
	# this primitive.
	dropped <- .ns_code_exprs(function(a, b = NULL, c = 1 + 2) NULL)
	expect_identical(length(dropped), 1L)
	expect_identical(dropped[[1]], quote(1 + 2))
})

# --- .ns_is_call_to -----------------------------------------------------------

test_that(".ns_is_call_to matches bare and namespace-qualified heads", {
	expect_true(.ns_is_call_to(quote(foo(1)), "foo"))
	# The head of `Matrix::crossprod(x)` is itself a call to `::`, so a naive
	# `is.symbol(expr[[1]])` test declines to match it and a `pkg::op()`
	# spelling escapes every guardrail built on this predicate.
	expect_true(.ns_is_call_to(str2lang("Matrix::crossprod(x)"), "crossprod"))
	expect_true(.ns_is_call_to(
		str2lang("fetwfe:::.floor_cluster_quad(q, 'site')"),
		".floor_cluster_quad"
	))
	# A character VECTOR: the operator and alias sets are both plural.
	expect_true(.ns_is_call_to(quote(a %*% b), c("%*%", "crossprod")))
	expect_true(.ns_is_call_to(quote(crossprod(a, b)), c("%*%", "crossprod")))
})

test_that(".ns_is_call_to declines a non-call and a different function", {
	expect_false(.ns_is_call_to(quote(x), "x"))
	expect_false(.ns_is_call_to(42, "foo"))
	expect_false(.ns_is_call_to("foo", "foo"))
	expect_false(.ns_is_call_to(quote(foo(1)), "bar"))
	# The qualified head resolves to the FUNCTION, so the package name is not
	# what matches.
	expect_false(.ns_is_call_to(str2lang("Matrix::crossprod(x)"), "Matrix"))
})

# --- .ns_subtree_has_symbol ---------------------------------------------------

test_that(".ns_subtree_has_symbol finds a symbol at any depth", {
	expect_true(.ns_subtree_has_symbol(
		str2lang("f(g(h(i(target))))"),
		"target"
	))
	# A character vector in the second argument.
	expect_true(.ns_subtree_has_symbol(quote(a + b), c("z", "b")))
	# A call's head is a symbol node like any other, and is reported as one.
	expect_true(.ns_subtree_has_symbol(quote(foo(1)), "foo"))
	# Reached only by recursing into a formals pairlist...
	expect_true(.ns_subtree_has_symbol(
		str2lang("function(a = target) NULL"),
		"target"
	))
	# ...and only by continuing past an empty symbol.
	expect_true(.ns_subtree_has_symbol(str2lang("x[, target]"), "target"))
})

test_that(".ns_subtree_has_symbol returns FALSE when the symbol is absent", {
	expect_false(.ns_subtree_has_symbol(
		str2lang("f(g(h(i(target))))"),
		"absent_symbol"
	))
	expect_false(.ns_subtree_has_symbol(quote(a + b), character(0)))
	# A string literal is not a symbol.
	expect_false(.ns_subtree_has_symbol(quote(f("target")), "target"))
})

# --- .ns_deparsed_code --------------------------------------------------------

test_that(".ns_deparsed_code's deparse contract drops comments and srcrefs", {
	old <- options(keep.source = TRUE)
	on.exit(options(old), add = TRUE)
	# Parsed with sources kept, so the fixture really does carry a srcref --
	# without that check this would pass vacuously, which is the failure shape
	# this whole file is about.
	f <- eval(parse(
		text = "function(a = 1) {\n\t# NSI_COMMENT_MARKER\n\ta + 1\n}",
		keep.source = TRUE
	))
	expect_false(is.null(attr(body(f), "srcref")))
	# ...and the comment really is recoverable from it, so the absence
	# asserted below is a measurement rather than a tautology. Measured: what
	# strips the comment is deparsing the BODY -- `useSource` reads the
	# FUNCTION's own srcref attribute, so `deparse(f, useSource)` round-trips
	# the comment while `deparse(body(f), useSource)` does not.
	expect_true(any(grepl(
		"NSI_COMMENT_MARKER",
		deparse(f, control = c("keepInteger", "keepNA", "useSource")),
		fixed = TRUE
	)))

	txt <- unlist(lapply(
		.ns_code_exprs(f),
		deparse,
		control = c("keepInteger", "keepNA")
	))
	expect_false(any(grepl("NSI_COMMENT_MARKER", txt, fixed = TRUE)))
	expect_false(any(grepl("srcref", txt, fixed = TRUE)))
	# The formal's default renders as its OWN element, not merely as a
	# substring of the body: `grepl("1", txt)` would also be satisfied by the
	# body's `a + 1` and so would pass with the default dropped.
	expect_true("1" %in% txt)
})

test_that(".ns_deparsed_code covers bodies and defaults for the namespace", {
	fns <- .ns_functions("fetwfe")
	code <- .ns_deparsed_code("fetwfe")
	expect_identical(names(code), names(fns))
	expect_true(all(vapply(code, is.character, logical(1))))
	expect_false(any(vapply(
		code,
		function(x) any(grepl("srcref", x, fixed = TRUE)),
		logical(1)
	)))

	# Every function with a code-carrying default renders MORE lines than its
	# body alone. No name and no count is pinned, so nothing here drifts; what
	# is pinned is that the defaults are still being scanned at all. Reverting
	# `.ns_code_exprs()` to `body()`-only makes every one of these equal
	# instead of greater, and this fails.
	has_code_default <- vapply(
		fns,
		function(f) {
			fmls <- formals(f)
			if (length(fmls) == 0L) {
				return(FALSE)
			}
			any(vapply(
				seq_along(fmls),
				function(i) {
					!identical(fmls[[i]], quote(expr = )) && !is.null(fmls[[i]])
				},
				logical(1)
			))
		},
		logical(1)
	)
	expect_true(any(has_code_default))
	body_lines <- vapply(
		names(fns),
		function(nm) {
			b <- body(fns[[nm]])
			if (is.null(b)) {
				0L
			} else {
				length(deparse(b, control = c("keepInteger", "keepNA")))
			}
		},
		integer(1)
	)
	code_lines <- vapply(code, length, integer(1))
	expect_true(all(
		code_lines[has_code_default] > body_lines[has_code_default]
	))
})

test_that(".ns_deparsed_code does not depend on the keep.source option", {
	# `control` is pinned rather than left to the default, so the rendering
	# does not move with the option in force at call time. That is the
	# in-process half of the property the guardrails rely on; the other half --
	# that the rendering is also the same under `R CMD check`, against an
	# installed package with no srcrefs -- is a cross-runner claim no
	# single-runner test can make, and is measured by running the suite under
	# both.
	old <- options(keep.source = TRUE)
	on.exit(options(old), add = TRUE)
	kept <- .ns_deparsed_code("fetwfe")
	options(keep.source = FALSE)
	dropped <- .ns_deparsed_code("fetwfe")
	expect_identical(kept, dropped)
})

# --- .ns_functions ------------------------------------------------------------

test_that(".ns_functions returns only functions, name-sorted", {
	fns <- .ns_functions()
	expect_true(length(fns) > 0L)
	expect_true(all(vapply(fns, is.function, logical(1))))
	expect_identical(names(fns), sort(names(fns)))
	expect_identical(fns, .ns_functions("fetwfe"))

	ns <- asNamespace("fetwfe")
	# `.packageName` is a character constant every package namespace carries,
	# so this probe is live in any package and cannot drift with this one. It
	# must be in the namespace and NOT in the result.
	expect_true(".packageName" %in% ls(ns, all.names = TRUE))
	expect_false(".packageName" %in% names(fns))
})

test_that(".ns_functions sees dot-prefixed internals, not just exports", {
	# `all.names = TRUE` is load-bearing: every internal this package's
	# guardrails care about starts with a dot, and dropping the argument would
	# make them all invisible while every set assertion still passed.
	fns <- .ns_functions("fetwfe")
	expect_true(".floor_cluster_quad" %in% names(fns))
	expect_true("fetwfe" %in% names(fns))
	expect_identical(
		fns[[".floor_cluster_quad"]],
		fetwfe:::.floor_cluster_quad
	)
})
