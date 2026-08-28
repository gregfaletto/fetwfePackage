# Package-agnostic primitives for guardrail tests that assert something about
# the SHAPE of the package's own source code.
#
# WHY THIS FILE EXISTS. The reflex when writing a "does the source still do X?"
# test is to `readLines()` an `R/*.R` file. That test runs on no automated
# machine: under `R CMD check` -- which is what every CI job and every CRAN
# machine runs -- the tests execute against the *installed* package, which has
# no `R/` directory at all, so the file-existence guard fires and the whole
# block skips. Issue #463 is exactly that failure: a guardrail that had been
# green for months had never actually executed outside a source tree.
#
# The package NAMESPACE, by contrast, exists identically under
# `devtools::load_all()` and under `test_check()` against an installed,
# byte-compiled package: `body()` still returns the AST, and `deparse()` of
# that AST is byte-identical between the two runners. So a guardrail that walks
# the namespace runs everywhere, and it sees every function the package
# defines -- including `@noRd` internals -- rather than whichever files someone
# remembered to list.
#
# testthat sources every `helper-*.R` in `tests/testthat/` once, before any
# `test-*.R` runs, under both runners, so these are in scope for every test
# file with no `source()` call. The `.ns_` prefix keeps them from colliding
# with the fixtures the other `helper-*.R` files put in the same environment.
#
# Everything here is package-agnostic on purpose; anything specific to one
# guardrail's domain belongs in that guardrail's own test file.

# Every function object in `pkg`'s namespace, as a name-sorted named list.
# Non-function objects (the namespace bookkeeping environments, `load_all()`'s
# injected `.__DEVTOOLS__`, package-level constants) are dropped.
.ns_functions <- function(pkg = "fetwfe") {
	ns <- asNamespace(pkg)
	nms <- sort(ls(ns, all.names = TRUE))
	objs <- mget(nms, envir = ns, inherits = FALSE)
	objs[vapply(objs, is.function, logical(1))]
}

# The deparsed body of every namespace function, as a named list of character
# vectors. `control` is pinned rather than left to the default so the rendering
# stays deterministic across R versions; `useSource` is deliberately absent, so
# comments never appear in the output even when a `srcref` is attached (which
# it is under `options(keep.source.pkgs = TRUE)` plus `load_all()`). A function
# with a `NULL` body contributes `character(0)`.
.ns_deparsed_bodies <- function(pkg = "fetwfe") {
	lapply(.ns_functions(pkg), function(f) {
		if (is.null(body(f))) {
			character(0)
		} else {
			deparse(body(f), control = c("keepInteger", "keepNA"))
		}
	})
}

# Depth-first walk over a parsed expression, calling `visit(node, ancestors)`
# at every node. `ancestors` is a list ordered outermost-first and excludes the
# node itself; `visit()`'s return value is ignored, so a visitor accumulates
# with `<<-`.
#
# The child is tested and recursed on INLINE, never through a local. R
# represents both an absent index (`x[, 1]`) and a defaulted-less formal with
# its *empty symbol*, and binding that to a local makes the local behave as a
# missing argument: the next mention of it raises
# `argument "k" is missing, with no default` from a line that looks like an
# ordinary recursive call. Testing `identical(kids[[i]], quote(expr = ))`
# inline is safe; `k <- kids[[i]]` is not.
.ns_walk_ast <- function(expr, visit, ancestors = list()) {
	visit(expr, ancestors)
	if (is.call(expr) || is.pairlist(expr) || is.expression(expr)) {
		kids <- as.list(expr)
		inner <- c(ancestors, list(expr))
		for (i in seq_along(kids)) {
			if (identical(kids[[i]], quote(expr = ))) {
				next
			}
			.ns_walk_ast(kids[[i]], visit, inner)
		}
	}
	invisible(NULL)
}

# Does the subtree rooted at `expr` mention any of the symbols named in `syms`
# (a character vector)?
.ns_subtree_has_symbol <- function(expr, syms) {
	found <- FALSE
	.ns_walk_ast(expr, function(node, ancestors) {
		if (!found && is.symbol(node) && as.character(node) %in% syms) {
			found <<- TRUE
		}
	})
	found
}

# Is `expr` a call to any of the function names in `fn` (a character vector)?
#
# Resolves a namespace-qualified head. The head of `Matrix::crossprod(...)` is
# itself a call -- ``::``(Matrix, crossprod) -- not a symbol, so the natural
# `is.symbol(expr[[1]])` test silently declines to match it, and a
# `pkg::op()` spelling escapes any guardrail built on this predicate. Measured:
# that is exactly how a `Matrix::crossprod()` site slips past a structural
# check while naming its operands literally.
.ns_is_call_to <- function(expr, fn) {
	if (!is.call(expr)) {
		return(FALSE)
	}
	# Element 1 of a call is its function part and is never the empty symbol,
	# so binding it here cannot hit the missing-argument trap above.
	head <- expr[[1]]
	if (
		is.call(head) &&
			length(head) == 3L &&
			is.symbol(head[[1]]) &&
			as.character(head[[1]]) %in% c("::", ":::")
	) {
		head <- head[[3]]
	}
	is.symbol(head) && as.character(head) %in% fn
}
