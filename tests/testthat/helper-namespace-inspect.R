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

# Every code expression belonging to a function: its body, plus each formal's
# default expression.
#
# `body()` ALONE IS NOT THE FUNCTION'S CODE. A formal default is arbitrary R --
# `f <- function(x, S, v = max(as.numeric(t(x) %*% S %*% x), 0))` computes a
# whole quadratic form and `body(f)` cannot see a character of it. A guardrail
# that read the file as raw text saw defaults for free; one that reads `body()`
# does not, and that difference is a real loss of coverage rather than a
# theoretical one -- it was measured against this package's previous
# source-text guardrail, which caught exactly this mutation while a
# `body()`-only walk passed it.
#
# A formal with no default is the empty symbol; it is tested INLINE, never
# bound to a local, for the reason `.ns_walk_ast()` documents below.
#
# One quirk, measured rather than intended: a formal whose default is the
# literal `NULL` is DROPPED, because `out[[length(out) + 1L]] <- NULL` at an
# index one past the end silently does nothing rather than appending. Harmless -- a `NULL` default has no
# subexpression any guardrail could scan -- and pinned by a test in
# `test-namespace-inspect-463.R` so it cannot change unnoticed. Written down
# because the line reads as an append and is not one.
.ns_code_exprs <- function(f) {
	out <- list()
	fn_body <- body(f)
	if (!is.null(fn_body)) {
		out[[length(out) + 1L]] <- fn_body
	}
	fmls <- formals(f)
	for (i in seq_along(fmls)) {
		if (identical(fmls[[i]], quote(expr = ))) {
			next
		}
		out[[length(out) + 1L]] <- fmls[[i]]
	}
	out
}

# The deparsed code of every namespace function -- body AND formal defaults,
# per `.ns_code_exprs()` -- as a named list of character vectors. `control` is
# pinned rather than left to the default so the rendering stays deterministic
# across R versions rather than tracking whatever the default control set
# becomes.
#
# Comments never appear in the output, but NOT because `useSource` is absent
# from `control` -- an earlier draft of this comment said so, and it is wrong.
# Measured: adding `useSource` here changes nothing, on any of this package's
# functions. `useSource` reproduces a FUNCTION's source ref, and everything
# deparsed below is a body or a default expression, never a function object;
# `deparse(f, ..., "useSource")` does leak comments, `deparse(body(f), ...)`
# does not. The guarantee comes from WHAT is deparsed, not from the control
# set. A function with no body and no defaulted formal contributes
# `character(0)`.
.ns_deparsed_code <- function(pkg = "fetwfe") {
	lapply(.ns_functions(pkg), function(f) {
		exprs <- .ns_code_exprs(f)
		if (length(exprs) == 0L) {
			return(character(0))
		}
		unlist(lapply(
			exprs,
			deparse,
			control = c("keepInteger", "keepNA")
		))
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

# The name of the function `expr` invokes, with any `::` / `:::` qualifier
# stripped, or `NA_character_` when there is no name to return: `expr` is not a
# call at all, or -- after the qualifier is stripped -- its head is not a symbol
# (`f()()`, or a string head).
#
# The unwrapping is the reason this is a primitive rather than an inline
# `expr[[1]]`. The head of `Matrix::crossprod(...)` is itself a call --
# ``::``(Matrix, crossprod) -- not a symbol, so the natural
# `is.symbol(expr[[1]])` test silently declines to match it, and a `pkg::op()`
# spelling escapes any guardrail built on it. Measured: that is exactly how a
# `Matrix::crossprod()` site slips past a structural check while naming its
# operands literally.
#
# Both consumers of that unwrapping go through here -- `.ns_is_call_to()` below
# and `test-cluster_floor.R`'s `.clf_callee_name()`. Adding a third copy is the
# thing this exists to stop, so extend it here rather than re-deriving it.
.ns_callee_name <- function(expr) {
	if (!is.call(expr)) {
		return(NA_character_)
	}
	# Element 1 of a call is its function part and is never the empty symbol,
	# so binding it here cannot hit the missing-argument trap `.ns_walk_ast()`
	# documents above. Named `fn_head` rather than `head` so it does not shadow
	# `utils::head()`.
	fn_head <- expr[[1]]
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

# Is `expr` a call to any of the function names in `fn` (a character vector)?
#
# A namespace-qualified head resolves to the FUNCTION's name, not the package's,
# because `.ns_callee_name()` above owns the unwrapping and the reason for it.
# A non-call, and a call whose head has no name, are both FALSE.
.ns_is_call_to <- function(expr, fn) {
	nm <- .ns_callee_name(expr)
	!is.na(nm) && nm %in% fn
}
