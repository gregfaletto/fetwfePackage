# Unit tests for the eight small fixes collected in issue #431 (the late
# post-execution audit of PRs #412-#421). One test_that() block per item, each
# named for its item. Every block was mutation-checked: the source change was
# reverted one item at a time and the named block confirmed to go red.
#
# HOUSE RULE FOR THIS FILE: every expect_error() / grepl() anchor passes
# `fixed = TRUE`. `expect_error`'s `regexp` is a regular expression, and every
# anchor the new guards can offer is a stopifnot() deparse full of parentheses
# and dots -- `stopifnot(is.character(sm), ...)` produces the message
# "is.character(sm) is not TRUE", whose `(sm)` reads as a capture group, so
# grepl(msg, msg) is FALSE while grepl(msg, msg, fixed = TRUE) is TRUE. The
# failure is loud (the anchor will not match its own message), but the natural
# repair under time pressure -- truncating the anchor to "is not TRUE" -- matches
# every stopifnot() failure in the frame and is genuinely blind.

# ------------------------------------------------------------------------------
# Item 1: .validate_gen_coefs_scalars() raises with call. = FALSE, so the error
# no longer prints `Error in .validate_gen_coefs_scalars(T = T, G = G, ...)`
# echoing formal names instead of the user's values.
#
# expect_null(conditionCall(e)) ALONE is a partially blind guardrail: it passes
# for any call-less error, including one raised before the validator is reached,
# and both entry points run .resolve_cohort_count_arg() and
# .validate_targeted_sparsity() first -- the latter already with call. = FALSE.
# (The same warning is recorded at test-recompute-gram-sandwich-400.R:148-151.)
# The cross-door identity and six-distinct-messages assertions restore the
# discrimination WITHOUT writing a single message literal, which keeps #429
# item 9 free to own the literals.
# ------------------------------------------------------------------------------
test_that("genCoefs()/genCoefsCore() scalar errors carry no call (#431 item 1)", {
	grab <- function(expr) tryCatch(expr, error = function(e) e)

	# One case per validator branch, in source order. Every other argument is
	# valid, so exactly the named branch fires.
	cases <- list(
		T_bad = list(G = 3, T = 1, d = 2, density = 0.5, eff_size = 1),
		G_bad = list(G = 0, T = 5, d = 2, density = 0.5, eff_size = 1),
		G_gt_T_minus_1 = list(G = 5, T = 5, d = 2, density = 0.5, eff_size = 1),
		d_bad = list(G = 3, T = 5, d = -1, density = 0.5, eff_size = 1),
		density_bad = list(G = 3, T = 5, d = 2, density = 0, eff_size = 1),
		eff_size_bad = list(G = 3, T = 5, d = 2, density = 0.5, eff_size = "a")
	)

	msgs <- character(0)
	for (nm in names(cases)) {
		e_gc <- grab(do.call(genCoefs, cases[[nm]]))
		e_core <- grab(do.call(genCoefsCore, cases[[nm]]))

		expect_s3_class(e_gc, "error")
		expect_s3_class(e_core, "error")

		# The fix itself: no call is attached to either condition.
		expect_null(conditionCall(e_gc), info = nm)
		expect_null(conditionCall(e_core), info = nm)

		# Both public doors reach the same single-sourced branch (#401), which
		# is what makes the expect_null() above informative rather than
		# satisfiable by any earlier call-less error.
		expect_identical(
			conditionMessage(e_gc),
			conditionMessage(e_core),
			info = nm
		)

		msgs <- c(msgs, conditionMessage(e_gc))
	}

	# All six branches are actually distinct and reachable: if a fixture ever
	# stopped reaching its own branch it would collide with another message here.
	expect_length(msgs, 6L)
	expect_length(unique(msgs), 6L)
})

# ------------------------------------------------------------------------------
# Item 2: .compute_att_pair() calls the treatment-effect function through
# .call_te() (by name, with symbol arguments) rather than do.call() on a
# function OBJECT. do.call(f, args) with `f` a closure builds a call whose first
# element is the entire function definition and whose remaining elements are the
# argument VALUES -- so every stopifnot() condition raised inside `f` captures a
# ~1.28 million character deparse and a ~570 KB condition object, and reports
# `Error in (function (sig_eps_sq, N, T, G, num_treats, ...)` instead of naming
# the function that crashed.
#
# ASSERTION ORDER IS LOAD-BEARING. Under the mutation the call's first element is
# the closure itself, and as.character() on a closure does not return a
# non-matching string -- it RAISES ("cannot coerce type 'closure' to vector of
# type 'character'"). Leading with the size bound and the is.name() check makes
# the mutation register as real FAILURES first rather than only as an error that
# aborts the block. The 2000-character bound is structural (5x headroom over the
# measured 391 on the fixed path, against 1,276,681 on the mutant), not a pinned
# exact length, and nothing stochastic is involved.
#
# EACH FIXTURE MUST REACH A stopifnot() IN THE NAMED FUNCTION'S OWN FRAME, not a
# conformability error on the way in: a fixture that trips a deeper frame names
# the deeper frame (`%*%`, a 23-character deparse), which PASSES the size
# assertion and fails the name assertion for a reason unrelated to the fix. The
# two triggers are different functions and are constructed separately.
# ------------------------------------------------------------------------------
test_that(".compute_att_pair() names the te function and stays small (#431 item 2)", {
	G <- 3L
	num_treats <- 6L
	# Large enough that a do.call()-inlined copy dwarfs the 2000-char bound.
	big_gram_inv <- matrix(0, nrow = 258L, ncol = 258L)

	# --- getTeResultsOLS(): stopifnot(nrow(psi_mat) == length(tes)) at
	# R/variance_machinery.R, inside `if (calc_ses)`. att_hat is computed just
	# above it, so cohort_tes %*% cohort_probs must still conform.
	ols_args <- list(
		sig_eps_sq = 1,
		N = 10L,
		T = 4L,
		G = G,
		num_treats = num_treats,
		cohort_tes = rep(0.1, G),
		psi_mat = matrix(0, nrow = 2L, ncol = G),
		gram_inv = big_gram_inv,
		tes = rep(0, 3), # length 3 != nrow(psi_mat) == 2
		calc_ses = TRUE
	)
	e_ols <- tryCatch(
		fetwfe:::.compute_att_pair(
			te_fn_name = "getTeResultsOLS",
			base_args = ols_args,
			in_sample_probs = list(
				cohort_probs = rep(1 / 3, G),
				cohort_probs_overall = rep(0.1, G)
			),
			indep_probs_args = list(),
			indep_count_data_available = FALSE
		),
		error = function(e) e
	)
	expect_s3_class(e_ols, "error")
	cl_ols <- conditionCall(e_ols)
	expect_lt(nchar(paste(deparse(cl_ols), collapse = "")), 2000)
	expect_true(is.name(cl_ols[[1]]))
	expect_identical(as.character(cl_ols[[1]]), "getTeResultsOLS")

	# --- getTeResults2(): stopifnot(all(!is.na(cohort_probs))), which sits
	# BEFORE calc_ses is consulted and before the first %*%. It needs a
	# one-element NA and no matrix-conformability reasoning at all; building a
	# psi_mat fixture for this one instead lands in .compute_att_var1() or %*%.
	bridge_args <- list(
		sig_eps_sq = 1,
		N = 10L,
		T = 4L,
		G = 1L,
		num_treats = num_treats,
		cohort_tes = 0.1,
		psi_mat = matrix(0, nrow = 1L, ncol = 1L),
		gram_inv = big_gram_inv,
		sel_treat_inds_shifted = 1L,
		d_inv_treat_sel = matrix(0, nrow = 1L, ncol = 1L),
		first_inds = 1L,
		theta_hat_treat_sel = 0,
		calc_ses = TRUE
	)
	e_bridge <- tryCatch(
		fetwfe:::.compute_att_pair(
			te_fn_name = "getTeResults2",
			base_args = bridge_args,
			in_sample_probs = list(
				cohort_probs = NA_real_,
				cohort_probs_overall = 0.1
			),
			indep_probs_args = list(),
			indep_count_data_available = FALSE
		),
		error = function(e) e
	)
	expect_s3_class(e_bridge, "error")
	cl_bridge <- conditionCall(e_bridge)
	expect_lt(nchar(paste(deparse(cl_bridge), collapse = "")), 2000)
	expect_true(is.name(cl_bridge[[1]]))
	expect_identical(as.character(cl_bridge[[1]]), "getTeResults2")
})

# ------------------------------------------------------------------------------
# Item 3: getP() and .collapse_design_for_twfe_covs() now call .base_cols()
# instead of re-deriving the pre-treatment column-count formula inline.
#
# WHAT THIS ADDS. Both fold sites are ALREADY covered by the existing suite --
# folding getP() wrongly produces 5 failures + 15 errors, and dropping
# `+ num_treats` from the input_prep.R fold produces 2 errors. What is new here
# is the exhaustive algebraic identity over the (G, T, d) domain, which nothing
# in the suite provides. This is a CONSOLIDATION-PARITY guardrail: it pins that
# the folded form equals the pre-fold form, not that either matches
# paper_arxiv.tex, and it cannot detect a formula that was always wrong. Given
# this PR's "no math changes" criterion, that is the right instrument.
#
# num_treats MUST BE BOUND AND NON-ZERO. The identity holds for all num_treats
# including 0, but the mutation this block guards (getP() folded as
# `num_treats * d` instead of `num_treats * (1 + d)`) differs from the correct
# value by exactly num_treats -- so at num_treats = 0 there are ZERO visible
# cells over the whole grid and the block stays green under its own mutation.
# getNumTreats(G, T) is positive everywhere on this grid (minimum
# getNumTreats(1, 2) == 1), and it is what arms the mutation.
#
# Loop variables are deliberately NOT named `T`: item 4's fixture below needs
# `T` to resolve to base::T, and although test_that() blocks do not leak, using
# distinct names removes the question entirely.
# ------------------------------------------------------------------------------
test_that(".base_cols() folds are algebraically identical to the inline forms (#431 item 3)", {
	cells <- 0L
	for (t_v in 2:20) {
		for (g_v in 1:(t_v - 1)) {
			for (d_v in 0:8) {
				nt <- fetwfe:::getNumTreats(g_v, t_v)
				expect_gt(nt, 0)

				# Oracle written out longhand; it does NOT call .base_cols(), so
				# the comparison is not tautological.
				oracle_p <- g_v +
					(t_v - 1) +
					d_v +
					d_v * g_v +
					d_v * (t_v - 1) +
					nt +
					nt * d_v
				expect_identical(
					fetwfe:::getP(g_v, t_v, d_v, nt),
					oracle_p,
					info = paste("getP", g_v, t_v, d_v)
				)

				# The .collapse_design_for_twfe_covs() slice bound.
				oracle_slice <- g_v + t_v - 1 + d_v * (1 + g_v + t_v - 1) + nt
				expect_identical(
					fetwfe:::.base_cols(g_v, t_v, d_v) + nt,
					oracle_slice,
					info = paste("slice", g_v, t_v, d_v)
				)

				cells <- cells + 1L
			}
		}
	}
	# The grid is non-empty and the loop actually ran (a silently-empty grid
	# would make every assertion above vacuous).
	expect_identical(cells, 1710L)
})

# ------------------------------------------------------------------------------
# Item 4: .base_cols() and getNumTreats() reject non-scalar / non-numeric
# arguments instead of silently returning a wrong count.
#
# THE RED SIDE, measured before the guards landed (recorded here rather than
# asserted -- after the guard those calls error, so there is no value left to
# compare against, and expect_false(identical(.base_cols(3, T, 2), 11)) would
# itself error rather than assert):
#
#     .base_cols(3, T, 2)      -> 11    (correct answer for T = 5, d = 2 is 23)
#     getNumTreats(3, T)       -> -3
#     .base_cols(3, 5, NULL)   -> numeric(0)
#     .base_cols(c(2, 3), 5, 2)-> c(20, 23)
#
# `T` MUST BE UNBOUND in this block, so that it resolves to base::T (i.e. TRUE)
# and is.numeric(TRUE) is FALSE. The fixture is self-revealing rather than
# fragile: if this file ever binds `T`, .base_cols(3, 5, 2) returns cleanly and
# the expect_error() fails loudly instead of passing for the wrong reason.
#
# The anchors buy discrimination against an UNRELATED error, not against the
# mutation -- a deleted guard returns a value rather than erroring.
# ------------------------------------------------------------------------------
test_that(".base_cols()/getNumTreats() reject bad arguments (#431 item 4)", {
	# base::T -- deliberately not bound in this block.
	expect_error(
		fetwfe:::.base_cols(3, T, 2),
		"is.numeric(T) is not TRUE",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.base_cols(3, 5, NULL),
		"is.numeric(d) is not TRUE",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.base_cols(c(2, 3), 5, 2),
		"length(G) == 1L is not TRUE",
		fixed = TRUE
	)

	expect_error(
		fetwfe:::getNumTreats(3, T),
		"is.numeric(T) is not TRUE",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::getNumTreats(NULL, 5),
		"is.numeric(G) is not TRUE",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::getNumTreats(c(2, 3), 5),
		"length(G) == 1L is not TRUE",
		fixed = TRUE
	)

	# Integer arguments still pass: is.numeric() accepts them, and
	# test-cleanup-batch-185.R:28 passes 3L / 6L. (`.base_cols()` returns a
	# double even on integer input, because the literal `1` in `G + T - 1` is a
	# double; `getNumTreats()` coerces with as.integer() per its @return
	# contract.)
	expect_identical(fetwfe:::.base_cols(3L, 6L, 2L), 26)
	expect_identical(fetwfe:::getNumTreats(3L, 6L), 12L)
})

# ------------------------------------------------------------------------------
# Item 5: .cv_grpreg_at_min()'s nearest-match fallback now asserts that the
# index it picked actually corresponds to cv.grpreg's lambda.min, instead of
# silently returning the nearest one.
#
# Neither case below is reachable from a real cv.grpreg() fit as fetwfe calls it
# -- cv.grpreg sets lambda.min by INDEXING fit$lambda with no arithmetic
# applied, so it is bit-identical to a grid element by construction, and
# setupLambda.gBridge's grid is a strictly monotone exponential sequence with
# distinct == length at every resolution tested (grpreg 3.5.0). Mocks are
# therefore the only way to reach the fallback; do not attempt to provoke it
# from a real fit.
#
# ONLY CASE (a) CARRIES MUTATION POWER. With the new stopifnot() deleted, (a)
# returns 1 silently and reddens. Case (b) returns 2 under both the mutated and
# the unmutated helper, so it is a no-regression pin on first-wins, not a
# mutation-covered assertion.
#
# Grids are written INCREASING to match the real setupLambda.gBridge output
# (all(diff(fit$lambda) > 0)); neither assertion depends on the order.
# ------------------------------------------------------------------------------
test_that(".cv_grpreg_at_min() rejects a lambda.min absent from the grid (#431 item 5)", {
	mk <- function(grid, lmin) {
		list(lambda.min = lmin, fit = list(lambda = grid))
	}

	# (a) lambda.min not in the grid at all -- grpreg would have broken its
	# contract. Returned 1 silently before this fix; must now error.
	expect_error(
		fetwfe:::.cv_grpreg_at_min(mk(c(0.1, 0.25, 0.5, 1), 99)),
		"isTRUE(all.equal(cv_fit$fit$lambda[lam_idx], lambda_star)) is not TRUE",
		fixed = TRUE
	)

	# (b) a duplicated grid value equal to lambda.min: which() returns length 2,
	# the fallback fires, which.min() takes the FIRST (which selects the
	# earliest column of fit$beta, since duplicates are equal lambdas and
	# "smallest lambda" would say nothing), and the new assertion passes.
	expect_identical(
		fetwfe:::.cv_grpreg_at_min(mk(c(0.1, 0.5, 0.5, 1), 0.5)),
		2L
	)

	# Sanity: the ordinary exact-match path is untouched and never reaches the
	# fallback.
	expect_identical(
		fetwfe:::.cv_grpreg_at_min(mk(c(0.1, 0.25, 0.5, 1), 0.5)),
		3L
	)
})

# ------------------------------------------------------------------------------
# Item 6: .check_c6_dims_toplevel() reads every slot with [[ rather than $, so a
# malformed object can no longer partial-match its way to a passing contract.
#
# FOUR ONE-SIBLING FIXTURES, ONE PER SLOT READ -- not one fixture in which all
# four reads diverge. .assert_contract() (R/class_helpers.R) stop()s on the
# FIRST violated contract, and this helper checks beta_hat first, so on an
# all-diverging fixture the [[ version always dies on contract 1 and reverting
# `y`, `X_ints`, or `calc_ses` to `$` raises a byte-identical message. Measured:
# 0 of 4 mutations detected with a bare expect_error(), 1 of 4 with an anchored
# one, and no assertion over that single fixture can observe the other three.
#
# Each fixture renames exactly one slot to a partial-match sibling, so that read
# has a prefix sibling and no exact slot while the other three stay exact and
# correct. Each is anchored on ITS OWN contract name, so a mutation that merely
# shifts which contract fires is still caught; `info = s` makes a red loop name
# the fixture that failed.
#
# N, T and p MUST STAY BOUND INSIDE THIS BLOCK. Item 4's fixture above needs `T`
# UNBOUND so that .base_cols(3, T, 2) resolves to base::T. test_that() blocks do
# not leak, so a block-local T <- 3L here is harmless -- but hoisted to file
# scope it would make .base_cols(3, 5, 2) return cleanly and item 4's
# expect_error() fail. The repo already carries this exact caution at
# test-recompute-gram-sandwich-400.R:93-94.
#
# RECORDED GAP: these fixtures pin the four SLOT reads. The p / N / T
# conversions in the same function are not mutation-covered by any test; their
# safety rests on the live-object slot inventory plus a green suite.
# ------------------------------------------------------------------------------
test_that(".check_c6_dims_toplevel() reads slots with [[ not $ (#431 item 6)", {
	N <- 4L
	T <- 3L
	p <- 5L
	base <- list(
		beta_hat = rnorm(p),
		p = p,
		y = rnorm(N * T),
		X_ints = matrix(0, N * T, p),
		calc_ses = TRUE,
		N = N,
		T = T
	)
	sib <- c(
		beta_hat = "beta_hat_raw",
		y = "y_mean",
		X_ints = "X_ints_full",
		calc_ses = "calc_ses_flag"
	)
	anchor <- c(
		beta_hat = "C6 length(beta_hat) == p",
		y = "C6 length(y) == N * T",
		X_ints = "C6 nrow(X_ints) == N * T",
		calc_ses = "C8 calc_ses is length-1 logical"
	)

	# Sanity: the well-formed object passes every contract, so the four errors
	# below are caused by the renamed slot and nothing else.
	expect_silent(fetwfe:::.check_c6_dims_toplevel(base, "etwfe"))

	for (s in names(sib)) {
		fx <- base
		v <- fx[[s]]
		fx[[s]] <- NULL
		fx[[sib[[s]]]] <- v
		expect_error(
			fetwfe:::.check_c6_dims_toplevel(fx, "etwfe"),
			anchor[[s]],
			fixed = TRUE,
			info = s
		)
	}
})

# ------------------------------------------------------------------------------
# Item 7: .recompute_gram_and_sandwich() rejects an unusable `stop_message` at
# ENTRY when on_singular = "stop", rather than raising a bare `Error:` with no
# text once the Gram turns out to be singular.
#
# THE WELL-CONDITIONED FIXTURE IS LOAD-BEARING. With the SINGULAR support and
# the guard deleted, stop(NULL, call. = FALSE) still raises and expect_error()
# passes, so the mutation goes undetected (expect_error(stop(NULL, call. =
# FALSE)) and expect_error(stop(c("a", "b"), call. = FALSE)) both pass,
# measured). With the well-conditioned Gram and the guard deleted, all four
# calls RETURN NORMALLY (calc_ses = TRUE), which is observable.
#
# CONSEQUENCE OF ENTRY PLACEMENT: a caller passing on_singular = "stop" with a
# bad stop_message and a WELL-CONDITIONED Gram used to return normally and now
# errors. No live caller does this -- all five call sites were enumerated, and
# only R/simultaneous_cis.R passes "stop", with a valid package-level constant.
# ------------------------------------------------------------------------------
test_that(".recompute_gram_and_sandwich() guards stop_message at entry (#431 item 7)", {
	N <- 4L
	T <- 2L
	n <- N * T
	x1 <- as.double(seq_len(n))
	# Well-conditioned support: x1 and x1^2 are not collinear, so the centered
	# Gram on {1, 2, 3} is invertible and calc_ses comes back TRUE.
	X_ok <- cbind(x1, x1^2, rep(c(1, 0), n / 2))

	call_it <- function(sm) {
		fetwfe:::.recompute_gram_and_sandwich(
			X_final = X_ok,
			y_final = as.double(seq_len(n)),
			N = N,
			T = T,
			treat_inds = c(1L, 2L),
			num_treats = 2L,
			sel_feat_inds = c(1L, 2L, 3L),
			sel_treat_inds_shifted = c(1L, 2L),
			se_type = "default",
			on_singular = "stop",
			stop_message = sm
		)
	}

	expect_error(
		call_it(NULL),
		"is.character(stop_message) is not TRUE",
		fixed = TRUE
	)
	expect_error(call_it(""), "nzchar(stop_message) is not TRUE", fixed = TRUE)
	expect_error(
		call_it(NA_character_),
		"!is.na(stop_message) is not TRUE",
		fixed = TRUE
	)
	expect_error(
		call_it(c("a", "b")),
		"length(stop_message) == 1L is not TRUE",
		fixed = TRUE
	)

	# A valid message still returns normally on this well-conditioned support --
	# so the four errors above are the guard, not a helper that always fails.
	ok <- call_it("a real message")
	expect_true(isTRUE(ok$calc_ses))
})

# ------------------------------------------------------------------------------
# Item 8: the singular-Gram "standard errors will not be calculated" message is
# a single package-level constant read by getGramInv()'s two warning() calls and
# by .compute_cluster_robust_sandwich()'s stop().
#
# BE HONEST ABOUT WHAT EACH ASSERTION BUYS. The expect_identical() assertions are
# tautologies with respect to the message's CONTENT -- after the fold, both
# sides read the same symbol, so they move together when the text is edited.
# They exist to catch a REINTRODUCED INLINE LITERAL THAT DIFFERS, not a rewrite
# of the constant. Only the two grepl(..., fixed = TRUE) anchors hold the
# wording. Same caveat as #429's constant (see R/simultaneous_cis.R).
# ------------------------------------------------------------------------------
test_that("the singular-Gram no-SE message is single-sourced (#431 item 8)", {
	msg <- fetwfe:::.SINGULAR_GRAM_NO_SE_MSG

	# The only assertions in this block that hold the wording itself.
	expect_true(grepl("not invertible", msg, fixed = TRUE))
	expect_true(grepl(
		"Standard errors will not be calculated",
		msg,
		fixed = TRUE
	))

	# getGramInv()'s warning path: a rank-deficient selected support (columns 1
	# and 2 identical) trips the minimum-eigenvalue branch, which degrades to
	# calc_ses = FALSE and warns with the constant.
	N <- 4L
	T <- 2L
	n <- N * T
	x1 <- as.double(seq_len(n))
	X_sing <- cbind(x1, x1, rep(c(1, 0), n / 2))
	w <- tryCatch(
		fetwfe:::getGramInv(
			N = N,
			T = T,
			X_final = X_sing,
			treat_inds = c(1L, 2L),
			num_treats = 2L,
			calc_ses = TRUE,
			sel_feat_inds = c(1L, 2L, 3L),
			sel_treat_inds_shifted = c(1L, 2L)
		),
		warning = function(w) w
	)
	expect_s3_class(w, "warning")
	expect_identical(conditionMessage(w), msg)
})
