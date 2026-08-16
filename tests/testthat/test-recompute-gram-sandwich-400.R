# Direct unit test for the extracted scaffold helper .recompute_gram_and_sandwich()
# (#400 Phase 2a). Its one genuine policy axis is `on_singular` (degrade vs stop
# on a singular recomputed Gram). This file exercises that axis SYNTHETICALLY, by
# feeding the helper a rank-deficient selected support directly. For which
# callers can reach the branch through a public door, and how, see the roxygen on
# .recompute_gram_and_sandwich() in R/variance_machinery.R -- it is the canonical
# statement and it names the end-to-end tests.

test_that(".recompute_gram_and_sandwich() degrades vs stops on a singular selected-support Gram (#400)", {
	# Rank-deficient support: feature columns 1 and 2 are identical, so the centered
	# Gram on the selected support {1,2,3} is singular (getGramInv's eigenvalue guard
	# fires and warns).
	N <- 4L
	T <- 2L
	n <- N * T
	x1 <- as.double(seq_len(n))
	X_sing <- cbind(x1, x1, rep(c(1, 0), n / 2))
	y_final <- rep(0, n)
	treat_inds <- c(1L, 2L)
	num_treats <- 2L
	sel_feat_inds <- c(1L, 2L, 3L)
	sel_treat_inds_shifted <- c(1L, 2L)

	# on_singular = "degrade" -> calc_ses = FALSE, no sandwich, no error.
	gs <- suppressWarnings(fetwfe:::.recompute_gram_and_sandwich(
		X_final = X_sing,
		y_final = y_final,
		N = N,
		T = T,
		treat_inds = treat_inds,
		num_treats = num_treats,
		sel_feat_inds = sel_feat_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		se_type = "cluster",
		on_singular = "degrade"
	))
	expect_false(isTRUE(gs$calc_ses))
	expect_null(gs$sandwich_full)
	expect_null(gs$treat_block_mask)

	# on_singular = "stop" -> errors with EXACTLY the stop_message it is given
	# (pins the helper's verbatim propagation). The string below is the very
	# object .simultaneous_cis_impl() passes at its live call site: since #429
	# both read the package-level constant, so the two cannot drift apart.
	msg <- fetwfe:::.SINGULAR_GRAM_ANALYTIC_STOP_MSG
	err <- tryCatch(
		suppressWarnings(fetwfe:::.recompute_gram_and_sandwich(
			X_final = X_sing,
			y_final = y_final,
			N = N,
			T = T,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_feat_inds = sel_feat_inds,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			se_type = "default",
			on_singular = "stop",
			stop_message = msg
		)),
		error = function(e) conditionMessage(e)
	)
	expect_identical(err, msg)

	# Sanity: a well-conditioned support returns calc_ses = TRUE + a real gram_inv,
	# and (cluster) a sandwich + mask -- so the degrade above is genuinely the
	# singular branch, not a helper that always fails.
	X_ok <- cbind(x1, x1^2, rep(c(1, 0), n / 2))
	gs_ok <- fetwfe:::.recompute_gram_and_sandwich(
		X_final = X_ok,
		y_final = as.double(seq_len(n)),
		N = N,
		T = T,
		treat_inds = treat_inds,
		num_treats = num_treats,
		sel_feat_inds = sel_feat_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		se_type = "cluster",
		on_singular = "stop",
		stop_message = msg
	)
	expect_true(isTRUE(gs_ok$calc_ses))
	expect_true(is.matrix(gs_ok$gram_inv))
	expect_false(is.null(gs_ok$sandwich_full))
	expect_length(gs_ok$treat_block_mask, length(sel_feat_inds))
})

# Two properties of the helper's own signature, both of which used to be
# invisible to the whole suite (#429 item 4): reordering the `on_singular`
# default and dropping `call. = FALSE` from its stop() each left every test
# green. This block rebuilds a rank-deficient support of the same shape as the
# one above -- feature columns 1 and 2 identical, so the centered Gram on the
# selected support {1, 2, 3} is singular. It is rebuilt rather than hoisted to
# file scope because `T <- 2L` there would shadow base R's `T` for every other
# block in the file. Nothing enforces that the two stay identical, and nothing
# needs to: each block only requires *a* singular support, so an edit to one
# does not invalidate the other.
test_that(".recompute_gram_and_sandwich() defaults to degrade and errors without a call (#429)", {
	N <- 4L
	T <- 2L
	n <- N * T
	x1 <- as.double(seq_len(n))
	X_sing <- cbind(x1, x1, rep(c(1, 0), n / 2))
	shared <- list(
		X_final = X_sing,
		y_final = rep(0, n),
		N = N,
		T = T,
		treat_inds = c(1L, 2L),
		num_treats = 2L,
		sel_feat_inds = c(1L, 2L, 3L),
		sel_treat_inds_shifted = c(1L, 2L),
		se_type = "default"
	)

	# match.arg() takes the first element, so the default IS the order of the
	# formal. The structural half names the intent...
	expect_identical(
		eval(formals(fetwfe:::.recompute_gram_and_sandwich)$on_singular),
		c("degrade", "stop")
	)
	# ...and the behavioral half is what catches a reorder: called on a singular
	# support WITHOUT on_singular, the helper must degrade rather than error.
	gs_default <- suppressWarnings(do.call(
		fetwfe:::.recompute_gram_and_sandwich,
		shared
	))
	expect_false(isTRUE(gs_default$calc_ses))

	# call. = FALSE on the "stop" branch. Without it R attaches this internal
	# helper's own call, so the user sees
	# `Error in .recompute_gram_and_sandwich(X_final = ..., ...)` with the whole
	# design matrix deparsed into it instead of the actionable message.
	e <- tryCatch(
		suppressWarnings(do.call(
			fetwfe:::.recompute_gram_and_sandwich,
			c(
				shared,
				list(
					on_singular = "stop",
					stop_message = fetwfe:::.SINGULAR_GRAM_ANALYTIC_STOP_MSG
				)
			)
		)),
		error = function(e) e
	)
	expect_s3_class(e, "error")
	expect_null(conditionCall(e))
	# ...and that it is the RIGHT error, not merely a call-less one. Without
	# this the block would pass on any error at all, including one raised
	# before the on_singular branch was reached.
	expect_identical(
		conditionMessage(e),
		fetwfe:::.SINGULAR_GRAM_ANALYTIC_STOP_MSG
	)
})
