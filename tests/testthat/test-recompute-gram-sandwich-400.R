# Direct unit test for the extracted scaffold helper .recompute_gram_and_sandwich()
# (#400 Phase 2a). Its one genuine policy axis -- `on_singular` (degrade vs stop
# on a singular recomputed Gram) -- is DEFENSIVE-ONLY and unreachable through any
# public fit (a fit with valid SEs already inverted its Gram on this support; see
# the note in test-cross-accessor-scaffold-guardrail-400.R). So the asymmetry can
# only be exercised here, by feeding the helper a rank-deficient selected support.

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
	# (pins the helper's verbatim propagation). NB: the live call-site literal in
	# .simultaneous_cis_impl() is a SEPARATE copy of this string, kept byte-identical
	# by hand -- that singular-Gram path is unreachable through any public fit, so no
	# integration test exercises it; edit one copy, edit both.
	msg <- paste0(
		"simultaneousCIs(): the Gram matrix on the selected support is not ",
		"invertible; the analytic method's assumptions are not satisfied. ",
		"For a high-dimensional (p >= NT) design, use method = 'bootstrap'."
	)
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
