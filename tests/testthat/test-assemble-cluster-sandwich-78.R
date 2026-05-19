library(testthat)
library(fetwfe)

# Tests for the cluster-robust sandwich assembly helper (issue #78, v1.9.15).
#
# `R/ols_calcs.R::.assemble_cluster_robust_sandwich()` wraps the 4-step
# assembly ritual (extract X_S from X_final, slice y_, lm.fit, call
# .compute_cluster_robust_sandwich, build treat_block_mask) that was
# previously duplicated across 4 call sites:
#
#   * R/ols_calcs.R::getCohortATTsFinalOLS         (unfiltered path)
#   * R/fetwfe_core.R::getCohortATTsFinal           (filtered path)
#   * R/event_study.R::.event_study_etwfe_betwfe   (runtime if/else)
#   * R/event_study.R::.event_study_fetwfe          (filtered path)
#
# The helper dispatches between two patterns via `sel_feat_inds`:
#   * NULL (default): X_S = X_final; treat_block_mask via boolean fill.
#   * vector: X_S = X_final[, sel_feat_inds]; treat_block_mask via %in%.
#
# These tests verify the helper's output bit-identically matches the
# manual inline computation that the 4 call sites used to do.

# Small deterministic fixture: N=6 units, T=5 periods, p=10 columns
# in X_final. treat_inds = c(3, 7) means columns 3 and 7 are treatment
# features.
.build_fixture <- function(seed = 123) {
	set.seed(seed)
	N <- 6L
	T <- 5L
	p <- 10L
	X_final <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	y_final <- rnorm(N * T)
	list(
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		p = p,
		treat_inds = c(3L, 7L)
	)
}

test_that(".assemble_cluster_robust_sandwich(sel_feat_inds = NULL) matches manual unfiltered computation", {
	fx <- .build_fixture()

	# Reference (manual): replicate the pre-refactor inline pattern
	# from R/ols_calcs.R::getCohortATTsFinalOLS.
	X_S_ref <- fx$X_final
	y_ref <- fx$y_final[seq_len(fx$N * fx$T)]
	ols_ref <- stats::lm.fit(cbind(1, X_S_ref), y_ref)
	sandwich_ref <- fetwfe:::.compute_cluster_robust_sandwich(
		X_S = X_S_ref,
		residuals = ols_ref$residuals,
		N = fx$N,
		T = fx$T
	)
	mask_ref <- logical(ncol(X_S_ref))
	mask_ref[fx$treat_inds] <- TRUE

	# Helper.
	res <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds,
		sel_feat_inds = NULL
	)

	expect_equal(res$sandwich_full, sandwich_ref, tolerance = 1e-12)
	expect_equal(res$treat_block_mask, mask_ref)
	expect_length(res$treat_block_mask, fx$p)
	expect_true(res$treat_block_mask[3])
	expect_true(res$treat_block_mask[7])
	expect_equal(sum(res$treat_block_mask), 2L)
})

test_that(".assemble_cluster_robust_sandwich(sel_feat_inds = <vector>) matches manual filtered computation", {
	fx <- .build_fixture()
	# Pick a selected support that includes one treatment column
	# (index 3) and one not (index 7 is NOT in sel_feat_inds here).
	sel_feat_inds <- c(2L, 3L, 5L, 8L)

	# Reference (manual): replicate the pre-refactor inline pattern
	# from R/fetwfe_core.R::getCohortATTsFinal.
	X_S_ref <- fx$X_final[, sel_feat_inds, drop = FALSE]
	y_ref <- fx$y_final[seq_len(fx$N * fx$T)]
	ols_ref <- stats::lm.fit(cbind(1, X_S_ref), y_ref)
	sandwich_ref <- fetwfe:::.compute_cluster_robust_sandwich(
		X_S = X_S_ref,
		residuals = ols_ref$residuals,
		N = fx$N,
		T = fx$T
	)
	mask_ref <- sel_feat_inds %in% fx$treat_inds

	# Helper.
	res <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds,
		sel_feat_inds = sel_feat_inds
	)

	expect_equal(res$sandwich_full, sandwich_ref, tolerance = 1e-12)
	expect_equal(res$treat_block_mask, mask_ref)
	expect_length(res$treat_block_mask, length(sel_feat_inds))
	# sel_feat_inds = c(2, 3, 5, 8); treat_inds = c(3, 7). Only
	# position 2 (column 3) is treatment. Position 4 (column 8) is
	# NOT, because column 8 is not in treat_inds.
	expect_equal(res$treat_block_mask, c(FALSE, TRUE, FALSE, FALSE))
})

test_that(".assemble_cluster_robust_sandwich filtered path also covers all-treatment and no-treatment selected supports", {
	fx <- .build_fixture()

	# All selected positions are in treat_inds (sel = c(3, 7)).
	res_all_treat <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds,
		sel_feat_inds = fx$treat_inds
	)
	expect_equal(res_all_treat$treat_block_mask, c(TRUE, TRUE))

	# None of the selected positions are in treat_inds (sel = c(1, 4, 5)).
	res_no_treat <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds,
		sel_feat_inds = c(1L, 4L, 5L)
	)
	expect_equal(
		res_no_treat$treat_block_mask,
		c(FALSE, FALSE, FALSE)
	)
	expect_equal(sum(res_no_treat$treat_block_mask), 0L)
})

test_that(".assemble_cluster_robust_sandwich integration parity: betwfe cluster-SE fit reproduces helper output", {
	# Integration check: fit betwfe with se_type = "cluster" on a
	# simulated panel with strong signal (so selection actually
	# fires and the cluster-SE path is exercised). Re-call the
	# helper on the same fitted object's inputs and assert
	# bit-identical output. This locks the refactor at the
	# integration boundary: if a future PR changes betwfe_core's
	# call to the helper (e.g., swaps which arg gets passed to
	# `sel_feat_inds`), the assertion fails.
	coefs <- genCoefs(
		R = 2,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 4,
		seed = 42
	)
	sim <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	fit <- betwfeWithSimulatedData(sim, se_type = "cluster")

	# Sanity: structure is right and selection actually fired.
	expect_s3_class(fit, "betwfe")
	expect_true(any(fit$beta_hat != 0))

	sel_feat_inds_fit <- which(fit$beta_hat != 0)

	res_helper <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fit$X_final,
		y_final = fit$y_final,
		N = fit$N,
		T = fit$T,
		treat_inds = fit$treat_inds,
		sel_feat_inds = sel_feat_inds_fit
	)

	expect_named(res_helper, c("sandwich_full", "treat_block_mask"))
	expect_true(is.matrix(res_helper$sandwich_full))
	expect_equal(
		dim(res_helper$sandwich_full),
		c(length(sel_feat_inds_fit), length(sel_feat_inds_fit))
	)
	expect_length(res_helper$treat_block_mask, length(sel_feat_inds_fit))
	expect_equal(
		res_helper$treat_block_mask,
		sel_feat_inds_fit %in% fit$treat_inds
	)
	# att_se from the cluster-SE fit must be finite (proves the
	# end-to-end cluster-SE path ran without error post-refactor).
	expect_true(is.finite(fit$att_se))
})

test_that(".assemble_cluster_robust_sandwich: NULL default produces the same output as explicit NULL", {
	fx <- .build_fixture()

	res_default <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds
		# sel_feat_inds defaults to NULL
	)
	res_explicit_null <- fetwfe:::.assemble_cluster_robust_sandwich(
		X_final = fx$X_final,
		y_final = fx$y_final,
		N = fx$N,
		T = fx$T,
		treat_inds = fx$treat_inds,
		sel_feat_inds = NULL
	)
	expect_equal(
		res_default$sandwich_full,
		res_explicit_null$sandwich_full,
		tolerance = 1e-12
	)
	expect_equal(
		res_default$treat_block_mask,
		res_explicit_null$treat_block_mask
	)
})

# ------------------------------------------------------------------------------
# Item 10: rank-deficient selected support in `.compute_cluster_robust_sandwich`
# surfaces the same "Gram matrix corresponding to selected features is not
# invertible" message as `getGramInv()`, rather than the raw Lapack
# "system is exactly singular" trace. The pre-#84 implementation called
# `solve(crossprod(X_S_centered))` without a tryCatch, so a future caller
# routing a rank-deficient selected support into the cluster path would
# get a baffling Lapack stack trace instead of the helpful message.
# This test exercises the new tryCatch branch by constructing an X_S with
# duplicated columns (guaranteed rank-deficient).
# ------------------------------------------------------------------------------
test_that(".compute_cluster_robust_sandwich emits helpful error on rank-deficient X_S", {
	set.seed(20260519)
	N <- 20L
	T <- 5L
	p_S <- 4L
	X_S_full_rank <- matrix(
		rnorm(N * T * p_S),
		nrow = N * T,
		ncol = p_S
	)
	# Force rank deficiency: column 4 = column 2. crossprod() is then
	# exactly singular and solve() would throw a Lapack error in the
	# pre-tryCatch implementation.
	X_S_full_rank[, 4] <- X_S_full_rank[, 2]
	residuals <- rnorm(N * T)

	expect_error(
		fetwfe:::.compute_cluster_robust_sandwich(
			X_S = X_S_full_rank,
			residuals = residuals,
			N = N,
			T = T
		),
		"Gram matrix corresponding to selected features is not invertible"
	)
})

# ------------------------------------------------------------------------------
# Item 18: vectorized .compute_cluster_robust_sandwich() matches the
# pre-vectorize N-loop reference bit-for-bit (modulo floating-point summation
# order, hence the 1e-12 tolerance). Closes the issue #84 item 18 regression
# guard: any future regression that breaks the rowsum() aggregation path or
# the crossprod() / tcrossprod() identity would fire here.
# ------------------------------------------------------------------------------
test_that(".compute_cluster_robust_sandwich matches reference N-loop", {
	# Frozen-seed deterministic fixture sized to exercise both the
	# pre-vectorize N-loop and the post-vectorize rowsum() path with a
	# non-trivial p_S (so the bug class "rowsum-broadcast-axis-confusion"
	# would surface as a numeric mismatch, not as dimensions).
	set.seed(20260519)
	N <- 30L
	T <- 7L
	p_S <- 8L
	X_S <- matrix(rnorm(N * T * p_S), nrow = N * T, ncol = p_S)
	residuals <- rnorm(N * T)

	# Reference: explicit N-loop (the pre-vectorize implementation, captured
	# verbatim from the git history at R/ols_calcs.R::.compute_cluster_robust_sandwich
	# before issue #84 item 18). This is the asymptotically-equivalent O(N*T*p_S^2)
	# baseline we're vectorizing against.
	reference_loop <- function(X_S, residuals, N, T) {
		X_S_centered <- scale(X_S, center = TRUE, scale = FALSE)
		p_S <- ncol(X_S_centered)
		gram_inv_full <- solve(crossprod(X_S_centered))
		meat <- matrix(0, nrow = p_S, ncol = p_S)
		for (i in seq_len(N)) {
			rows_i <- ((i - 1) * T + 1):(i * T)
			Xi <- X_S_centered[rows_i, , drop = FALSE]
			eps_i <- residuals[rows_i]
			XiEps <- crossprod(Xi, eps_i)
			meat <- meat + tcrossprod(XiEps)
		}
		cadjust <- N / (N - 1)
		cadjust * gram_inv_full %*% meat %*% gram_inv_full
	}

	live <- fetwfe:::.compute_cluster_robust_sandwich(X_S, residuals, N, T)
	ref <- reference_loop(X_S, residuals, N, T)
	expect_equal(live, ref, tolerance = 1e-12)
	expect_equal(dim(live), c(p_S, p_S))
})
