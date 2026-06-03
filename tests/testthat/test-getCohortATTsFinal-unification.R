library(testthat)
library(fetwfe)

# Unit tests for the PR #118 unification of getCohortATTsFinal() and the
# (now deleted) getCohortATTsFinalOLS(). The four call paths through the
# unified function correspond to the four `fused × include_selected`
# quadrants:
#
#   fused = FALSE, include_selected = FALSE  →  ETWFE / twfeCovs
#   fused = TRUE,  include_selected = TRUE   →  FETWFE
#   fused = FALSE, include_selected = TRUE   →  BETWFE
#   fused = TRUE,  include_selected = FALSE  →  no real caller; API lockdown
#
# Test 1 freezes the legacy OLS implementation as an inline closure so the
# refactored function can be asserted byte-identical against it. The closure
# is the pre-#118 source of getCohortATTsFinalOLS() copied here verbatim
# (with comments preserved). Future maintainers can re-execute the closure
# against the unified function to confirm the byte-identity claim still
# holds even after later refactors.

# ---- Fixture --------------------------------------------------------------

.build_unif_fixture <- function(seed = 118L) {
	set.seed(seed)
	N <- 8L
	T <- 4L
	G <- 2L
	d <- 2L
	num_treats <- 3L + 2L
	p <- (G + T - 1L + d) + num_treats + 4L
	X_final <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	y_final <- rnorm(N * T)
	treat_inds <- seq.int(G + T - 1L + d + 1L, G + T - 1L + d + num_treats)
	first_inds <- c(1L, 4L)
	c_names <- c("Cohort1", "Cohort2")
	tes <- rnorm(num_treats)
	list(
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		G = G,
		p = p,
		num_treats = num_treats,
		treat_inds = treat_inds,
		first_inds = first_inds,
		c_names = c_names,
		tes = tes,
		sig_eps_sq = 0.7
	)
}

# ---- Test 1: frozen-snapshot equivalence against legacy OLS ---------------

# Inline copy of the pre-#118 getCohortATTsFinalOLS() body. Kept here so the
# byte-equivalence assertion is reproducible without git archaeology.
.legacy_getCohortATTsFinalOLS <- function(
	X_final,
	treat_inds,
	num_treats,
	first_inds,
	c_names,
	tes,
	sig_eps_sq,
	G,
	N,
	T,
	p,
	alpha = 0.05,
	se_type = "default",
	y_final = NULL
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

	stopifnot(length(tes) <= num_treats)
	stopifnot(all(!is.na(tes)))
	stopifnot(nrow(X_final) == N * T)
	stopifnot(ncol(X_final) == p)
	X_to_pass <- X_final

	res <- fetwfe:::getGramInv(
		N = N,
		T = T,
		X_final = X_to_pass,
		treat_inds = treat_inds,
		num_treats = num_treats,
		calc_ses = TRUE
	)
	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses

	sandwich_full <- NULL
	treat_block_mask <- NULL

	if (identical(se_type, "cluster") && calc_ses) {
		stopifnot(!is.null(y_final))
		stopifnot(length(y_final) >= N * T)
		res <- fetwfe:::.assemble_cluster_robust_sandwich(
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			treat_inds = treat_inds,
			sel_feat_inds = NULL
		)
		sandwich_full <- res$sandwich_full
		treat_block_mask <- res$treat_block_mask
	}

	cohort_tes <- rep(as.numeric(NA), G)
	cohort_te_ses <- rep(as.numeric(NA), G)
	psi_mat <- matrix(0, num_treats, G)

	for (g in 1:G) {
		inds_g <- fetwfe:::.cohort_block_inds(g, G, first_inds, num_treats)
		first_ind_g <- inds_g[1]
		last_ind_g <- inds_g[length(inds_g)]

		stopifnot(last_ind_g >= first_ind_g)
		stopifnot(all(first_ind_g:last_ind_g %in% 1:num_treats))

		cohort_tes[g] <- mean(tes[first_ind_g:last_ind_g])

		psi_g <- fetwfe:::getPsiGUnfused(
			first_ind_g,
			last_ind_g,
			sel_treat_inds_shifted = 1:num_treats
		)
		stopifnot(length(psi_g) == num_treats)
		psi_mat[, g] <- psi_g

		if (calc_ses) {
			if (identical(se_type, "cluster")) {
				psi_g_full <- numeric(length(treat_block_mask))
				psi_g_full[treat_block_mask] <- psi_g
				cohort_te_ses[g] <- sqrt(max(
					as.numeric(
						t(psi_g_full) %*%
							sandwich_full %*%
							psi_g_full
					),
					0
				))
			} else {
				cohort_te_ses[g] <- sqrt(
					sig_eps_sq *
						as.numeric(t(psi_g) %*% gram_inv %*% psi_g) /
						(N * T)
				)
			}
		}
	}

	stopifnot(length(c_names) == G)
	stopifnot(length(cohort_tes) == G)

	if (all(!is.na(gram_inv))) {
		stopifnot(length(cohort_te_ses) == G)
		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			cohort_te_ses,
			cohort_tes - stats::qnorm(1 - alpha / 2) * cohort_te_ses,
			cohort_tes + stats::qnorm(1 - alpha / 2) * cohort_te_ses,
			fetwfe:::.compute_p_values(cohort_tes, cohort_te_ses)
		)
		names(cohort_te_ses) <- c_names
		names(cohort_tes) <- c_names
	} else {
		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			rep(NA, G),
			rep(NA, G),
			rep(NA, G),
			rep(NA_real_, G)
		)
	}

	colnames(cohort_te_df) <- c(
		"cohort",
		"estimate",
		"se",
		"ci_low",
		"ci_high",
		"p_value"
	)
	class(cohort_te_df) <- c("catt_df", "data.frame")

	stopifnot(length(tes) == nrow(psi_mat))

	list(
		cohort_te_df = cohort_te_df,
		cohort_tes = cohort_tes,
		cohort_te_ses = cohort_te_ses,
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		calc_ses = calc_ses,
		sandwich_full = sandwich_full,
		treat_block_mask = treat_block_mask
	)
}

test_that("unified getCohortATTsFinal(fused=FALSE, include_selected=FALSE) matches the pre-#118 OLS implementation byte-for-byte", {
	fx <- .build_unif_fixture()

	# Reference: legacy OLS closure baked into this file.
	expected <- .legacy_getCohortATTsFinalOLS(
		X_final = fx$X_final,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		p = fx$p
	)

	# Unified function on the same inputs.
	actual <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = NULL,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = FALSE
	)

	# Byte-identity on every observable slot. d_inv_treat_sel is absent in
	# both (legacy OLS never produces it; unified function omits it under
	# fused = FALSE).
	expect_equal(actual$cohort_te_df, expected$cohort_te_df, tolerance = 0)
	expect_equal(actual$cohort_tes, expected$cohort_tes, tolerance = 0)
	expect_equal(actual$cohort_te_ses, expected$cohort_te_ses, tolerance = 0)
	expect_equal(actual$psi_mat, expected$psi_mat, tolerance = 0)
	expect_equal(actual$gram_inv, expected$gram_inv, tolerance = 0)
	expect_equal(actual$calc_ses, expected$calc_ses)
	expect_null(actual$sandwich_full)
	expect_null(expected$sandwich_full)
	expect_null(actual$treat_block_mask)
	expect_null(expected$treat_block_mask)
	expect_null(actual$d_inv_treat_sel)
})

test_that("unified getCohortATTsFinal(fused=FALSE, include_selected=FALSE) matches the pre-#118 OLS implementation under se_type='cluster'", {
	fx <- .build_unif_fixture()

	expected <- .legacy_getCohortATTsFinalOLS(
		X_final = fx$X_final,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		p = fx$p,
		se_type = "cluster",
		y_final = fx$y_final
	)

	actual <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = NULL,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = FALSE,
		se_type = "cluster",
		y_final = fx$y_final
	)

	expect_equal(actual$cohort_te_df, expected$cohort_te_df, tolerance = 0)
	expect_equal(actual$cohort_te_ses, expected$cohort_te_ses, tolerance = 0)
	expect_equal(actual$sandwich_full, expected$sandwich_full, tolerance = 0)
	expect_equal(actual$treat_block_mask, expected$treat_block_mask)
})

# ---- Test 2: include_selected = TRUE appends a `selected` column ----------

test_that("getCohortATTsFinal(fused=TRUE, include_selected=TRUE) attaches a `selected` column matching cohort_tes != 0 (FETWFE)", {
	fx <- .build_unif_fixture()
	# Mimic FETWFE's call shape: sel_feat_inds spans all features (no
	# bridge zeroing on this fixture), sel_treat_inds_shifted spans
	# all treatment effects.
	out <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = seq_len(fx$p),
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = TRUE,
		calc_ses = TRUE,
		include_selected = TRUE
	)

	expect_true("selected" %in% colnames(out$cohort_te_df))
	# cohort_tes is named (the function attaches c_names); the data.frame
	# column carries the same logicals without names. Strip names before
	# comparing.
	expect_identical(out$cohort_te_df$selected, unname(out$cohort_tes != 0))
	expect_true(!is.null(out$d_inv_treat_sel))
	expect_true(is.matrix(out$d_inv_treat_sel))
	# psi_mat dimensions track sel_treat_inds_shifted, not num_treats.
	expect_equal(nrow(out$psi_mat), fx$num_treats)
	expect_equal(ncol(out$psi_mat), fx$G)
})

test_that("getCohortATTsFinal(fused=FALSE, include_selected=TRUE) attaches a `selected` column (BETWFE) but no d_inv_treat_sel", {
	fx <- .build_unif_fixture()
	out <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = seq_len(fx$p),
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = TRUE
	)

	expect_true("selected" %in% colnames(out$cohort_te_df))
	# cohort_tes is named (the function attaches c_names); the data.frame
	# column carries the same logicals without names. Strip names before
	# comparing.
	expect_identical(out$cohort_te_df$selected, unname(out$cohort_tes != 0))
	expect_null(out$d_inv_treat_sel)
})

# ---- Test 3: include_selected = FALSE quadrants ---------------------------

test_that("getCohortATTsFinal(fused=TRUE, include_selected=FALSE) omits the `selected` column but still produces d_inv_treat_sel (API-lockdown)", {
	# This quadrant has no real-life caller — FETWFE always wants the
	# selection indicator. The test exists to lock in that the flags
	# are orthogonal: a future maintainer adding a fused-but-no-
	# selection-indicator caller shouldn't have to rewrite the unified
	# function. Without this test, the FALSE branch could be deleted as
	# "dead" during a future refactor.
	fx <- .build_unif_fixture()
	out <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = seq_len(fx$p),
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = TRUE,
		calc_ses = TRUE,
		include_selected = FALSE
	)

	expect_false("selected" %in% colnames(out$cohort_te_df))
	expect_equal(ncol(out$cohort_te_df), 6L)
	expect_true(!is.null(out$d_inv_treat_sel))
})

# ---- Test 4: psi_mat row-count contract across calling conventions --------

# Issue #116 Gap 4: getCohortATTsFinal() builds psi_mat with one row per
# entry of sel_treat_inds_shifted. Under the OLS (ETWFE / twfeCovs)
# convention the caller passes seq_len(num_treats), so nrow(psi_mat) ==
# num_treats; under the bridge (FETWFE / BETWFE) convention with PARTIAL
# selection it passes a strict subset, so nrow(psi_mat) < num_treats. The
# internal stopifnot()s only bound nrow(psi_mat) <= num_treats -- the
# partial-selection row count is otherwise unguarded, so a silent swap of
# the two conventions would go uncaught. This test pins both shapes.
test_that("getCohortATTsFinal psi_mat row count tracks the calling convention (OLS vs bridge)", {
	fx <- .build_unif_fixture()

	# --- OLS convention: full sel_treat_inds_shifted, sel_feat_inds = NULL.
	# nrow(psi_mat) must equal num_treats. (Cheap insurance: the fused = TRUE
	# block above already covers the full-set case; the bridge-convention
	# assertion below is the genuinely new coverage.)
	out_ols <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = NULL,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = seq_len(fx$num_treats),
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = FALSE
	)
	expect_equal(nrow(out_ols$psi_mat), fx$num_treats)
	expect_equal(ncol(out_ols$psi_mat), fx$G)

	# --- Bridge convention: strict-subset sel_treat_inds_shifted.
	# getGramInv() enforces sum(sel_feat_inds %in% treat_inds) ==
	# length(sel_treat_inds_shifted), so sel_feat_inds and
	# sel_treat_inds_shifted must be co-derived from the SAME selection (as
	# the real FETWFE caller in R/fetwfe_core.R is). Build sel_feat by
	# keeping every non-treatment feature plus only the treatment columns
	# named by the subset.
	sel_subset <- c(1L, 3L, 5L) # strict subset, length 3 < num_treats (5)
	non_treat_inds <- setdiff(seq_len(fx$p), fx$treat_inds)
	sel_feat <- sort(c(non_treat_inds, fx$treat_inds[sel_subset]))

	out_bridge <- fetwfe:::getCohortATTsFinal(
		X_final = fx$X_final,
		sel_feat_inds = sel_feat,
		treat_inds = fx$treat_inds,
		num_treats = fx$num_treats,
		first_inds = fx$first_inds,
		sel_treat_inds_shifted = sel_subset,
		c_names = fx$c_names,
		tes = fx$tes,
		sig_eps_sq = fx$sig_eps_sq,
		G = fx$G,
		N = fx$N,
		T = fx$T,
		fused = TRUE,
		calc_ses = TRUE,
		include_selected = TRUE
	)
	# Under partial bridge selection nrow(psi_mat) tracks the subset length,
	# strictly below num_treats -- the contract a convention swap would break.
	expect_equal(nrow(out_bridge$psi_mat), length(sel_subset))
	expect_lt(nrow(out_bridge$psi_mat), fx$num_treats)
	expect_equal(ncol(out_bridge$psi_mat), fx$G)
})
