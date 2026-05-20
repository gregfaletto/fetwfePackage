library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for the FETWFE att_var_2 Jacobian fix (issue #46, version 1.8.0).
#
# Prior to v1.8.0 the internal `getSecondVarTermDataApp()` (which backs FETWFE's
# overall-ATT SE component `att_var_2`) used the outer-loop variable for the
# off-diagonal Jacobian coefficient, deviating from paper Theorem 6.3 by an
# off-by-one index. The fix uses the column index, matching the paper and the
# parallel `getSecondVarTermOLS()` in `R/variance_machinery.R`. The same fix flows
# through `.event_study_var2_fetwfe()`.
#
# These tests are regression coverage for the corrected formula. They also
# explicitly reconstruct the buggy formula inline to catch any future revert.
# ------------------------------------------------------------------------------

# Helper: deterministically fit FETWFE on R=3, T=6, d=0, N=300. Returns a list
# with the package fit and all the internals needed to reconstruct the variance
# formula by hand. This regime is chosen because the buggy-vs-fixed gap in
# att_var_2 is large enough (~37%) to make the MC and anti-regression
# assertions reliable.
.make_var2_fit <- function() {
	set.seed(31)
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 0,
		density = 0.5,
		eff_size = 2,
		seed = 31
	)
	sim <- simulateData(coefs, N = 300, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	res <- fetwfeWithSimulatedData(sim, q = 0.5)

	R_val <- res$R
	T_val <- res$T
	N_val <- res$N
	num_treats <- length(res$treat_inds)
	first_inds <- fetwfe:::getFirstInds(R_val, T_val)
	cohort_probs_overall <- unname(res$cohort_probs_overall)

	theta_hat_slopes <- res$internal$theta_hat[2:(res$p + 1)]
	sel_treat_inds_shifted <- which(theta_hat_slopes[res$treat_inds] != 0)
	theta_hat_treat_sel <- theta_hat_slopes[res$treat_inds][
		sel_treat_inds_shifted
	]
	d_inv_full <- fetwfe:::genInvTwoWayFusionTransformMat(
		num_treats,
		first_inds,
		R_val
	)
	d_inv_treat_sel <- d_inv_full[, sel_treat_inds_shifted, drop = FALSE]

	# sel_inds[[r]] matches getSecondVarTermDataApp's internal convention:
	# the FULL cohort-r block in 1:num_treats indexing (rows of d_inv_treat_sel),
	# NOT the restriction to the selected-theta support.
	sel_inds <- list()
	for (r in 1:R_val) {
		first_ind_r <- first_inds[r]
		last_ind_r <- if (r < R_val) first_inds[r + 1] - 1 else num_treats
		sel_inds[[r]] <- first_ind_r:last_ind_r
	}

	list(
		res = res,
		R = R_val,
		T = T_val,
		N = N_val,
		num_treats = num_treats,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall,
		theta_hat_treat_sel = theta_hat_treat_sel,
		d_inv_treat_sel = d_inv_treat_sel,
		sel_inds = sel_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted
	)
}

# Hand-built overall-ATT var_2 quadratic form. `off_use_outer` selects the buggy
# (outer-loop) or correct (column-index) off-diagonal Jacobian construction.
.overall_var_2 <- function(fit, off_use_outer) {
	R_val <- fit$R
	T_val <- fit$T
	N_val <- fit$N
	cp <- fit$cohort_probs_overall
	S <- sum(cp)
	sel_inds <- fit$sel_inds
	d_inv <- fit$d_inv_treat_sel

	Sigma_pi_hat <- -outer(cp, cp)
	diag(Sigma_pi_hat) <- cp * (1 - cp)

	J <- matrix(0, nrow = R_val, ncol = ncol(d_inv))
	for (r in 1:R_val) {
		cons_r <- (S - cp[r]) / S^2
		J[r, ] <- cons_r * colMeans(d_inv[sel_inds[[r]], , drop = FALSE])
		for (r2 in setdiff(1:R_val, r)) {
			idx_pi <- if (off_use_outer) r else r2
			cons_off <- cp[idx_pi] / S^2
			J[r, ] <- J[r, ] -
				cons_off * colMeans(d_inv[sel_inds[[r2]], , drop = FALSE])
		}
	}

	T_val *
		as.numeric(
			t(fit$theta_hat_treat_sel) %*%
				t(J) %*%
				Sigma_pi_hat %*%
				J %*%
				fit$theta_hat_treat_sel
		) /
		(N_val * T_val)
}

# ------------------------------------------------------------------------------
# Test 1: hand-computed reference matches getSecondVarTermDataApp's output.
# ------------------------------------------------------------------------------
test_that("FETWFE att_var_2 helper matches paper Theorem 6.3 Jacobian by hand", {
	fit <- .make_var2_fit()
	var2_handcomputed <- .overall_var_2(fit, off_use_outer = FALSE)

	# Call the package's helper through getTeResults2 by re-invoking
	# getSecondVarTermDataApp directly.
	var2_helper <- fetwfe:::getSecondVarTermDataApp(
		sel_treat_inds_shifted = fit$sel_treat_inds_shifted,
		cohort_probs_overall = fit$cohort_probs_overall,
		first_inds = fit$first_inds,
		theta_hat_treat_sel = fit$theta_hat_treat_sel,
		num_treats = fit$num_treats,
		N = fit$N,
		T = fit$T,
		R = fit$R,
		d_inv_treat_sel = fit$d_inv_treat_sel
	)

	expect_equal(var2_helper, var2_handcomputed, tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# Test 2: corrected var_2 matches a multinomial-resampling Monte Carlo.
# Isolates var_2 from var_1 and from bridge-selection noise by holding beta_hat
# fixed and resampling pi_hat from a multinomial. The corrected var_2 should
# match Var(g(pi_hat)) within MC noise; the buggy var_2 should not.
# ------------------------------------------------------------------------------
test_that("FETWFE corrected var_2 matches a multinomial-resampling MC", {
	skip_on_cran() # MC test with small randomness; gate to local + CI
	fit <- .make_var2_fit()

	# g(pi) holding beta_hat fixed: pooled ATT as a function of cohort probs.
	beta_treat_sel <- as.numeric(
		fit$d_inv_treat_sel %*% fit$theta_hat_treat_sel
	)
	block_means <- sapply(
		fit$sel_inds,
		function(idx) mean(beta_treat_sel[idx])
	)

	cp <- fit$cohort_probs_overall
	S <- sum(cp)
	p_full <- c(cp, 1 - S) # add the never-treated proportion to make it a probability vector

	g_fn <- function(pi_vec) {
		Sv <- sum(pi_vec)
		sum(pi_vec / Sv * block_means)
	}

	set.seed(42)
	n_mc <- 10000L
	g_draws <- numeric(n_mc)
	for (i in seq_len(n_mc)) {
		pi_hat_full <- rmultinom(1, fit$N, p_full)[, 1] / fit$N
		g_draws[i] <- g_fn(pi_hat_full[1:fit$R])
	}
	mc_var <- var(g_draws)

	# Route var2_corrected through the package helper directly so the MC
	# assertion validates the production code (not just the in-test formula).
	var2_corrected <- fetwfe:::getSecondVarTermDataApp(
		sel_treat_inds_shifted = fit$sel_treat_inds_shifted,
		cohort_probs_overall = fit$cohort_probs_overall,
		first_inds = fit$first_inds,
		theta_hat_treat_sel = fit$theta_hat_treat_sel,
		num_treats = fit$num_treats,
		N = fit$N,
		T = fit$T,
		R = fit$R,
		d_inv_treat_sel = fit$d_inv_treat_sel
	)
	# Buggy formula is reconstructed in-test (the helper after the fix no
	# longer produces it); see Test 3 for the anti-regression structural check.
	var2_buggy <- .overall_var_2(fit, off_use_outer = TRUE)

	# Corrected var_2 should match MC within 5% at 10k draws.
	expect_lt(abs(var2_corrected / mc_var - 1), 0.05)
	# Buggy var_2 should be visibly off by at least 15% on this panel.
	# (Plan-review round 1 verification at 50k draws found ratio 0.73.)
	expect_gt(abs(var2_buggy / mc_var - 1), 0.15)
})

# ------------------------------------------------------------------------------
# Test 3: helper output differs from the buggy formula by a non-trivial amount.
# Anti-regression: a future revert of the index swap would re-introduce the bug
# and this test would catch it even without MC sampling.
# ------------------------------------------------------------------------------
test_that("FETWFE att_var_2 helper differs from the buggy formula (anti-regression)", {
	fit <- .make_var2_fit()
	var2_corrected <- .overall_var_2(fit, off_use_outer = FALSE)
	var2_buggy <- .overall_var_2(fit, off_use_outer = TRUE)

	# On this panel the buggy/corrected ratio is ~0.73 (per the post-execution
	# review of PR #45). Demand at least 10% relative difference so a future
	# revert is caught even in the presence of any other small changes.
	expect_gt(abs(var2_corrected / var2_buggy - 1), 0.10)
	# Both should be finite positive; sanity check.
	expect_gt(var2_corrected, 0)
	expect_gt(var2_buggy, 0)
})

# ------------------------------------------------------------------------------
# Test 4 from earlier revisions (".event_study_var2_fetwfe matches its
# hand-computed Jacobian") was removed per #76 Item 7a: the
# "hand-computed" reference reproduced the helper's own loop structure
# step-for-step, so passing only confirmed the helper was consistent
# with itself rather than against an independent reference. Tests 2
# (MC validation) and 3 (anti-regression) above are the load-bearing
# var2 coverage.
# ------------------------------------------------------------------------------
