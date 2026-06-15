# Multiplier-bootstrap simultaneous bands for simultaneousCIs() (#142, Phase 1).
#
# An alternative to the analytic `mvtnorm::qmvnorm()` critical value: perturb the
# per-unit influence-function (IF) matrix `F` (N_units x K) with random
# multipliers and read the sup-t critical value off the bootstrap distribution
# (Chernozhukov-Chetverikov-Kato 2013). Asymptotically equivalent to the analytic
# band under the paper's tight-Gaussianity, with three practical edges: it scales
# to large families (where `qmvnorm` strains), it is heteroskedasticity/cluster
# robust by construction (the empirical per-unit IF), and it is the path that
# generalizes to the high-dimensional regime (Phase 2, where the analytic Gram
# inverse does not exist).
#
# Phase 1 covers the regression-channel families (`cohort`,
# `all_post_treatment`, `custom`), where the cohort-probability variance term
# `Sigma_2` is zero, so `F = F_reg`. The `event_study` family (which needs the
# per-unit propensity IF `F_pi`) and the high-dimensional regime are Phase 2.

#' Draw multiplier-bootstrap weights
#'
#' @description `N x B` matrix of mean-zero, unit-variance multipliers, one column
#'   per bootstrap replicate. `"rademacher"` draws `+/-1` with probability 1/2;
#'   `"mammen"` draws the two-point distribution of Mammen (1993) (mean 0,
#'   variance 1, third moment 1).
#' @keywords internal
#' @noRd
.draw_multipliers <- function(N, B, type = c("rademacher", "mammen")) {
	type <- match.arg(type)
	if (type == "rademacher") {
		matrix(
			sample(c(-1, 1), size = N * B, replace = TRUE),
			nrow = N,
			ncol = B
		)
	} else {
		# Mammen two-point: values a1 < 0 < a2 with P(a1) = p1.
		phi <- (1 + sqrt(5)) / 2
		a1 <- 1 - phi # = -(sqrt(5)-1)/2
		a2 <- phi # = (sqrt(5)+1)/2
		p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
		draws <- ifelse(
			matrix(stats::runif(N * B), nrow = N, ncol = B) < p1,
			a1,
			a2
		)
		matrix(draws, nrow = N, ncol = B)
	}
}

#' Per-unit regression influence-function matrix `F_reg` (N_units x K)
#'
#' @description
#' For each effect `k` in the family, the per-unit regression IF is
#' `F_reg[i, k] = sum_{t} (x_{it}' v_k) * resid_{it}`, where `v_k` is the
#' debiasing direction for effect `k` and `resid` are the OLS-refit residuals on
#' the selected support (the same residuals the cluster-robust sandwich uses).
#' Vectorized as `F_reg = XEps %*% V`, with `XEps = rowsum(X_sel_c * resid, unit)`
#' (the per-unit aggregated score, `N x p_sel`) and `V` the `p_sel x K` matrix of
#' debiasing directions. This is `debiasedATT()`'s `(X %*% v) * resid` per-unit
#' score generalized to `K` columns, in `getGramInv()`'s centered / `n = NT`
#' convention. By construction `crossprod(F_reg) / (N*T)^2` equals the
#' cluster-robust `Sigma_1` divided by the finite-sample adjustment `N/(N-1)`.
#'
#' @param X_sel Numeric `NT x p_sel`; the SELECTED-support design (uncentered).
#' @param y Numeric; the (post-GLS) response, length `>= NT`.
#' @param N,T Integers; units and periods.
#' @param Psi_full Numeric `p_sel x K`; effect contrasts mapped to the selected
#'   support and zero-padded over nuisance columns (column `k` is effect `k`'s
#'   direction in selected coefficient space).
#' @return `F_reg`, an `N x K` matrix.
#' @keywords internal
#' @noRd
.build_regression_if <- function(X_sel, y, N, T, Psi_full) {
	n <- N * T
	X_sel_c <- scale(X_sel, center = TRUE, scale = FALSE)
	# OLS-refit residuals on the selected support (matches
	# .compute_cluster_robust_sandwich; the intercept is absorbed by centering).
	resid <- stats::lm.fit(cbind(1, X_sel), y[seq_len(n)])$residuals
	unit <- rep(seq_len(N), each = T)
	XEps <- rowsum(X_sel_c * resid, group = unit, reorder = FALSE) # N x p_sel
	Sig <- crossprod(X_sel_c) / n
	K <- ncol(Psi_full)
	V <- matrix(0, ncol(X_sel_c), K)
	for (k in seq_len(K)) {
		V[, k] <- solve(Sig, Psi_full[, k])
	}
	XEps %*% V
}

#' Studentized sup-t multiplier-bootstrap critical value
#'
#' @description
#' Given the per-unit IF matrix `F` (`N x K`), draws `B` multiplier-bootstrap
#' replicates of the studentized max-statistic
#' `T_b = max_k |sum_i xi_i F[i,k]| / sd_k` (`sd_k = sqrt(crossprod(F)[k,k])`) and
#' returns the `(1 - alpha)` quantile as the simultaneous critical value, along
#' with the per-effect empirical standard errors `se_k = sd_k / (N*T)` (the
#' cluster-robust SEs implied by `F`). Degenerate (zero-variance) effects are
#' excluded from the max (they would divide by ~0); their CIs collapse to the
#' point estimate.
#'
#' @return A list with `crit` (scalar), `ses` (length-K), and `nondeg`
#'   (logical length-K).
#' @keywords internal
#' @noRd
.simultaneous_bootstrap_crit <- function(
	F_mat,
	n,
	alpha,
	B,
	multiplier = c("rademacher", "mammen"),
	seed = NULL
) {
	multiplier <- match.arg(multiplier)
	N <- nrow(F_mat)
	K <- ncol(F_mat)
	n_obs_sq <- n^2 # (N*T)^2
	col_ss <- colSums(F_mat^2) # crossprod diagonal = sd_k^2
	var_tol <- .Machine$double.eps^0.5 * max(col_ss, 1)
	nondeg <- col_ss > var_tol
	sd_k <- sqrt(col_ss)
	# Per-effect cluster-robust SE: crossprod(F)/n^2 = cluster Sigma / cadjust, so
	# the matching SE restores the finite-sample factor cadjust = N/(N-1). (It
	# cancels in the studentized `crit` below, so it only scales band widths.)
	cadjust <- N / (N - 1)
	ses <- sqrt(cadjust * col_ss / n_obs_sq)

	draw_boot <- function() {
		Xi <- .draw_multipliers(N, B, multiplier) # N x B
		G <- crossprod(F_mat[, nondeg, drop = FALSE], Xi) # K_nd x B
		apply(abs(G) / sd_k[nondeg], 2, max) # length-B sup-t sample T_b
	}
	# With 0 or 1 non-degenerate effect the sup-t statistic reduces to the
	# pointwise |Z|/sd, whose 1-alpha quantile is exactly qnorm(1 - alpha/2); use
	# it directly rather than the Monte-Carlo estimate (which can dip below the
	# pointwise value and make a K=1 band anti-conservative) -- and skip the draw
	# entirely (drawing would also take max() over an empty set when every effect
	# is degenerate). `boot_max = NULL` signals this branch to the caller.
	if (sum(nondeg) <= 1L) {
		crit <- stats::qnorm(1 - alpha / 2)
		boot_max <- NULL
	} else {
		# Conv A: a non-NULL seed gives reproducible draws and is restored
		# afterward; a NULL seed draws from the ambient RNG (do NOT wrap -- that
		# would re-seed).
		boot_max <- if (is.null(seed)) {
			draw_boot()
		} else {
			.with_preserved_rng(seed, draw_boot())
		}
		crit <- stats::quantile(
			boot_max,
			probs = 1 - alpha,
			names = FALSE,
			type = 7
		)
	}

	list(crit = crit, ses = ses, nondeg = nondeg, boot_max = boot_max)
}

#' Bootstrap path for `simultaneousCIs()` (Phase 1: regression-channel families)
#'
#' @description
#' Called by `.simultaneous_cis_impl()` when `method = "bootstrap"`. Builds the
#' per-unit IF matrix from the already-resolved selected support + effect
#' contrasts, runs the multiplier bootstrap, and returns the same
#' `"simultaneous_cis"` S3 object the analytic path returns, plus `method` / `B`
#' / `seed` / `multiplier` fields. Phase 1 supports `cohort` /
#' `all_post_treatment` / `custom` (the cohort-probability variance term is zero,
#' so `F = F_reg`); `event_study` needs the per-unit propensity IF and is a
#' follow-up.
#' @keywords internal
#' @noRd
.simultaneous_cis_bootstrap <- function(
	family,
	alpha,
	B,
	seed,
	multiplier,
	X_final,
	y_final,
	N,
	T,
	treat_inds,
	sel_feat_inds,
	sel_treat_inds_shifted,
	Psi,
	K,
	estimates,
	effect_labels,
	pointwise_crit,
	bonferroni_crit
) {
	if (identical(family, "event_study")) {
		stop(
			"simultaneousCIs(): method = 'bootstrap' does not yet support ",
			"family = 'event_study' (its cohort-probability variance channel ",
			"needs a per-unit propensity influence function, planned as a ",
			"follow-up). Use method = 'analytic' for event_study, or method = ",
			"'bootstrap' with family = 'cohort', 'all_post_treatment', or ",
			"'custom'.",
			call. = FALSE
		)
	}
	n <- N * T

	# X_sel + zero-padded Psi_full over the selected support. FETWFE/BETWFE select
	# a subset (`sel_feat_inds`); etwfe/twfeCovs use the full design.
	if (any(!is.na(sel_feat_inds))) {
		X_sel <- X_final[, sel_feat_inds, drop = FALSE]
		Psi_full <- matrix(0, length(sel_feat_inds), K)
		treat_pos <- match(treat_inds[sel_treat_inds_shifted], sel_feat_inds)
		Psi_full[treat_pos, ] <- Psi
	} else {
		X_sel <- X_final
		Psi_full <- matrix(0, ncol(X_final), K)
		Psi_full[treat_inds[sel_treat_inds_shifted], ] <- Psi
	}

	F_mat <- .build_regression_if(X_sel, y_final, N, T, Psi_full)
	bc <- .simultaneous_bootstrap_crit(
		F_mat,
		n = n,
		alpha = alpha,
		B = B,
		multiplier = multiplier,
		seed = seed
	)
	crit <- bc$crit
	ses <- bc$ses

	# Single-step max-T adjusted p-values, the exact dual of the band. With > 1
	# non-degenerate effect: from the bootstrap sup-t distribution,
	# p_k = mean_b{ T_b >= |estimate_k| / se_k }. With a single non-degenerate
	# effect `boot_max` is NULL and the band uses the exact qnorm crit, so the
	# matching dual is the two-sided normal p-value 2 * pnorm(-|z|) (mirrors the
	# analytic K = 1 path; keeps "outside band iff adjusted p < alpha" exact).
	t_stat <- ifelse(bc$nondeg & ses > 0, abs(estimates) / ses, NA_real_)
	adjusted_p_values <- if (is.null(bc$boot_max)) {
		2 * stats::pnorm(-abs(t_stat))
	} else {
		vapply(
			t_stat,
			function(t) if (is.na(t)) NA_real_ else mean(bc$boot_max >= t),
			numeric(1)
		)
	}

	ci <- data.frame(
		effect = effect_labels,
		estimate = estimates,
		simultaneous_ci_low = estimates - crit * ses,
		simultaneous_ci_high = estimates + crit * ses,
		pointwise_ci_low = estimates - pointwise_crit * ses,
		pointwise_ci_high = estimates + pointwise_crit * ses,
		stringsAsFactors = FALSE
	)
	out <- list(
		ci = ci,
		adjusted_p_values = adjusted_p_values,
		critical_value = crit,
		pointwise_critical_value = pointwise_crit,
		bonferroni_critical_value = bonferroni_crit,
		family = family,
		alpha = alpha,
		K = K,
		method = "bootstrap",
		B = B,
		seed = seed,
		multiplier = multiplier
	)
	class(out) <- "simultaneous_cis"
	out
}
