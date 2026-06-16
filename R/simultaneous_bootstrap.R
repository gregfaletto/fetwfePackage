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

#' Per-unit regression influence-function matrix `F_reg` (N_units x K), fixed-p
#'
#' @description
#' For each effect `k`, the per-unit regression IF is
#' `F_reg[i, k] = sum_t (x_{it}' v_k) * resid_{it}`, vectorized as
#' `F_reg = XEps %*% V` with `XEps = rowsum(X_sel_c * resid, unit)` and
#' `v_k = solve(Sig, Psi_full[,k])`, `Sig = X_sel_c'X_sel_c / n`. `resid` are the
#' OLS-refit residuals on the (centered) selected support, so by construction
#' `crossprod(F_reg)/(NT)^2` equals the cluster-robust `Sigma_1` up to the
#' finite-sample factor `N/(N-1)`. This is `debiasedATT()`'s `(X v) * resid` score
#' generalized to `K` columns. Fixed-`p` (`p_sel < NT`) only; see
#' `.build_regression_if_highdim()` for the `p >= NT` regime.
#'
#' @param X_sel Numeric `NT x p_sel`; the SELECTED-support design (uncentered).
#' @param y Numeric; the (post-GLS) response, length `>= NT`.
#' @param N,T Integers; units and periods.
#' @param Psi_full Numeric `p_sel x K`; effect contrasts mapped to the selected
#'   support and zero-padded over nuisance columns.
#' @return `F_reg`, an `N x K` matrix, with attribute `highdim = FALSE`.
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
	F_mat <- XEps %*% V
	attr(F_mat, "highdim") <- FALSE
	F_mat
}

#' Per-unit regression IF, high-dimensional (p >= NT) regime (#142 Phase 2)
#'
#' @description
#' The full-design desparsified construction of `debiasedATT()` (#31), generalized
#' to `K` effects. With `p >= NT` the OLS refit overfits and the Gram `X'X/n` is
#' singular, so this uses the *full* (uncentered) design `X`, the regularized
#' **bridge** residuals `resid`, and a per-effect **nodewise (desparsified-lasso)**
#' direction `v_k = riesz_lasso(Sig, targets[,k], lambda_node_k)`. Then
#' `F[i,k] = sum_t (x_{it}' v_k) * resid_{it}` -- exactly `debiasedATT()`'s
#' high-dim per-unit score, one direction per effect. `targets[,k]` is effect
#' `k`'s direction in the full theta-space (`A' a_beta_k`). NOT post-selection:
#' uniformly valid under the desparsified-lasso conditions of paper Theorem
#' `debiased.highdim.thm` (experimental). `lambda_node = lambda_c * max(|target|)
#' * sqrt(log p / N)`.
#'
#' @param X Numeric `NT x p`; the FULL (uncentered) design.
#' @param resid Numeric length `NT`; bridge residuals
#'   `y - theta_hat[1] - X theta_hat[-1]`.
#' @param N,T Integers.
#' @param targets Numeric `p x K`; per-effect full theta-space directions.
#' @param lambda_c,riesz_max_iter,riesz_tol Nodewise-solver controls.
#' @return `F`, `N x K`, with attributes `highdim = TRUE` and `diagnostics`
#'   (per-effect feasibility / converged / lambda_node).
#' @keywords internal
#' @noRd
.build_regression_if_highdim <- function(
	X,
	resid,
	N,
	T,
	targets,
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	n <- N * T
	p <- ncol(X)
	Sig <- crossprod(X) / n # full, uncentered (debiasedATT convention)
	unit <- rep(seq_len(N), each = T)
	K <- ncol(targets)
	F_mat <- matrix(0, N, K)
	feasibility <- numeric(K)
	converged <- logical(K)
	lambda_node <- numeric(K)
	for (k in seq_len(K)) {
		a_k <- targets[, k]
		lambda_node[k] <- lambda_node_default(
			p = p,
			N = N,
			c = lambda_c,
			scale = max(abs(a_k))
		)
		v_k <- riesz_lasso(
			Sig,
			a_k,
			lambda_node[k],
			max_iter = riesz_max_iter,
			tol = riesz_tol
		)
		score_k <- as.numeric((X %*% v_k) * resid)
		F_mat[, k] <- rowsum(score_k, group = unit, reorder = FALSE)
		feasibility[k] <- attr(v_k, "feasibility")
		converged[k] <- attr(v_k, "converged")
	}
	attr(F_mat, "highdim") <- TRUE
	attr(F_mat, "diagnostics") <- list(
		feasibility = feasibility,
		converged = converged,
		lambda_node = lambda_node
	)
	F_mat
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
	# cancels in the studentized `crit` below, so it only scales band widths.) The
	# high-dimensional channel reuses the same cadjust for uniformity -- making its
	# bands marginally (N/(N-1)) wider than debiasedATT's regression SE, which omits
	# it; harmless, and -> 1 as N grows.
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
	lambda_c,
	riesz_max_iter,
	riesz_tol,
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
	bonferroni_crit,
	theta_hat_full = NULL,
	targets = NULL
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

	# Regime dispatch. `targets` (the full-design per-effect theta-space
	# directions) is non-NULL iff the caller detected the high-dimensional full
	# `p >= NT` regime: use the full-design desparsified construction (#31
	# generalized to K effects -- uniformly valid, NOT post-selection). Otherwise
	# the fixed-p selected-support construction (Phase 1).
	if (!is.null(targets)) {
		# Bridge residuals (debiasedATT convention) over the FULL design.
		resid_bridge <- as.numeric(
			y_final[seq_len(n)] -
				theta_hat_full[1] -
				X_final %*% theta_hat_full[-1]
		)
		F_mat <- .build_regression_if_highdim(
			X_final,
			resid_bridge,
			N,
			T,
			targets,
			lambda_c = lambda_c,
			riesz_max_iter = riesz_max_iter,
			riesz_tol = riesz_tol
		)
	} else {
		# Fixed-p: selected support + zero-padded Psi_full. FETWFE/BETWFE select a
		# subset (`sel_feat_inds`); etwfe/twfeCovs use the full design.
		if (any(!is.na(sel_feat_inds))) {
			X_sel <- X_final[, sel_feat_inds, drop = FALSE]
			Psi_full <- matrix(0, length(sel_feat_inds), K)
			treat_pos <- match(
				treat_inds[sel_treat_inds_shifted],
				sel_feat_inds
			)
			Psi_full[treat_pos, ] <- Psi
		} else {
			X_sel <- X_final
			Psi_full <- matrix(0, ncol(X_final), K)
			Psi_full[treat_inds[sel_treat_inds_shifted], ] <- Psi
		}
		F_mat <- .build_regression_if(X_sel, y_final, N, T, Psi_full)
	}

	# High-dimensional regime is EXPERIMENTAL: warn if any per-effect nodewise
	# direction did not meet its KKT feasibility constraint or converge (its band
	# may be unreliable; the returned diagnostics flag which effects).
	if (isTRUE(attr(F_mat, "highdim"))) {
		dg <- attr(F_mat, "diagnostics")
		bad <- (dg$feasibility > dg$lambda_node * (1 + riesz_tol)) |
			!dg$converged
		if (any(bad)) {
			warning(
				"simultaneousCIs(): the high-dimensional nodewise debiasing ",
				"direction for ",
				sum(bad),
				" of ",
				length(bad),
				" effect(s) did not meet its feasibility constraint / converge; ",
				"those simultaneous bands may be unreliable. Raise `riesz_max_iter` ",
				"and/or `lambda_c` and inspect the returned `feasibility` / ",
				"`converged` diagnostics. The p >= NT path is experimental.",
				call. = FALSE
			)
		}
	}
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
		multiplier = multiplier,
		regime = if (isTRUE(attr(F_mat, "highdim"))) {
			"high-dimensional"
		} else {
			"fixed-p"
		}
	)
	# High-dimensional fits append the per-effect nodewise-direction diagnostics
	# (feasibility = ||Sig v - a||_inf, the KKT certificate; convergence flags;
	# the penalties used).
	if (isTRUE(attr(F_mat, "highdim"))) {
		d <- attr(F_mat, "diagnostics")
		out$feasibility <- d$feasibility
		out$converged <- d$converged
		out$lambda_node <- d$lambda_node
	}
	class(out) <- "simultaneous_cis"
	out
}
