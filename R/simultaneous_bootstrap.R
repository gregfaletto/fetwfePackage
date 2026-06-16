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
# Phase 1 covered the regression-channel families (`cohort`,
# `all_post_treatment`, `custom`), where the cohort-probability variance term
# `Sigma_2` is zero, so `F = F_reg`. Phase 2 added the `event_study` family (the
# per-unit propensity IF `F_pi`) and the high-dimensional (`p >= NT`) regime --
# including their combination (high-dim `event_study` uses BOTH channels), and
# the high-dim band center is debiased (Theorem 6.6) rather than the raw bridge.

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

#' Per-unit cohort-probability (propensity) influence-function matrix `F_pi`
#'
#' @description
#' The `event_study` family is the one family whose simultaneous-band variance
#' carries a non-zero cohort-probability (propensity) channel `Sigma_2`: each
#' event-time effect pools across the cohorts treated by that event time, with
#' weights equal to the estimated cohort-membership probabilities
#' `pi_hat_g = N_g / N`, so estimating those probabilities contributes variance.
#' This builds the per-unit influence function of that channel, the delta-method
#' image of the cohort-sample-proportion IF `xi_i = e_{W_i} - psi*` (paper
#' `\eqref{psi.if.assum}`, line 2072; `e_{W_i}` is the one-hot of unit `i`'s
#' cohort, `0` for never-treated):
#'
#'   `F_pi[i, k] = T * a_k' (e_{W_i} - pi_hat)`,  `a_k = J_k theta_sel`  (length `G`).
#'
#' Per-effect this reproduces the analytic `V_2 = v_psi = M' grad-f' Sigma_M
#' grad-f M` (line 2925), so `colSums(F_pi^2) / (N*T)^2 == diag(Sigma_2)` exactly
#' (because `sum_i (e_{W_i} - pi_hat)(e_{W_i} - pi_hat)' = N * Sigma_pi`). The
#' leading `c = T` is PINNED by matching that anchor against
#' `.assemble_joint_cov_var2`'s `diag(Sigma_2)_k = T * a_k' Sigma_pi a_k / (N*T)
#' = a_k' Sigma_pi a_k / N` -- no factor-of-N/T/NT freedom.
#'
#' This channel is perturbed with its OWN independent multiplier stream (`eta`),
#' separate from the regression channel's (`xi`): the paper's per-unit IF is
#' `f_i = f_reg,i + f_pi,i` with `Var(f_i) = V_1 + V_2` and ZERO cross-term
#' (proof of Thm C.1, "Step 5 -- verify the cross-covariance vanishes", line
#' 3641: `E[IF^{(beta)}_{(it)} xi_i'] = 0` since `w_{it}` and `xi_i` are both
#' `(W_i, X_i)`-measurable while the GLS residual is conditionally mean-zero).
#' A SHARED multiplier would inject the empirical `sum_i F_reg[i,k] F_pi[i,k]`
#' cross-term (absent from the analytic `Sigma`), mis-studentizing the sup-t.
#'
#' Per-unit ordering is irrelevant (the two channels use independent multipliers,
#' so there is no cross-channel pairing to preserve), so `F_pi` is assembled from
#' cohort COUNTS alone: with `A = [a_1 ... a_K]` (`G x K`) and the constant row
#' `c0 = pi_hat' A` (length `K`), stack `n_g = round(N * pi_hat_g)` copies of the
#' centered one-hot `T * (A[g, ] - c0)` for `g = 1..G`, then `n_never = N - sum
#' n_g` copies of `-T * c0` (the never-treated row, `e = 0`). `Sigma_2 != 0` only
#' for `event_study`; the other families have zero Jacobians, hence `F_pi == 0`.
#'
#' @param J_list Length-`K` list of per-effect cohort-weight Jacobians from
#'   `.build_j_list_for_family(family = "event_study", ...)`; each is `G x p_sel`.
#' @param theta_sel Numeric length `p_sel`; the selected-support treatment-effect
#'   parameters (theta-space for FETWFE, beta-space for the OLS family / BETWFE).
#' @param cohort_probs_overall Numeric; cohort-membership probabilities
#'   `pi_hat_g = N_g / N` (length `>= G`; `[1:G]` are the treated cohorts, the
#'   residual mass is never-treated).
#' @param G Integer; number of treated cohorts.
#' @param N,T Integers; units and periods. (Period arg named `T` to match the
#'   sibling variance helpers `.assemble_joint_cov_var2()` / `.build_jacobian()`.)
#' @return `F_pi`, an `N x K` numeric matrix (`K = length(J_list)`).
#' @keywords internal
#' @noRd
.build_propensity_if <- function(
	J_list,
	theta_sel,
	cohort_probs_overall,
	G,
	N,
	T
) {
	K <- length(J_list)
	pi_hat <- cohort_probs_overall[seq_len(G)]
	# A = [a_1 ... a_K], forced to G x K via matrix(., nrow = G): the bare
	# vapply(..., numeric(G)) collapses to a length-K vector when G == 1, which
	# would break A[g, ] indexing for a single treated cohort.
	A <- matrix(
		vapply(
			seq_len(K),
			function(k) as.numeric(J_list[[k]] %*% theta_sel),
			numeric(G)
		),
		nrow = G,
		ncol = K
	)
	# Centered one-hot e_{W_i} - pi_hat (e = 0 for never-treated). Constant
	# subtracted row c0 = pi_hat' A (length K); cohort-g rows carry A[g, ] - c0
	# (one-hot e_g minus pi_hat), never-treated rows carry -c0.
	c0 <- as.numeric(crossprod(pi_hat, A))
	n_g <- round(N * pi_hat)
	n_never <- N - sum(n_g)
	if (n_never < 0L || sum(n_g) + n_never != N) {
		stop(
			"simultaneousCIs(): could not recover integer cohort counts from ",
			"`cohort_probs_overall` for the event_study propensity influence ",
			"function.",
			call. = FALSE
		)
	}
	blocks <- vector("list", G + 1L)
	for (g in seq_len(G)) {
		blocks[[g]] <- matrix(
			rep(T * (A[g, ] - c0), each = n_g[g]),
			nrow = n_g[g],
			ncol = K
		)
	}
	blocks[[G + 1L]] <- matrix(
		rep(-T * c0, each = n_never),
		nrow = n_never,
		ncol = K
	)
	do.call(rbind, blocks)
}

#' Studentized sup-t multiplier-bootstrap critical value
#'
#' @description
#' Given the per-unit regression IF matrix `F` (`N x K`) and, optionally, the
#' per-unit propensity IF matrix `F_pi` (`N x K`, `event_study` only), draws `B`
#' multiplier-bootstrap replicates of the studentized max-statistic
#' `T_b = max_k |sum_i xi_i F[i,k] + sum_i eta_i F_pi[i,k]| / sd_k`,
#' `sd_k = sqrt(css_reg_k + css_pi_k)` (`css_* = colSums(F_*^2)`), and returns the
#' `(1 - alpha)` quantile as the simultaneous critical value, along with the
#' per-effect standard errors. The two channels use INDEPENDENT multiplier
#' streams (`xi` for regression, `eta` for propensity), reproducing the paper's
#' zero-cross-term `Var(f_i) = V_1 + V_2` (Thm C.1 Step 5). The reported SE
#' applies the finite-sample factor `cadjust = N/(N-1)` to the regression channel
#' ONLY (the cluster sandwich bakes `cadjust` into `Sigma_1`; the multinomial
#' `Sigma_2` carries none): `se_k = sqrt((cadjust*css_reg_k + css_pi_k) / n^2)`
#' (`cadjust` cancels in the studentized `crit`). Degenerate (zero combined
#' variance) effects are excluded from the max; their CIs collapse to the point
#' estimate.
#'
#' When `F_pi_mat` is `NULL` (all non-event_study families, both regimes) the
#' `eta` stream is NOT drawn and the draw is byte-identical to the single-channel
#' Phase-1/2 behavior.
#'
#' @return A list with `crit` (scalar), `ses` (length-K), `nondeg`
#'   (logical length-K), and `boot_max` (length-B, or `NULL` in the degenerate
#'   single-effect branch).
#' @keywords internal
#' @noRd
.simultaneous_bootstrap_crit <- function(
	F_mat,
	F_pi_mat = NULL,
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
	css_reg <- colSums(F_mat^2) # regression-channel crossprod diagonal
	css_pi <- if (is.null(F_pi_mat)) {
		rep(0, K) # no propensity channel (non-event_study families)
	} else {
		colSums(F_pi_mat^2)
	}
	col_ss <- css_reg + css_pi # combined sd_k^2 = diag(Sigma_1) + diag(Sigma_2)
	# Degeneracy uses the COMBINED column sums (relative to max(col_ss, 1)), so a
	# propensity-only-variance effect is not dropped.
	var_tol <- .Machine$double.eps^0.5 * max(col_ss, 1)
	nondeg <- col_ss > var_tol
	sd_k <- sqrt(col_ss)
	# Per-effect SE. The regression channel: crossprod(F)/n^2 = cluster Sigma_1 /
	# cadjust, so the matching SE restores the finite-sample factor
	# cadjust = N/(N-1). The propensity channel's multinomial Sigma_2 carries NO
	# cadjust (the anchor colSums(F_pi^2)/n^2 == diag(Sigma_2) is exact), so
	# cadjust is applied to the regression channel ONLY -- the unique scaling that
	# makes the event_study bootstrap SE equal the analytic event_study SE to
	# machine precision. cadjust is deliberately NOT applied to the studentization
	# (`col_ss` above): the sup-t `crit` below calibrates on the empirical per-unit
	# IF correlation cov2cor(Sigma_1 + Sigma_2). For the Sigma_2 = 0 families that
	# equals the analytic correlation (a uniform scale cancels in cov2cor); for
	# event_study (Sigma_2 != 0) it differs from the analytic qmvnorm's
	# cov2cor(cadjust*Sigma_1 + Sigma_2) by an O(1/N) reweighting -- within MC
	# error, marginally conservative (see #302 for the optional exact-match draw
	# scaling). The high-dimensional channel reuses the same cadjust for
	# uniformity -- making its bands marginally (N/(N-1)) wider than debiasedATT's
	# regression SE, which omits it; harmless, and -> 1 as N grows.
	cadjust <- N / (N - 1)
	ses <- sqrt((cadjust * css_reg + css_pi) / n_obs_sq)

	draw_boot <- function() {
		Xi <- .draw_multipliers(N, B, multiplier) # N x B
		G <- crossprod(F_mat[, nondeg, drop = FALSE], Xi) # K_nd x B
		# Independent propensity-channel multiplier stream `eta`, drawn ONLY when
		# F_pi is present and placed AFTER the `xi` draw so the single-channel
		# (F_pi = NULL) RNG stream is byte-identical to Phase 1/2.
		if (!is.null(F_pi_mat)) {
			Eta <- .draw_multipliers(N, B, multiplier) # N x B, independent of Xi
			G <- G + crossprod(F_pi_mat[, nondeg, drop = FALSE], Eta)
		}
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

#' Bootstrap path for `simultaneousCIs()`
#'
#' @description
#' Called by `.simultaneous_cis_impl()` when `method = "bootstrap"`. Builds the
#' per-unit IF matrix from the already-resolved selected support + effect
#' contrasts, runs the multiplier bootstrap, and returns the same
#' `"simultaneous_cis"` S3 object the analytic path returns, plus `method` / `B`
#' / `seed` / `multiplier` fields. The regression-channel families (`cohort` /
#' `all_post_treatment` / `custom`) have zero cohort-probability variance, so
#' `F = F_reg` and `F_pi = NULL`. The `event_study` family additionally carries
#' the per-unit propensity IF `F_pi` (`.build_propensity_if()`), perturbed with
#' its own independent multiplier stream (the two-channel `Sigma = Sigma_1 +
#' Sigma_2` bootstrap); it needs `J_list` / `theta_sel` / `cohort_probs_overall`
#' / `G`. The propensity channel composes regime-agnostically, so high-dim
#' `event_study` (the desparsified `targets` path) carries `F_pi` too. In the
#' high-dim regime the band center is the debiased estimate
#' (`estimates + colSums(F_mat)/(N*T)`, the Theorem 6.6 realization matching
#' `debiasedATT()`); fixed-p centers on the (unbiased) bridge estimate unchanged.
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
	targets = NULL,
	J_list = NULL,
	theta_sel = NULL,
	cohort_probs_overall = NULL,
	G = NULL
) {
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
		# Debias the band center (Theorem 6.6). In the high-dim regime the bridge
		# (post-selection) estimate is biased; the per-effect desparsified
		# correction `colSums(F_mat)/(N*T)` is exactly `debiasedATT()`'s
		# `mean(score)` (the per-effect generalization of `att_db = sum(a*theta) +
		# mean(score)`), so `estimates_used` here equals `debiasedATT()$att` for the
		# matching overall-ATT contrast. Centering the band on the bridge estimate
		# instead would carry the post-selection bias (variance -> 0 with N, bias
		# does not) and undercover.
		estimates_used <- estimates + colSums(F_mat) / (N * T)
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
		# Fixed-p: the selected estimate is unbiased (selection consistency), so the
		# band centers on the bridge estimate unchanged (byte-identical to Phase 1/2).
		estimates_used <- estimates
	}

	# Propensity (cohort-probability) channel `F_pi`: non-zero only for the
	# event_study family, whose pooled event-time effects weight cohorts by the
	# estimated probabilities `pi_hat_g = N_g/N`. Built for event_study in BOTH
	# regimes (it depends only on cohort counts + `theta_sel`, so it composes
	# regime-agnostically with the high-dim desparsified regression channel) from
	# the same `J_list` / `Sigma_pi` the analytic `Sigma_2` uses, so the bootstrap
	# and analytic `Sigma_2` agree to machine precision. NULL for the other
	# families (`Sigma_2 = 0`), keeping their bootstrap draw byte-identical.
	F_pi_mat <- if (identical(family, "event_study")) {
		.build_propensity_if(
			J_list = J_list,
			theta_sel = theta_sel,
			cohort_probs_overall = cohort_probs_overall,
			G = G,
			N = N,
			T = T
		)
	} else {
		NULL
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
		F_pi_mat = F_pi_mat,
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
	t_stat <- ifelse(bc$nondeg & ses > 0, abs(estimates_used) / ses, NA_real_)
	adjusted_p_values <- if (is.null(bc$boot_max)) {
		2 * stats::pnorm(-abs(t_stat))
	} else {
		vapply(
			t_stat,
			function(t) if (is.na(t)) NA_real_ else mean(bc$boot_max >= t),
			numeric(1)
		)
	}

	# Band center `estimates_used`: the bridge `estimates` in fixed-p, the debiased
	# `estimates + colSums(F_mat)/(N*T)` in the high-dim regime (see the regime
	# dispatch above). Both the band endpoints AND the adjusted-p `t_stat` use it,
	# so "outside band iff adjusted p < alpha" stays exact.
	ci <- data.frame(
		effect = effect_labels,
		estimate = estimates_used,
		simultaneous_ci_low = estimates_used - crit * ses,
		simultaneous_ci_high = estimates_used + crit * ses,
		pointwise_ci_low = estimates_used - pointwise_crit * ses,
		pointwise_ci_high = estimates_used + pointwise_crit * ses,
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
