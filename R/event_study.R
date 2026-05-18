# Event-study aggregation and plotting for fetwfe / etwfe / betwfe.
# See vignettes/etwfe_betwfe_vignette.Rmd for usage; see
# `.plans/feat-event-study-plot/plan.md` for the design rationale.

# Bare column names referenced inside ggplot2::aes() in .plot_event_study()
# are intentional (NSE); declare them here so R CMD check doesn't emit a
# "no visible binding for global variable" NOTE.
utils::globalVariables(c("event_time", "estimate", "ci_low", "ci_high"))

#' Compute pooled event-time treatment-effect estimates
#'
#' @description
#' For a fitted object from `fetwfe()`, `etwfe()`, or `betwfe()`, computes the
#' pooled event-time treatment-effect estimates `tau_E(e)`, defined as
#' cohort-weighted averages of the cell-level treatment-effect estimates at
#' each post-treatment event time `e = t - r` (where `t` is calendar time and
#' `r` is the cohort's first-treated calendar time). Weights are
#' sample-cohort-size weights (matching `did::aggte(type = "dynamic")`
#' convention).
#'
#' Standard errors combine two terms, mirroring the package's existing
#' overall-ATT SE machinery: `var_1(e)` from regression-coefficient noise
#' (computed via the same `gram_inv` machinery the package uses for cohort
#' SEs, or the cluster-robust sandwich under `se_type = "cluster"`), and
#' `var_2(e)` from cohort-probability noise (analog of the existing
#' `getSecondVarTermOLS` / `getSecondVarTermDataApp` machinery, with the
#' multinomial Jacobian restricted to cohorts valid at event time `e`).
#' Combined as `sqrt(var_1 + var_2)` when `indep_counts` was supplied to the
#' fit (asymptotically exact), else the conservative Cauchy-Schwarz bound
#' `sqrt(var_1 + var_2 + 2 sqrt(var_1 * var_2))`.
#'
#' @param x A fitted object of class `"fetwfe"`, `"etwfe"`, or `"betwfe"`.
#' @param alpha (Optional) Significance level for confidence intervals.
#'   Defaults to `x$alpha` (the alpha used at fit time).
#' @return A data frame with class `c("fetwfe_event_study", "data.frame")` and
#'   columns:
#'   \describe{
#'     \item{event_time}{Integer; event time `e = t - r`, ranging from 0
#'       to `T - 2`.}
#'     \item{n_cohorts}{Integer; number of cohorts contributing to the
#'       pooled estimate at event time `e`.}
#'     \item{estimate}{Numeric; the pooled event-time ATT estimate.}
#'     \item{se}{Numeric; combined standard error.}
#'     \item{ci_low}{Numeric; lower bound of the (1 - alpha) Wald CI.}
#'     \item{ci_high}{Numeric; upper bound of the (1 - alpha) Wald CI.}
#'     \item{p_value}{Numeric; two-sided Wald p-value
#'       (`2 * pnorm(-|estimate / se|)`), `NA` when `se` is `0` or `NA`.}
#'   }
#'   Only post-treatment event times (`e >= 0`) are included; pre-treatment
#'   placebo periods would require an extended regression specification and
#'   are out of scope for this initial release.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- fetwfeWithSimulatedData(dat)
#'   event_study(res)
#' }
#' @export
event_study <- function(x, alpha = NULL) {
	if (inherits(x, "fetwfe")) {
		return(.event_study_fetwfe(x, alpha))
	}
	if (inherits(x, "etwfe") || inherits(x, "betwfe")) {
		return(.event_study_etwfe_betwfe(x, alpha))
	}
	stop(
		"event_study() requires an object of class 'fetwfe', 'etwfe', or 'betwfe'."
	)
}

#' Event-study aggregation for ETWFE / BETWFE
#' @keywords internal
#' @noRd
.event_study_etwfe_betwfe <- function(x, alpha = NULL) {
	# Method-entry precondition (#86). Validates the fitted object and
	# derives `has_valid_ses` — the contract gate that fixes #73
	# (event_study reporting finite SEs when the fit's att_se is NA
	# for q >= 1).
	contract <- .check_for_event_study(x)

	if (is.null(alpha)) {
		alpha <- x$alpha
	}
	z <- stats::qnorm(1 - alpha / 2)

	if (is.null(x$cohort_probs_overall)) {
		stop(
			"event_study(): cohort_probs_overall missing from object. Re-fit with etwfe() / betwfe() at version 1.7.0 or later."
		)
	}

	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	N <- x$N
	T <- x$T
	R <- x$R
	sig_eps_sq <- x$sig_eps_sq
	cohort_probs_overall <- x$cohort_probs_overall
	X_final <- x$X_final
	y_final <- x$y_final
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	first_inds <- getFirstInds(R = R, T = T)
	tes <- beta_hat[treat_inds]

	# Determine the selected support
	if (inherits(x, "etwfe")) {
		sel_feat_inds <- NA
		sel_treat_inds_shifted <- seq_len(num_treats)
	} else {
		# betwfe: selection is in beta-space
		sel_feat_inds <- which(beta_hat != 0)
		sel_treat_inds_shifted <- which(beta_hat[treat_inds] != 0)
	}

	calc_ses <- contract$has_valid_ses
	gram_inv <- NULL
	if (calc_ses && length(sel_treat_inds_shifted) > 0) {
		res_gram <- getGramInv(
			N = N,
			T = T,
			X_final = X_final,
			sel_feat_inds = sel_feat_inds,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			calc_ses = TRUE
		)
		gram_inv <- res_gram$gram_inv
		calc_ses <- res_gram$calc_ses
	} else {
		calc_ses <- FALSE
	}

	# Cluster-robust sandwich (recomputed from existing slots)
	sandwich_full <- NULL
	treat_block_mask <- NULL
	if (identical(se_type, "cluster") && calc_ses) {
		if (any(!is.na(sel_feat_inds))) {
			X_S <- X_final[, sel_feat_inds, drop = FALSE]
			treat_block_mask <- sel_feat_inds %in% treat_inds
		} else {
			X_S <- X_final
			treat_block_mask <- logical(ncol(X_S))
			treat_block_mask[treat_inds] <- TRUE
		}
		y_ <- y_final[seq_len(N * T)]
		ols_fit <- stats::lm.fit(cbind(1, X_S), y_)
		sandwich_full <- .compute_cluster_robust_sandwich(
			X_S = X_S,
			residuals = ols_fit$residuals,
			N = N,
			T = T
		)
	}

	max_event <- T - 2L
	event_times <- 0:max_event
	estimates <- numeric(length(event_times))
	ses <- numeric(length(event_times))
	n_cohorts <- integer(length(event_times))

	for (k in seq_along(event_times)) {
		e <- event_times[k]
		V_e <- which(seq_len(R) <= T - 1L - e)
		n_cohorts[k] <- length(V_e)
		if (length(V_e) == 0L) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}

		# Cohort-size weights, renormalized over V_e (sample-cohort-size convention)
		probs_Ve <- cohort_probs_overall[V_e]
		S_V <- sum(probs_Ve)
		if (S_V <= 0) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}
		weights_Ve <- probs_Ve / S_V

		# psi_e_tes (length num_treats): weight at idx(r, e) for r in V_e, 0 elsewhere
		psi_e_tes <- numeric(num_treats)
		for (j in seq_along(V_e)) {
			psi_e_tes[first_inds[V_e[j]] + e] <- weights_Ve[j]
		}

		# Point estimate
		estimates[k] <- as.numeric(crossprod(psi_e_tes, tes))

		if (!calc_ses) {
			ses[k] <- NA_real_
			next
		}

		# var_1(e) on the selected support
		if (identical(se_type, "cluster")) {
			psi_full <- numeric(ncol(sandwich_full))
			psi_e_sel <- psi_e_tes[sel_treat_inds_shifted]
			psi_full[treat_block_mask] <- psi_e_sel
			var_1_e <- as.numeric(
				t(psi_full) %*% sandwich_full %*% psi_full
			)
		} else {
			psi_e_sel <- psi_e_tes[sel_treat_inds_shifted]
			var_1_e <- sig_eps_sq *
				as.numeric(t(psi_e_sel) %*% gram_inv %*% psi_e_sel) /
				(N * T)
		}

		# var_2(e) via existing getSecondVarTermOLS with single-cell psi_mat and
		# masked cohort_probs_overall (zeros outside V_e).
		psi_e_mat <- matrix(0, nrow = num_treats, ncol = R)
		for (r in V_e) {
			psi_e_mat[first_inds[r] + e, r] <- 1
		}
		masked_probs <- numeric(R)
		masked_probs[V_e] <- cohort_probs_overall[V_e]
		if (length(V_e) == 1L) {
			# Single-cohort: var_2 vanishes by construction (weight = 1)
			var_2_e <- 0
		} else {
			var_2_e <- getSecondVarTermOLS(
				psi_mat = psi_e_mat,
				tes = tes,
				cohort_probs_overall = masked_probs,
				first_inds = first_inds,
				num_treats = num_treats,
				N = N,
				T = T,
				R = R
			)
		}

		# Combine
		if (is_indep) {
			ses[k] <- sqrt(var_1_e + var_2_e)
		} else {
			ses[k] <- sqrt(var_1_e + var_2_e + 2 * sqrt(var_1_e * var_2_e))
		}
	}

	.assemble_event_study_df(event_times, n_cohorts, estimates, ses, z)
}

#' Event-study aggregation for FETWFE
#' @keywords internal
#' @noRd
.event_study_fetwfe <- function(x, alpha = NULL) {
	# Method-entry precondition (#86). Validates the fitted object and
	# derives `has_valid_ses` — the contract gate that fixes #73
	# (event_study reporting finite SEs when the fit's att_se is NA
	# for q >= 1).
	contract <- .check_for_event_study(x)

	if (is.null(alpha)) {
		alpha <- x$alpha
	}
	z <- stats::qnorm(1 - alpha / 2)

	if (is.null(x$internal) || is.null(x$internal$theta_hat)) {
		stop(
			"event_study(): theta_hat missing from x$internal. Re-fit with fetwfe() at version 1.7.0 or later."
		)
	}

	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	N <- x$N
	T <- x$T
	R <- x$R
	d <- x$d
	p <- x$p
	sig_eps_sq <- x$sig_eps_sq
	cohort_probs_overall <- x$cohort_probs_overall
	X_final <- x$internal$X_final
	y_final <- x$internal$y_final
	theta_hat_full <- x$internal$theta_hat
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	first_inds <- getFirstInds(R = R, T = T)
	tes <- beta_hat[treat_inds]

	# Selected support in theta-space (slopes only; drop intercept)
	theta_hat_slopes <- theta_hat_full[2:(p + 1)]
	sel_feat_inds <- which(theta_hat_slopes != 0)
	sel_treat_inds_shifted <- which(theta_hat_slopes[treat_inds] != 0)
	theta_hat_treat_sel <- theta_hat_slopes[treat_inds][sel_treat_inds_shifted]

	# d_inv_treat: the treatment-block of D^{-1}, then restrict columns to selected
	d_inv_treat <- genInvTwoWayFusionTransformMat(
		n_vars = num_treats,
		first_inds = first_inds,
		R = R
	)
	if (length(sel_treat_inds_shifted) > 0) {
		d_inv_treat_sel <- d_inv_treat[,
			sel_treat_inds_shifted,
			drop = FALSE
		]
	} else {
		d_inv_treat_sel <- NULL
	}

	calc_ses <- contract$has_valid_ses
	gram_inv <- NULL
	if (calc_ses && length(sel_treat_inds_shifted) > 0) {
		res_gram <- getGramInv(
			N = N,
			T = T,
			X_final = X_final,
			sel_feat_inds = sel_feat_inds,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			calc_ses = TRUE
		)
		gram_inv <- res_gram$gram_inv
		calc_ses <- res_gram$calc_ses
	} else {
		calc_ses <- FALSE
	}

	# Cluster-robust sandwich (recomputed in theta-space on the selected support)
	sandwich_full <- NULL
	treat_block_mask <- NULL
	if (identical(se_type, "cluster") && calc_ses) {
		X_S <- X_final[, sel_feat_inds, drop = FALSE]
		treat_block_mask <- sel_feat_inds %in% treat_inds
		y_ <- y_final[seq_len(N * T)]
		ols_fit <- stats::lm.fit(cbind(1, X_S), y_)
		sandwich_full <- .compute_cluster_robust_sandwich(
			X_S = X_S,
			residuals = ols_fit$residuals,
			N = N,
			T = T
		)
	}

	max_event <- T - 2L
	event_times <- 0:max_event
	estimates <- numeric(length(event_times))
	ses <- numeric(length(event_times))
	n_cohorts <- integer(length(event_times))

	for (k in seq_along(event_times)) {
		e <- event_times[k]
		V_e <- which(seq_len(R) <= T - 1L - e)
		n_cohorts[k] <- length(V_e)
		if (length(V_e) == 0L) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}

		probs_Ve <- cohort_probs_overall[V_e]
		S_V <- sum(probs_Ve)
		if (S_V <= 0) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}
		weights_Ve <- probs_Ve / S_V

		# Point estimate: weighted sum of per-cell tes (which already
		# incorporates d_inv_treat times theta_hat under the FETWFE fit).
		psi_e_tes <- numeric(num_treats)
		for (j in seq_along(V_e)) {
			psi_e_tes[first_inds[V_e[j]] + e] <- weights_Ve[j]
		}
		estimates[k] <- as.numeric(crossprod(psi_e_tes, tes))

		if (!calc_ses) {
			ses[k] <- NA_real_
			next
		}

		# psi_e_theta on the selected theta-space treatment support:
		# psi_e_theta[col] = sum_{r in V_e} weights_Ve[r] * d_inv_treat_sel[idx(r,e), col]
		psi_e_theta <- numeric(length(sel_treat_inds_shifted))
		for (j in seq_along(V_e)) {
			r <- V_e[j]
			psi_e_theta <- psi_e_theta +
				weights_Ve[j] * d_inv_treat_sel[first_inds[r] + e, ]
		}

		# var_1(e)
		if (identical(se_type, "cluster")) {
			psi_full <- numeric(ncol(sandwich_full))
			psi_full[treat_block_mask] <- psi_e_theta
			var_1_e <- as.numeric(
				t(psi_full) %*% sandwich_full %*% psi_full
			)
		} else {
			var_1_e <- sig_eps_sq *
				as.numeric(
					t(psi_e_theta) %*% gram_inv %*% psi_e_theta
				) /
				(N * T)
		}

		# var_2(e) via the FETWFE-specific helper
		var_2_e <- .event_study_var2_fetwfe(
			e = e,
			V_e = V_e,
			theta_hat_treat_sel = theta_hat_treat_sel,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = cohort_probs_overall,
			first_inds = first_inds,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R
		)

		if (is_indep) {
			ses[k] <- sqrt(var_1_e + var_2_e)
		} else {
			ses[k] <- sqrt(var_1_e + var_2_e + 2 * sqrt(var_1_e * var_2_e))
		}
	}

	.assemble_event_study_df(event_times, n_cohorts, estimates, ses, z)
}

#' Per-event-time second-variance term for FETWFE event-study
#'
#' Structurally parallel to `getSecondVarTermDataApp` but with the
#' cohort-block time-averaging Jacobian replaced by a single-row
#' selection at `idx(r, e)`, and `cohort_probs_overall` masked to zero
#' outside `V_e`. When `|V_e| = 1` the Jacobian rows vanish exactly and
#' `var_2(e) = 0`.
#'
#' Design choice: this helper is structurally parallel to
#' `getSecondVarTermDataApp` in `R/fetwfe_core.R` (the cohort-block time-
#' averaging is replaced by a single-row selection at `idx(r, e)`, and
#' `cohort_probs_overall` is masked to zero outside `V_e`). Both helpers
#' implement paper Theorem 6.3's Jacobian formula (`paper_arxiv.tex:2577-
#' 2592`) where the off-diagonal coefficient `J_{rs} = -pi_s / S^2` uses the
#' column-index marginal cohort probability. (Prior to v1.8.0 the off-
#' diagonal coefficient was indexed by the outer-loop row; see issue #46.)
#' @keywords internal
#' @noRd
.event_study_var2_fetwfe <- function(
	e,
	V_e,
	theta_hat_treat_sel,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	num_treats,
	N,
	T,
	R
) {
	# Masked cohort_probs_overall: zero outside V_e
	masked <- numeric(R)
	masked[V_e] <- cohort_probs_overall[V_e]
	S_V <- sum(masked)
	if (S_V <= 0 || length(V_e) <= 1L) {
		# |V_e| = 1: single-cohort weight is 1; Jacobian rows vanish.
		return(0)
	}

	# Multinomial Sigma_pi_hat on the masked cohort probabilities
	Sigma_pi_hat <- -outer(masked, masked)
	diag(Sigma_pi_hat) <- masked * (1 - masked)

	# Per-event-time Jacobian: R rows, length(sel_treat_inds_shifted) cols.
	# Rows for r not in V_e are zero (and are zero-killed by Sigma_pi_hat).
	# Diagonal: (S_V - pi_r)/S_V^2 * d_inv_treat_sel[idx(r, e), ]
	# Off-diagonal: subtract pi_{r'}/S_V^2 * d_inv_treat_sel[idx(r', e), ]
	# for r' in V_e \ {r}. Off-diagonal coefficient uses the COLUMN-index
	# pi_{r'}, matching paper Theorem 6.3.
	jacobian_e <- matrix(0, nrow = R, ncol = ncol(d_inv_treat_sel))
	for (r in V_e) {
		cons_r <- (S_V - masked[r]) / S_V^2
		jacobian_e[r, ] <- cons_r * d_inv_treat_sel[first_inds[r] + e, ]
		for (r_prime in setdiff(V_e, r)) {
			cons_off <- masked[r_prime] / S_V^2
			jacobian_e[r, ] <- jacobian_e[r, ] -
				cons_off * d_inv_treat_sel[first_inds[r_prime] + e, ]
		}
	}

	T *
		as.numeric(
			t(theta_hat_treat_sel) %*%
				t(jacobian_e) %*%
				Sigma_pi_hat %*%
				jacobian_e %*%
				theta_hat_treat_sel
		) /
		(N * T)
}

#' Assemble the event-study output data frame
#' @keywords internal
#' @noRd
.assemble_event_study_df <- function(
	event_times,
	n_cohorts,
	estimates,
	ses,
	z
) {
	ci_low <- estimates - z * ses
	ci_high <- estimates + z * ses
	p_value <- .compute_p_values(estimates, ses)
	out <- data.frame(
		event_time = as.integer(event_times),
		n_cohorts = as.integer(n_cohorts),
		estimate = estimates,
		se = ses,
		ci_low = ci_low,
		ci_high = ci_high,
		p_value = p_value
	)
	class(out) <- c("fetwfe_event_study", "data.frame")
	out
}

#' @export
plot.fetwfe <- function(x, ...) {
	.plot_event_study(x, ...)
}

#' @export
plot.etwfe <- function(x, ...) {
	.plot_event_study(x, ...)
}

#' @export
plot.betwfe <- function(x, ...) {
	.plot_event_study(x, ...)
}

#' Render the event-study plot via ggplot2
#' @keywords internal
#' @noRd
.plot_event_study <- function(x, alpha = NULL, ...) {
	.check_for_plot(x)
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop(
			"Install ggplot2 to use plot.fetwfe() / plot.etwfe() / plot.betwfe(): install.packages('ggplot2')"
		)
	}
	es_df <- event_study(x, alpha = alpha)
	ggplot2::ggplot(
		es_df,
		ggplot2::aes(x = event_time, y = estimate)
	) +
		ggplot2::geom_hline(
			yintercept = 0,
			linetype = "dashed",
			color = "gray50"
		) +
		ggplot2::geom_errorbar(
			ggplot2::aes(ymin = ci_low, ymax = ci_high),
			width = 0.2
		) +
		ggplot2::geom_point(size = 2) +
		ggplot2::labs(
			x = "Event time (t - r)",
			y = "Pooled ATT",
			title = "Event-study estimates"
		)
}
