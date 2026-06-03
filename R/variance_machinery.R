# .combine_att_variance
#' @title Combine the two ATT variance pieces into the headline variance
#' @description Returns the tight Gaussian sum `att_var_1 + att_var_2`
#'   (independent-sample case, or same-data (Psi-IF) case under
#'   `se_type != "conservative"`) or the Cauchy-Schwarz upper bound
#'   `att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2)` (same-data
#'   non-(Psi-IF) case via `se_type = "conservative"`). Returns a VARIANCE;
#'   callers apply `sqrt()` for a standard error. Theorem (c') / (c),
#'   Assumption (Psi-IF); see `paper_arxiv.tex` ~ line 1233 and ~ line 2013.
#'   Consolidates five copies of this selection (#209): `getTeResultsOLS()`,
#'   `getTeResults2()`, `.build_variance_components()`, and the two
#'   `.event_study_*` helpers.
#' @param att_var_1,att_var_2 Numeric; the regression-coefficient and
#'   cohort-probability variance pieces (applied elementwise).
#' @param indep Logical; `TRUE` in the independent-sample case (each caller's
#'   local flag: `indep_probs` / `indep_counts_used` / `is_indep`).
#' @param se_type Character; the fit's `se_type`.
#' @return Numeric; the combined variance. Apply `sqrt()` for a standard error.
#' @keywords internal
#' @noRd
.combine_att_variance <- function(att_var_1, att_var_2, indep, se_type) {
	if (indep || !identical(se_type, "conservative")) {
		att_var_1 + att_var_2
	} else {
		att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2)
	}
}

# getTeResultsOLS
#' @title Calculate Overall ATT and its Standard Error for ETWFE
#' @description Computes the overall Average Treatment Effect on the Treated
#'   (ATT) by taking a weighted average of cohort-specific ATTs. It also
#'   calculates the standard error for this overall ATT, potentially considering
#'   variability from estimated cohort probabilities.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param cohort_tes Numeric vector; estimated ATTs for each of the `R` cohorts.
#' (Simple average of all of the estimated treatment effects for each cohort
#' across time.)
#' @param cohort_probs Numeric vector; weights for each cohort, typically
#'   estimated probabilities of belonging to cohort `g` conditional on being
#'   treated. Length `R`. Sums to 1.
#' @param psi_mat Numeric matrix; matrix where column `g` is `psi_g` (from
#'   `getCohortATTsFinal`). Dimensions: `length(sel_treat_inds_shifted)` x `R`.
#' @param gram_inv Numeric matrix; inverse of the Gram matrix for selected
#'   treatment effect features.
#' @param tes Numeric vector; all `num_treats` estimated treatment effects
#'   (original parameterization).
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort P(W=g). Length `R`.
#' @param calc_ses Logical; if `TRUE`, calculate standard errors.
#' @param indep_probs Logical; if `TRUE`, assumes `cohort_probs` (and
#'   `cohort_probs_overall`) were estimated from an independent sample, leading
#'   to the asymptotically-exact two-sample SE formula `sqrt(att_var_1 +
#'   att_var_2)`. When `FALSE`, the same-data SE combination is governed by
#'   `se_type`: `"default"` returns the tight Gaussian SE `sqrt(att_var_1 +
#'   att_var_2)` from Theorem (c$'$) under (Psi-IF); `"conservative"` returns
#'   the Cauchy-Schwarz upper bound from part (c).
#' @param se_type Character; one of `"default"` (tight Gaussian variance under
#'   (Psi-IF), Theorem (c$'$)), `"conservative"` (Cauchy-Schwarz upper bound
#'   from Theorem (c) for non-(Psi-IF) propensity estimators), or `"cluster"`
#'   (experimental unit-clustered Liang-Zeger sandwich SE; overrides the
#'   model-based `V_1` only, still uses the tight `V_1 + V_2` combination).
#'   Default `"default"`.
#' @param sandwich_full Numeric matrix (or `NULL`); the full cluster-robust
#'   sandwich variance on the selected support, as returned by
#'   `getCohortATTsFinal()`. Required when `se_type = "cluster"` and
#'   `calc_ses = TRUE`.
#' @param treat_block_mask Logical vector (or `NULL`) of length
#'   `nrow(sandwich_full)`, marking which columns correspond to treatment-effect
#'   features. Required when `se_type = "cluster"` and `calc_ses = TRUE`.
#' @return A list containing:
#'   \item{att_hat}{Numeric scalar; the estimated overall ATT.}
#'   \item{att_te_se}{Numeric scalar; the standard error for `att_hat`. NA if
#'     `calc_ses` is `FALSE`.}
#'   \item{att_te_se_no_prob}{Numeric scalar; standard error for `att_hat`
#'     ignoring variability from estimating cohort probabilities (i.e., only
#'     `att_var_1`). NA if `calc_ses` is `FALSE`.}
#'   \item{att_var_1}{Numeric scalar; the Kock-2013 regression-coefficient
#'     variance contribution. Maps to paper's `V_1 := sigma^2 * alpha' * Sigma^{-1} * alpha`
#'     up to `Var(hat_T_N) approx V_1 / N` (so `V_1 = N * att_var_1` at the
#'     per-unit scale paper notation uses). NA if `calc_ses` is `FALSE`.}
#'   \item{att_var_2}{Numeric scalar; the propensity-score variance
#'     contribution. Maps to paper's `V_2 := v_psi(beta_0, S)` analogously.
#'     NA if `calc_ses` is `FALSE`.}
#' @details The overall ATT (`att_hat`) is `cohort_tes %*% cohort_probs`.
#'   If `calc_ses` is `TRUE`:
#'   - `att_var_1` (variance from `theta_hat` estimation) is computed using
#'     `psi_att = psi_mat %*% cohort_probs` and `gram_inv` (or, under
#'     `se_type = "cluster"`, using the zero-padded `psi_att_full` against
#'     `sandwich_full`).
#'   - `att_var_2` (variance from cohort probability estimation) is computed by
#'     calling `getSecondVarTermOLS`. It is unchanged across `se_type` values.
#'   - When `indep_probs = TRUE`, the two pieces are independent and
#'     `att_te_se = sqrt(att_var_1 + att_var_2)` is asymptotically exact.
#'   - When `indep_probs = FALSE` (the common same-data case), the
#'     combination depends on `se_type`:
#'     * `"default"` (or `"cluster"`) uses the tight Gaussian variance
#'       `sqrt(att_var_1 + att_var_2)`, asymptotically exact under the paper's
#'       Theorem `te.asym.norm.thm`(c$'$) and Assumption (Psi-IF). Paper line
#'       1233 onwards. (Psi-IF) is satisfied by the package's default cohort
#'       sample-proportions estimator `hat_pi_g = N_g / N`, by multinomial
#'       logit, by any GLM on `W | X`, and by kernel/series regression of
#'       `1{W = g}` on `X`.
#'     * `"conservative"` uses the Cauchy-Schwarz upper bound
#'       `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
#'       Theorem (c) of the paper. This is the right tool only when the
#'       propensity-score estimator violates (Psi-IF) -- e.g.,
#'       Robins-Rotnitzky-augmented doubly-robust estimators that augment the
#'       propensity score with outcome residuals -- which the package does not
#'       currently implement.
#'   - `att_te_se_no_prob` is `sqrt(att_var_1)`.
#' @keywords internal
#' @noRd
getTeResultsOLS <- function(
	# model,
	sig_eps_sq,
	N,
	T,
	R,
	num_treats,
	cohort_tes,
	cohort_probs,
	psi_mat,
	gram_inv,
	tes,
	cohort_probs_overall,
	calc_ses,
	indep_probs = FALSE,
	se_type = "default",
	sandwich_full = NULL,
	treat_block_mask = NULL
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)

	att_hat <- as.numeric(cohort_tes %*% cohort_probs)

	if (calc_ses) {
		stopifnot(nrow(psi_mat) == length(tes))
		stopifnot(nrow(psi_mat) <= num_treats)
		stopifnot(ncol(psi_mat) == R)

		# Get ATT standard error
		# first variance term: convergence of theta
		psi_att <- psi_mat %*% cohort_probs

		if (identical(se_type, "cluster")) {
			stopifnot(!is.null(sandwich_full))
			stopifnot(!is.null(treat_block_mask))
			stopifnot(is.logical(treat_block_mask))
			stopifnot(length(treat_block_mask) == nrow(sandwich_full))
			stopifnot(sum(treat_block_mask) == nrow(psi_mat))

			psi_att_full <- numeric(length(treat_block_mask))
			psi_att_full[treat_block_mask] <- psi_att

			# Issue #84 item 9: floor the cluster-sandwich quadratic form
			# at zero. See the matching guard in `getTeResults2`
			# (same file, R/variance_machinery.R) for rationale. Issue
			# #139 layers a two-tier (warning / error) diagnostic on top
			# of the floor via `.floor_cluster_quad()`; on well-conditioned
			# data the behavior is unchanged.
			att_var_1 <- .floor_cluster_quad(
				as.numeric(
					t(psi_att_full) %*% sandwich_full %*% psi_att_full
				),
				"getTeResultsOLS/att_var_1"
			)
		} else {
			att_var_1 <- sig_eps_sq *
				as.numeric(t(psi_att) %*% gram_inv %*% psi_att) /
				(N * T)
		}

		stopifnot(length(tes) <= num_treats)

		stopifnot(length(tes) == nrow(psi_mat))

		# Second variance term: convergence of cohort membership probabilities
		att_var_2 <- getSecondVarTermOLS(
			psi_mat = psi_mat,
			tes = tes,
			cohort_probs_overall = cohort_probs_overall,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R
		)

		# Combine the two variance pieces. The independent-sample case
		# (`indep_probs = TRUE`) and the same-data (Psi-IF)-satisfying case
		# (`se_type %in% c("default", "cluster")`) both use the tight
		# Gaussian formula `att_var_1 + att_var_2`. The conservative
		# Cauchy-Schwarz fallback is reserved for the same-data
		# non-(Psi-IF) case via `se_type = "conservative"`. See
		# Theorem `te.asym.norm.thm`(c$'$) and Assumption (Psi-IF)
		# (`paper_arxiv.tex` ~ line 1233 and ~ line 2013).
		att_te_se <- sqrt(
			.combine_att_variance(att_var_1, att_var_2, indep_probs, se_type)
		)

		att_te_se_no_prob <- sqrt(att_var_1)
	} else {
		att_te_se <- NA
		att_te_se_no_prob <- NA
		att_var_1 <- NA
		att_var_2 <- NA
	}

	return(list(
		att_hat = att_hat,
		att_te_se = att_te_se,
		att_te_se_no_prob = att_te_se_no_prob,
		att_var_1 = att_var_1,
		att_var_2 = att_var_2
	))
}

# getSecondVarTermOLS
#' @title Calculate Second Variance Term for ATT Standard Error for ETWFE
#' @description Computes the second component of the variance for the Average
#'   Treatment Effect on the Treated (ATT). This component accounts for the
#'   variability due to the estimation of cohort membership probabilities.
#' @param psi_mat Numeric matrix; a matrix where each column `g` is the `psi_g`
#'   vector used in calculating the ATT for cohort `g`. Dimensions:
#'   `length(sel_treat_inds_shifted)` x `R`.
#' @param tes Numeric vector; the estimated treatment effects for all
#'   `num_treats` possible cohort-time combinations.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort (P(W=g)). Length `R`.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @return A numeric scalar representing the second variance component for the
#'   ATT.
#' @details This function calculates `Sigma_pi_hat`, the covariance matrix of
#'   the cohort assignment indicators, and a Jacobian matrix. These are then
#'   combined with `theta_hat_treat_sel` to compute the variance term as
#'   `T * t(theta_hat_treat_sel) %*% t(jacobian_mat) %*% Sigma_pi_hat %*% jacobian_mat %*%
#'   theta_hat_treat_sel / (N * T)`. The construction of the Jacobian involves averaging parts of
#'   `d_inv_treat_sel` corresponding to different cohorts.
#' @keywords internal
#' @noRd
getSecondVarTermOLS <- function(
	psi_mat,
	tes,
	cohort_probs_overall,
	num_treats,
	N,
	T,
	R
) {
	stopifnot(length(tes) == nrow(psi_mat))
	# Get Sigma_pi_hat, the (sample-estimated) covariance matrix for the
	# sample proportions (derived from the multinomial distribution).
	Sigma_pi_hat <- .multinomial_cov(cohort_probs_overall[1:R])

	stopifnot(nrow(Sigma_pi_hat) == R)
	stopifnot(ncol(Sigma_pi_hat) == R)

	# Construct Jacobian matrix corresponding to the mapping from the
	# individual cohort probabilities to the proportions calculated for the ATT.
	# See Proof of Theorem 6.1 for form. Consolidated into `.build_jacobian()`
	# (#192, WORKFLOW_LESSONS §14 Class A); the OLS-family path (d_inv = NULL)
	# returns the R x R column-indexed Jacobian byte-identically.
	jacobian_mat <- .build_jacobian(
		cohort_probs_overall = cohort_probs_overall,
		R = R,
		d_inv_treat_sel = NULL,
		mode = "scalar_block_avg"
	)

	stopifnot(nrow(psi_mat) <= num_treats)
	stopifnot(ncol(psi_mat) == R)

	stopifnot(length(tes) == nrow(psi_mat))

	# Finally, calculate variance term for ATT from Sigma_pi_hat, Jacobian
	# matrix, psi_mat, and beta_hat. See proof of Theorem D.2 for details.

	# Issue #127: floor `att_var_2` at zero, mirroring the `att_var_1`
	# floor PR #111 added at the analogous sibling sites in this file. The
	# `Sigma_pi_hat` quadratic form is PSD in exact arithmetic, so
	# `att_var_2 >= 0` always — but floating-point cancellation can leave
	# a near-zero value marginally negative, NaN-ing the conservative
	# overall-ATT-SE branch `2 * sqrt(att_var_1 * att_var_2)` downstream.
	att_var_2 <- max(
		T *
			as.numeric(
				t(tes) %*%
					psi_mat %*%
					t(jacobian_mat) %*%
					Sigma_pi_hat %*%
					jacobian_mat %*%
					t(psi_mat) %*%
					tes
			) /
			(N * T),
		0
	)

	return(att_var_2)
}


# getPsiRUnfused
#' @title Calculate Psi Vector for Cohort ATT (Unfused Case)
#' @description Computes the `psi_g` vector for a specific cohort `g` when no
#'   fusion penalization is applied to the treatment effects (or when calculating
#'   SEs as if it were an OLS on selected variables). This vector is used in
#'   standard error calculations for the cohort's Average Treatment Effect on
#'   the Treated (ATT). Specifically, `psi_g` places the constant
#'   `1 / k_full` (where `k_full = last_ind_g - first_ind_g + 1` is the
#'   cohort's full treatment-block size) at each position corresponding to a
#'   selected treatment effect in cohort `g`, and zero elsewhere --- so that
#'   `t(psi_g) %*% theta_hat_treat_sel == cohort_tes[g]`.
#' @param first_ind_g Integer; the index of the first treatment effect for
#'   cohort `g` within the `num_treats` block of treatment effects.
#' @param last_ind_g Integer; the index of the last treatment effect for
#'   cohort `g` within the `num_treats` block.
#' @param sel_treat_inds_shifted Integer vector; indices of all selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @return A numeric vector `psi_g` of length
#'   `length(sel_treat_inds_shifted)`. Positions corresponding to selected
#'   treatment effects in cohort `g` carry `1 / k_full`; other positions are
#'   zero. Returns a zero vector if no coefficients in cohort `g`'s block are
#'   selected.
#' @details The function identifies which of `sel_treat_inds_shifted` fall
#'   within the cohort-g block `[first_ind_g, last_ind_g]`. For those
#'   identified indices, `psi_g` is set to `1 / k_full`, where
#'   `k_full = last_ind_g - first_ind_g + 1` is the cohort's full
#'   treatment-block size. This normalization (rather than dividing by the
#'   *selected* count `k_sel`) is what makes
#'   `t(psi_g) %*% theta_hat_treat_sel == cohort_tes[g]`, since the
#'   cohort point estimate `cohort_tes[g] = mean(tes[first_ind_g:last_ind_g])`
#'   averages over the full block (unselected entries are exact zeros
#'   post-bridge). If no coefficients in cohort g's block are selected,
#'   `psi_g` is a zero vector. For ETWFE / `twfeCovs()` callers
#'   (which pass `sel_treat_inds_shifted = 1:num_treats`, so `k_sel = k_full`),
#'   this is bit-identical to the previous divide-by-sum normalization.
#' @keywords internal
#' @noRd
getPsiRUnfused <- function(
	first_ind_g,
	last_ind_g,
	sel_treat_inds_shifted
) {
	k_full <- last_ind_g - first_ind_g + 1
	stopifnot(k_full >= 1)

	which_inds_ir <- sel_treat_inds_shifted %in% (first_ind_g:last_ind_g)

	psi_g <- rep(0, length(sel_treat_inds_shifted))

	if (sum(which_inds_ir) > 0) {
		inds_g <- which(which_inds_ir)

		stopifnot(is.integer(inds_g) | is.numeric(inds_g))
		stopifnot(identical(inds_g, as.integer(round(inds_g))))
		stopifnot(length(inds_g) >= 1)
		stopifnot(length(inds_g) == length(unique(inds_g)))
		stopifnot(length(inds_g) <= length(sel_treat_inds_shifted))
		stopifnot(all(inds_g %in% 1:length(sel_treat_inds_shifted)))

		psi_g[inds_g] <- 1 / k_full
	}

	return(psi_g)
}


#' @title Compute the unit-clustered sandwich variance on the selected support
#' @description Internal helper used when `se_type = "cluster"`. Computes the
#'   full Liang-Zeger CR1 sandwich variance matrix
#'   `(X_S' X_S)^{-1} %*% sum_i X_{i.S}' eps_i eps_i' X_{i.S} %*% (X_S' X_S)^{-1}`
#'   with units `i = 1, ..., N` as clusters and an `N/(N-1)` cluster-count
#'   adjustment (matching `sandwich::vcovCL(cadjust = TRUE, type = "HC0")`).
#'   The design `X_S` is internally centered (consistent with how
#'   `getGramInv` centers the bread). The residuals must come from OLS on
#'   the selected support, not bridge residuals; the two differ in finite
#'   samples and only the OLS residuals are justified under the
#'   oracle-property argument.
#' @param X_S Numeric matrix, NT-by-p_S; the design matrix restricted to the
#'   selected support, in the coordinate system the regression was solved in
#'   (GLS-transformed for ETWFE/twfeCovs/BETWFE in beta-space; fusion-then-
#'   GLS-transformed for FETWFE in theta-space). Will be centered internally.
#' @param residuals Numeric vector of length NT; residuals from OLS on the
#'   selected support.
#' @param N Integer; number of units.
#' @param T Integer; number of time periods (panel is balanced).
#' @return A p_S-by-p_S symmetric matrix giving the cluster-robust sandwich
#'   variance on the selected support.
#' @references
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#'   generalized linear models. *Biometrika* 73(1), 13-22.
#' Bertrand, M., Duflo, E., & Mullainathan, S. (2004). How much should we
#'   trust differences-in-differences estimates? *Quarterly Journal of
#'   Economics* 119(1), 249-275.
#' @keywords internal
#' @noRd
.compute_cluster_robust_sandwich <- function(X_S, residuals, N, T) {
	stopifnot(nrow(X_S) == N * T)
	stopifnot(length(residuals) == N * T)
	X_S_centered <- scale(X_S, center = TRUE, scale = FALSE)
	p_S <- ncol(X_S_centered)
	# Issue #84 item 10: convert the LAPACK rank-deficiency error from
	# `solve(crossprod(X_S_centered))` into the same user-facing message
	# that `getGramInv()` emits, so a cluster-SE fit on a rank-deficient
	# selected support produces a single coherent surface across both
	# variance routes rather than an obscure "Lapack routine dgesv: system
	# is exactly singular" trace.
	gram_inv_full <- tryCatch(
		solve(crossprod(X_S_centered)),
		error = function(e) {
			stop(
				"Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated."
			)
		}
	)
	# Vectorized assembly of the cluster meat. The original N-loop built up
	#   meat = sum_i (X_i' eps_i) (X_i' eps_i)'
	# one unit at a time. Stacking the column vectors `X_i' eps_i` as the
	# rows of an N-by-p_S matrix `XEps` lets us write the same sum as
	#   meat = t(XEps) %*% XEps == crossprod(XEps).
	# Each row of `XEps` is itself a sum over the unit's T rows of the
	# response-weighted design (`X_S_centered[rows_i, ] * eps_i`); a single
	# `rowsum()` over a unit-id grouping vector does that aggregation across
	# all units at once (issue #84 item 18). Equivalent to the pre-vectorize
	# loop bit-for-bit modulo floating-point summation order.
	weighted <- X_S_centered * residuals
	unit_id <- rep(seq_len(N), each = T)
	XEps <- rowsum(weighted, group = unit_id, reorder = FALSE)
	meat <- crossprod(XEps)
	cadjust <- N / (N - 1)
	cadjust * gram_inv_full %*% meat %*% gram_inv_full
}

#' @title Assemble the cluster-robust sandwich for the OLS-selected support
#' @description
#' Wraps the 4-step assembly ritual used by `getCohortATTsFinal()` and the
#' two `.event_study_*` dispatchers
#' (`.event_study_etwfe_betwfe`, `.event_study_fetwfe`): extract `X_S`
#' from `X_final` (optionally filtered by `sel_feat_inds`), run
#' `stats::lm.fit(cbind(1, X_S), y_final[seq_len(N * T)])`, compute the
#' cluster-robust sandwich via `.compute_cluster_robust_sandwich()`, and
#' build the corresponding `treat_block_mask`. Returns both as a list so
#' the caller can use them downstream.
#'
#' Resolves GitHub #78. The drift class this addresses is the v1.9.5
#' `10e-6` threshold typo: same logic written across 4 siblings, one
#' site could get a fix while the others quietly didn't.
#' @param X_final Numeric matrix. The post-GLS, post-augmentation design.
#' @param y_final Numeric vector. The post-GLS response (length >=
#'   `N * T`; only the first `N * T` entries are used).
#' @param N Integer. Number of units.
#' @param T Integer. Number of time periods.
#' @param treat_inds Integer vector. The (1-based) indices of treatment
#'   features in the full `X_final` columns. Used in both branches'
#'   `treat_block_mask` construction.
#' @param sel_feat_inds Integer vector or `NULL`. If `NULL` (the default),
#'   the full `X_final` is used as `X_S` and `treat_block_mask` is built
#'   via `logical(ncol(X_S)); [treat_inds] <- TRUE`. If a vector, `X_S`
#'   is `X_final[, sel_feat_inds, drop = FALSE]` and `treat_block_mask`
#'   is `sel_feat_inds %in% treat_inds`.
#' @return A list with two elements:
#'   \describe{
#'     \item{`sandwich_full`}{Numeric matrix; the Liang-Zeger
#'       cluster-robust sandwich from `.compute_cluster_robust_sandwich()`.}
#'     \item{`treat_block_mask`}{Logical vector of length `ncol(X_S)`;
#'       TRUE at positions corresponding to treatment features.}
#'   }
#' @keywords internal
#' @noRd
.assemble_cluster_robust_sandwich <- function(
	X_final,
	y_final,
	N,
	T,
	treat_inds,
	sel_feat_inds = NULL
) {
	X_S <- if (is.null(sel_feat_inds)) {
		X_final
	} else {
		X_final[, sel_feat_inds, drop = FALSE]
	}
	y_ <- y_final[seq_len(N * T)]
	ols_fit <- stats::lm.fit(cbind(1, X_S), y_)
	sandwich_full <- .compute_cluster_robust_sandwich(
		X_S = X_S,
		residuals = ols_fit$residuals,
		N = N,
		T = T
	)
	treat_block_mask <- if (is.null(sel_feat_inds)) {
		mask <- logical(ncol(X_S))
		mask[treat_inds] <- TRUE
		mask
	} else {
		sel_feat_inds %in% treat_inds
	}
	list(sandwich_full = sandwich_full, treat_block_mask = treat_block_mask)
}


# getTeResults2
#' @title Calculate Overall ATT and its Standard Error (Version 2)
#' @description Computes the overall Average Treatment Effect on the Treated
#'   (ATT) by taking a weighted average of cohort-specific ATTs. It also
#'   calculates the standard error for this overall ATT, potentially considering
#'   variability from estimated cohort probabilities.
#'
#' The overall Average Treatment Effect on the Treated is
#' \deqn{\widehat{\text{ATT}}
#'       \;=\;
#'       \sum_{g=1}^{R}\widehat{\tau}_{\text{ATT},g}\,
#'                     \widehat{\pi}_{g\mid\tau},}
#' a weighted mean of cohort-specific ATTs with weights
#' \(\widehat{\pi}_{g\mid\tau}=N_g/N_\tau\).
#' The variance has **two additive pieces**
#'
#' * *Term 1* - randomness from \(\widehat\theta\):
#'   \(\sigma^{2}(NT)^{-1}\psi_{\text{att}}^{\!\top}\widehat G^{-1}\psi_{\text{att}}\)
#'   with \(\psi_{\text{att}}=\Psi\,\widehat\pi\) and
#'   \(\Psi=(\psi_1,\dots,\psi_R)\).
#' * *Term 2* - randomness from \(\widehat\pi\) itself.
#'   It equals \(N^{-1}\theta_{\text{sel}}^{\!\top}J^{\!\top}
#'   \widehat\Sigma_\pi J\theta_{\text{sel}}\) when cohort counts and outcomes
#'   come from independent samples, and it enters a conservative sum of
#'   variances otherwise.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param cohort_tes Numeric vector; estimated ATTs for each of the `R` cohorts.
#' (Simple average of all of the estimated treatment effects for each cohort
#' across time.)
#' @param cohort_probs Numeric vector; weights for each cohort, typically
#'   estimated probabilities of belonging to cohort `g` conditional on being
#'   treated. Length `R`. Sums to 1.
#' @param psi_mat Numeric matrix; matrix where column `g` is `psi_g` (from
#'   `getCohortATTsFinal`). Dimensions: `length(sel_treat_inds_shifted)` x `R`.
#' @param gram_inv Numeric matrix; inverse of the Gram matrix for selected
#'   treatment effect features.
#' @param sel_treat_inds_shifted Integer vector; indices of selected treatment
#'   effects within the `num_treats` block (shifted to start from 1).
#' @param d_inv_treat_sel Numeric matrix; block of the inverse fusion matrix for
#'   selected treatment effects.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort P(W=g). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in
#'   transformed (fused) space for selected treatment effects.
#' @param calc_ses Logical; if `TRUE`, calculate standard errors.
#' @param indep_probs Logical; if `TRUE`, assumes `cohort_probs` (and
#'   `cohort_probs_overall`) were estimated from an independent sample, leading
#'   to a different SE formula (sum of variances) compared to when they are
#'   estimated from the same sample (conservative SE including a covariance term).
#' @param se_type Character; one of `"default"` (tight Gaussian variance under
#'   Assumption (Psi-IF), Theorem (c$'$)), `"conservative"` (Cauchy-Schwarz
#'   upper bound from Theorem (c) for non-(Psi-IF) propensity estimators), or
#'   `"cluster"` (experimental unit-clustered Liang-Zeger sandwich SE;
#'   overrides the model-based `V_1` only, still uses the tight `V_1 + V_2`
#'   combination). Default `"default"`.
#' @param sandwich_full Numeric matrix (or `NULL`); the full cluster-robust
#'   sandwich variance on the selected support, as returned by
#'   `getCohortATTsFinal()`. Required when `se_type = "cluster"` and
#'   `calc_ses = TRUE`.
#' @param treat_block_mask Logical vector (or `NULL`) of length
#'   `nrow(sandwich_full)`, marking which selected columns correspond to
#'   treatment-effect features. Required when `se_type = "cluster"` and
#'   `calc_ses = TRUE`.
#' @return A list containing:
#'   \item{att_hat}{Numeric scalar; the estimated overall ATT.}
#'   \item{att_te_se}{Numeric scalar; the standard error for `att_hat`. NA if
#'     `calc_ses` is `FALSE`.}
#'   \item{att_te_se_no_prob}{Numeric scalar; standard error for `att_hat`
#'     ignoring variability from estimating cohort probabilities (i.e., only
#'     `att_var_1`). NA if `calc_ses` is `FALSE`.}
#'   \item{att_var_1}{Numeric scalar; the Kock-2013 regression-coefficient
#'     variance contribution. Maps to paper's `V_1 := sigma^2 * alpha' * Sigma^{-1} * alpha`
#'     (per-unit-scale `V_1 = N * att_var_1`). NA if `calc_ses` is `FALSE`.}
#'   \item{att_var_2}{Numeric scalar; the propensity-score variance
#'     contribution. Maps to paper's `V_2 := v_psi(beta_0, S)`. NA if
#'     `calc_ses` is `FALSE`.}
#' @details The overall ATT (`att_hat`) is `cohort_tes %*% cohort_probs`.
#'   If `calc_ses` is `TRUE`:
#'   - `att_var_1` (variance from `theta_hat` estimation) is computed using
#'     `psi_att = psi_mat %*% cohort_probs` and `gram_inv` (or, under
#'     `se_type = "cluster"`, via the zero-padded `psi_att_full` against
#'     `sandwich_full`).
#'   - `att_var_2` (variance from cohort probability estimation) is computed by
#'     calling `getSecondVarTermDataApp`. It is unchanged across `se_type`
#'     values.
#'   - When `indep_probs = TRUE`, the two pieces are independent and
#'     `att_te_se = sqrt(att_var_1 + att_var_2)` is asymptotically exact.
#'   - When `indep_probs = FALSE` (the common same-data case), the
#'     combination depends on `se_type`:
#'     * `"default"` (or `"cluster"`) uses the tight Gaussian variance
#'       `sqrt(att_var_1 + att_var_2)`, asymptotically exact under the paper's
#'       Theorem `te.asym.norm.thm`(c$'$) and Assumption (Psi-IF). Paper line
#'       1233 onwards.
#'     * `"conservative"` uses the Cauchy-Schwarz upper bound
#'       `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
#'       Theorem (c). This is the right tool only when the propensity-score
#'       estimator violates (Psi-IF) -- e.g., Robins-Rotnitzky-augmented
#'       doubly-robust estimators -- which the package does not currently
#'       implement.
#'   - `att_te_se_no_prob` is `sqrt(att_var_1)`.
#'
#' Prior to v1.12.0, the same-data path returned the conservative
#' Cauchy-Schwarz upper bound as the default. The headline change in v1.12.0
#' is that the tight Gaussian formula is now the default; the Cauchy-Schwarz
#' bound is reachable via `se_type = "conservative"`.
#'
#' All matrices required for Term 2 are produced upstream:
#' * `psi_mat` from \code{getCohortATTsFinal()}
#' * `d_inv_treat_sel` from the same routine
#' * the Jacobian \(J\) is built in
#'   \code{getSecondVarTermDataApp()} using `d_inv_treat_sel`
#' @seealso \code{getSecondVarTermDataApp()}
#' @keywords internal
#' @noRd
getTeResults2 <- function(
	# model,
	sig_eps_sq,
	N,
	T,
	R,
	num_treats,
	cohort_tes,
	cohort_probs,
	psi_mat,
	gram_inv,
	sel_treat_inds_shifted,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	calc_ses,
	indep_probs = FALSE,
	se_type = "default",
	sandwich_full = NULL,
	treat_block_mask = NULL
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)

	stopifnot(all(!is.na(cohort_probs)))

	# point estimate ATT
	att_hat <- as.numeric(cohort_tes %*% cohort_probs)

	if (calc_ses) {
		stopifnot(all(!is.na(psi_mat)))
		stopifnot(!is.na(sig_eps_sq))
		stopifnot(all(!is.na(gram_inv)))
		# Get ATT standard error

		# first variance term: convergence of theta
		psi_att <- psi_mat %*% cohort_probs

		if (identical(se_type, "cluster")) {
			stopifnot(!is.null(sandwich_full))
			stopifnot(!is.null(treat_block_mask))
			stopifnot(is.logical(treat_block_mask))
			stopifnot(length(treat_block_mask) == nrow(sandwich_full))
			stopifnot(sum(treat_block_mask) == nrow(psi_mat))

			psi_att_full <- numeric(length(treat_block_mask))
			psi_att_full[treat_block_mask] <- psi_att

			# Issue #84 item 9: floor the cluster-sandwich quadratic form
			# at zero. The Liang-Zeger sandwich is PSD in exact arithmetic,
			# so `t(psi) %*% sandwich %*% psi >= 0` always — but the
			# now-vectorized `rowsum` assembly (see item 18) and the
			# subsequent matrix triple-product accumulate enough
			# floating-point rounding that a near-zero quadratic form
			# could come out negative by a few ulps, NaN-ing downstream
			# sqrt() into a baffling "NA SE" surface. Floor at zero so
			# any rounding-induced negative collapses to the
			# mathematically correct value. Issue #139 layers a two-tier
			# (warning / error) diagnostic on top of the floor via
			# `.floor_cluster_quad()`; on well-conditioned data the
			# behavior is unchanged.
			att_var_1 <- .floor_cluster_quad(
				as.numeric(
					t(psi_att_full) %*% sandwich_full %*% psi_att_full
				),
				"getTeResults2/att_var_1"
			)
		} else {
			att_var_1 <- sig_eps_sq *
				as.numeric(t(psi_att) %*% gram_inv %*% psi_att) /
				(N * T)
		}

		stopifnot(!is.na(att_var_1))

		# Second variance term: convergence of cohort membership probabilities
		att_var_2 <- getSecondVarTermDataApp(
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = cohort_probs_overall,
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R
		)

		stopifnot(!is.na(att_var_2))

		# Combine the two variance pieces. The independent-sample case
		# (`indep_probs = TRUE`) and the same-data (Psi-IF)-satisfying case
		# (`se_type %in% c("default", "cluster")`) both use the tight
		# Gaussian formula `att_var_1 + att_var_2`. The conservative
		# Cauchy-Schwarz fallback is reserved for the same-data
		# non-(Psi-IF) case via `se_type = "conservative"`. See
		# Theorem `te.asym.norm.thm`(c$'$) and Assumption (Psi-IF)
		# (`paper_arxiv.tex` ~ line 1233 and ~ line 2013). Prior to v1.12.0
		# the same-data path defaulted to the Cauchy-Schwarz bound regardless
		# of (Psi-IF); v1.12.0 inverted that default.
		att_te_se <- sqrt(
			.combine_att_variance(att_var_1, att_var_2, indep_probs, se_type)
		)

		att_te_se_no_prob <- sqrt(att_var_1)
	} else {
		att_te_se <- NA
		att_te_se_no_prob <- NA
		att_var_1 <- NA
		att_var_2 <- NA
	}

	return(list(
		att_hat = att_hat,
		att_te_se = att_te_se,
		att_te_se_no_prob = att_te_se_no_prob,
		att_var_1 = att_var_1,
		att_var_2 = att_var_2
	))
}


# getPsiRFused
#' @title Calculate Psi Vector and D-inverse Block for Cohort ATT (Fused Case)
#' @description Computes the `psi_g` vector and the relevant block of the
#'   inverse fusion transformation matrix (`d_inv_treat_sel`) for a specific
#'   cohort `g`. These are used in standard error calculations for the cohort's
#'   Average Treatment Effect on the Treated (ATT) when fusion penalization
#'   has been applied.
#'
#' For a treated cohort \(g\) let
#' \(M_g=T-g+1\) be the number of post-treatment periods.
#' The CATT is the *simple average*
#' \(M_g^{-1}\sum_{t=g}^{T}\tau_{\text{ATT}}(g,t)\).
#' In the fused basis each \(\tau_{\text{ATT}}(g,t)\) is a row of
#' \(D^{(2)}(\mathcal R)^{-1}\widehat\theta\),
#' so the averaging weights are the **row-means** of the corresponding block
#' of the inverse fusion matrix.
#' This helper:
#'
#' * extracts those rows (`first_ind_g:last_ind_g`) from
#'   `d_inv_treat`,
#' * averages them column-wise to form \(\psi_g\),
#' * returns the block itself (`d_inv_treat_sel`) because later routines need
#'   it for probability-variance propagation.
#' @param first_ind_g Integer; the index of the first treatment effect parameter
#'   for cohort `g` within the original `num_treats` block (1-based).
#' @param last_ind_g Integer; the index of the last treatment effect parameter
#'   for cohort `g` within the original `num_treats` block (1-based).
#' @param sel_treat_inds_shifted Integer vector; indices (1-based) of the
#'   treatment effects that were selected by the model, relative to the start
#'   of the `num_treats` block. E.g., if original indices 5, 7 were selected
#'   from a block starting at index 1, this would be c(5, 7).
#' @param d_inv_treat Numeric matrix; the full inverse two-way fusion
#'   transformation matrix for all `num_treats` treatment effects. Dimensions:
#'   `num_treats` x `num_treats`.
#' @return A list containing:
#'   \item{psi_g}{Numeric vector. It's the column means of the sub-matrix of
#'     `d_inv_treat` corresponding to rows `first_ind_g:last_ind_g` and columns
#'     specified by `sel_treat_inds_shifted`. If `first_ind_g == last_ind_g`,
#'     it's just that specific row of `d_inv_treat` (subsetted by selected columns).}
#'   \item{d_inv_treat_sel}{Numeric matrix. The sub-matrix of `d_inv_treat`
#'     with rows `first_ind_g:last_ind_g` and columns corresponding to
#'     `sel_treat_inds_shifted`.}
#' @details `psi_g` effectively averages the rows of `d_inv_treat` (that correspond
#'   to cohort `g`'s treatment effects) for the columns that were actually
#'   selected by the model. `d_inv_treat_sel` is this specific block of the
#'   `d_inv_treat` matrix.
#'
#' * If only one transformed treatment coefficient was selected,
#'   both the mean and the returned block are forced to the correct
#'   1 x 1 or 1 x *`k`* shape so that higher-level code can
#'   `rbind()` the blocks without special cases.
#' @inheritParams getPsiRFused
#' @seealso \code{getCohortATTsFinal()}
#' @keywords internal
#' @noRd
getPsiRFused <- function(
	first_ind_g,
	last_ind_g,
	sel_treat_inds_shifted,
	d_inv_treat
) {
	stopifnot(length(sel_treat_inds_shifted) >= 0)
	stopifnot(last_ind_g >= first_ind_g)
	# Get psi vector: the part of D inverse that we need to look at is the
	# block corresponding to the treatment effect estimates, which is the
	# num_treats x num_treats matrix yielded by
	# genInvTwoWayFusionTransformMat(num_treats, first_inds).

	# Correct rows of matrix

	## psi_g := column-wise mean of those rows  (weights for average treatment
	## effect for cohort g)
	##
	## * If |S| > 1, result is a vector length |S|
	## * If |S| == 1, treat the scalar mean as length-1 vector

	if (last_ind_g > first_ind_g) {
		if (length(sel_treat_inds_shifted) > 1) {
			psi_g <- colMeans(d_inv_treat[
				first_ind_g:last_ind_g,
				sel_treat_inds_shifted
			])
		} else {
			psi_g <- mean(d_inv_treat[
				first_ind_g:last_ind_g,
				sel_treat_inds_shifted
			])
		}
		if (length(sel_treat_inds_shifted) == 1) {
			# Need to coerce this object to be a matrix with one column so it
			# works smoothly with rbind() later
			d_inv_treat_sel <- matrix(
				d_inv_treat[first_ind_g:last_ind_g, sel_treat_inds_shifted],
				ncol = 1
			)
		} else {
			d_inv_treat_sel <- d_inv_treat[
				first_ind_g:last_ind_g,
				sel_treat_inds_shifted
			]
		}
	} else {
		psi_g <- d_inv_treat[first_ind_g:last_ind_g, sel_treat_inds_shifted]
		# Since first_ind_g and last_ind_g are the same, need to coerce this
		# object to be a matrix with one row so that it works smoothly with
		# rbind() later

		## Block of D^{-1} used later for probability-variance term
		d_inv_treat_sel <- matrix(
			d_inv_treat[first_ind_g:last_ind_g, sel_treat_inds_shifted],
			nrow = 1
		)
	}

	stopifnot(is.matrix(d_inv_treat_sel))

	return(list(psi_g = psi_g, d_inv_treat_sel = d_inv_treat_sel))
}


# getSecondVarTermDataApp
#' @title Calculate Second Variance Term for ATT Standard Error (Data Application)
#' @description Computes the second component of the variance for the Average
#'   Treatment Effect on the Treated (ATT). This component accounts for the
#'   variability due to the estimation of cohort membership probabilities.
#' Computes **Term 2** of the overall-ATT variance:
#' \deqn{\frac{1}{N}\,
#'       \hat\theta_{\!\text{sel}}^{\!\top}
#'       J^{\!\top}\widehat\Sigma_\pi J
#'       \hat\theta_{\!\text{sel}},}
#' where
#' \itemize{
#'   \item \(\widehat\Sigma_\pi\) is the multinomial covariance of the
#'         cohort-count vector
#'         (\(\widehat\pi_g(1-\widehat\pi_g)\) on the diagonal,
#'         \(-\widehat\pi_g\widehat\pi_s\) off-diagonal).
#'   \item \(J\) is the Jacobian of the weighting function
#'         \(f_g(\pi)=\pi_g/\sum_{k}\pi_k\) evaluated at
#'         \(\widehat\pi\).  In matrix form
#'         \eqn{J_{rs}=
#'           \begin{cases}
#'             (1-\widehat\pi_g)/S^2,& g=s,\\[4pt]
#'             -\,\widehat\pi_g    /S^2,& g\neq s,
#'           \end{cases}} with \(S=\sum_k\widehat\pi_k\).
#'   \item Each column of `d_inv_treat_sel` is a selected **transformed**
#'         treatment coefficient; row-averaging over the rows belonging to
#'         cohort \(s\) produces the part of \(J\) that multiplies
#'         \(\theta_{\text{sel}}\).
#' }
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort (P(W=g)). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the `num_treats` block.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in the
#'   transformed (fused) space, corresponding only to the selected treatment
#'   effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param d_inv_treat_sel Numeric matrix; the relevant block of the inverse
#'   two-way fusion transformation matrix corresponding to selected treatment
#'   effects. Dimensions: `num_treats` (or fewer if selection occurs) x
#'   `length(sel_treat_inds_shifted)`.
#' @return A numeric scalar representing the second variance component for the
#'   ATT.
#' @details This function calculates `Sigma_pi_hat`, the covariance matrix of
#'   the cohort assignment indicators, and a Jacobian matrix. These are then
#'   combined with `theta_hat_treat_sel` to compute the variance term as
#'   `T * t(theta_hat_treat_sel) %*% t(jacobian_mat) %*% Sigma_pi_hat %*% jacobian_mat %*%
#'   theta_hat_treat_sel / (N * T)`. The construction of the Jacobian involves averaging parts of
#'   `d_inv_treat_sel` corresponding to different cohorts.
#'
#' It first reconstructs \(J\) from `d_inv_treat_sel` and the cohort
#' probabilities, then plugs everything into the quadratic form above and
#' finally rescales by \(T/(N T)=1/N\).
#' @inheritParams getSecondVarTermDataApp
#' @seealso \code{getTeResults2()}
#' @keywords internal
#' @noRd
getSecondVarTermDataApp <- function(
	sel_treat_inds_shifted,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	num_treats,
	N,
	T,
	R,
	d_inv_treat_sel
) {
	stopifnot(all(!is.na(d_inv_treat_sel)))
	stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
	stopifnot(length(theta_hat_treat_sel) == length(sel_treat_inds_shifted))

	# Get Sigma_pi_hat, the (sample-estimated) covariance matrix for the
	# sample proportions (derived from the multinomial distribution).
	Sigma_pi_hat <- .multinomial_cov(cohort_probs_overall[1:R])

	stopifnot(nrow(Sigma_pi_hat) == R)
	stopifnot(ncol(Sigma_pi_hat) == R)

	# Jacobian (paper Theorem 6.3, paper_arxiv.tex:2577-2592):
	##
	## J_{rs} =  (S-pi_g)/S^2      if g = s
	##           -pi_s   /S^2      if g != s    (off-diagonal uses COLUMN index s)
	##
	## where S := sum(cohort_probs_overall).
	##
	## Note: prior to v1.8.0 (issue #46) the off-diagonal coefficient was
	## indexed by the outer-loop row g instead of the column s; this matched
	## the textbook delta-method gradient only when cohort probabilities
	## were uniform. The fix uses the column index, matching the paper and
	## ETWFE's parallel `getSecondVarTermOLS()`.
	##
	## Each column block of d_inv_treat_sel belongs to a selected theta coordinate;
	## averaging the rows of cohorts gives the vector needed to multiply theta_sel.
	##
	## Consolidated into `.build_jacobian()` (#192, WORKFLOW_LESSONS §14 Class
	## A); the FETWFE scalar-block-averaging path (d_inv non-NULL) returns the
	## R x p_sel block-averaged Jacobian byte-identically.
	jacobian_mat <- .build_jacobian(
		cohort_probs_overall = cohort_probs_overall,
		R = R,
		d_inv_treat_sel = d_inv_treat_sel,
		mode = "scalar_block_avg",
		first_inds = first_inds,
		num_treats = num_treats
	)

	## variance term: theta_sel' J' sum_pi J theta_sel / N
	# Issue #127: floor `att_var_2` at zero, mirroring the `att_var_1`
	# floor PR #111 added at the analogous sibling sites in this file. The
	# `Sigma_pi_hat` quadratic form is PSD in exact arithmetic, so
	# `att_var_2 >= 0` always — but floating-point cancellation can leave
	# a near-zero value marginally negative, NaN-ing the conservative
	# overall-ATT-SE branch `2 * sqrt(att_var_1 * att_var_2)` downstream.
	att_var_2 <- max(
		T *
			as.numeric(
				t(theta_hat_treat_sel) %*%
					t(jacobian_mat) %*%
					Sigma_pi_hat %*%
					jacobian_mat %*%
					theta_hat_treat_sel
			) /
			(N * T),
		0
	)

	return(att_var_2)
}


# .build_jacobian
#' @title Unified multinomial-weight Jacobian builder
#' @description
#' Single home for the cohort-weight-function Jacobian \eqn{J} used by the
#' package's Term-2 (cohort-probability) variance machinery. Consolidates the
#' three pre-#192 inline Jacobian-build sites
#' (`getSecondVarTermOLS()`, `getSecondVarTermDataApp()`, and
#' `.event_study_var2_fetwfe()` in `R/event_study.R`) plus the new
#' `simultaneousCIs()` call site into one helper (WORKFLOW_LESSONS section 14 Class A:
#' the third copy triggers a refactor). All four sites implement the same
#' delta-method gradient of \eqn{f_g(\pi)=\pi_g/\sum_k \pi_k}; the column-index
#' off-diagonal coefficient \eqn{J_{rs}=-\pi_s/S^2} matches paper Theorem 6.3
#' (`paper_arxiv.tex:2577-2592`). Prior to v1.8.0 the off-diagonal used the
#' outer-loop row index; see issue #46.
#' @param cohort_probs_overall Numeric vector of length `R`; marginal cohort
#'   probabilities P(W = g). For `mode = "per_effect_masked"` the caller passes
#'   the full (un-masked) vector; this helper masks to `V_e` internally.
#' @param R Integer; number of treated cohorts.
#' @param d_inv_treat_sel Numeric matrix (`num_treats` x `p_sel`) or `NULL`.
#'   `NULL` selects the OLS-family R x R identity path (reproduces
#'   `getSecondVarTermOLS()`); non-`NULL` selects the FETWFE R x p_sel path
#'   (block-averaged for `scalar_block_avg`, single-row for
#'   `per_effect_masked`).
#' @param mode Character; one of `"scalar_block_avg"` (cohort-block averaging,
#'   used by the two overall-ATT Term-2 helpers) or `"per_effect_masked"`
#'   (per-event-time single-row masked construction, used by the event-study
#'   Term-2 helper and `simultaneousCIs()`).
#' @param V_e Integer vector; required for `mode = "per_effect_masked"`. The
#'   cohorts valid at the event time (`which(cohort_offsets_int <= T - e)`).
#' @param first_inds Integer vector of length `R`; required for FETWFE paths.
#'   Index of the first treatment effect for each cohort within the
#'   `num_treats` block.
#' @param e Integer; required for `mode = "per_effect_masked"`. The event-time
#'   offset.
#' @param num_treats Integer; total number of base treatment effect parameters.
#'   Required for the FETWFE `scalar_block_avg` path (cohort-block indexing).
#' @return A numeric matrix: R x R when `d_inv_treat_sel` is `NULL`, otherwise
#'   R x `ncol(d_inv_treat_sel)`. For `mode = "per_effect_masked"` with
#'   `|V_e| <= 1` (or zero masked mass) the all-zero matrix is returned (the
#'   single-cohort weight is 1, so the Jacobian rows vanish).
#' @keywords internal
#' @noRd
.build_jacobian <- function(
	cohort_probs_overall,
	R,
	d_inv_treat_sel = NULL,
	mode = c("scalar_block_avg", "per_effect_masked"),
	V_e = NULL,
	first_inds = NULL,
	e = NULL,
	num_treats = NULL
) {
	mode <- match.arg(mode)

	if (mode == "scalar_block_avg") {
		stopifnot(length(cohort_probs_overall) == R)
		stopifnot(sum(cohort_probs_overall) < 1 - 1e-6)

		if (is.null(d_inv_treat_sel)) {
			# OLS-family path (reproduces getSecondVarTermOLS L252-276):
			# R x R Jacobian. Column g is -pi_g/S^2 everywhere except the
			# diagonal (S - pi_g)/S^2.
			jacobian_mat <- matrix(as.numeric(NA), nrow = R, ncol = R)
			for (g in 1:R) {
				col_g_val <- -cohort_probs_overall[g] /
					sum(cohort_probs_overall)^2
				jacobian_mat[, g] <- rep(col_g_val, R)
				cons_g <- (sum(cohort_probs_overall) -
					cohort_probs_overall[g]) /
					sum(cohort_probs_overall)^2
				jacobian_mat[g, g] <- cons_g
			}
			stopifnot(all(!is.na(jacobian_mat)))
			return(jacobian_mat)
		}

		# FETWFE path (reproduces getSecondVarTermDataApp L990-1041): R x
		# p_sel Jacobian. Row g averages the d_inv_treat_sel rows of each
		# cohort block, weighted by the delta-method gradient.
		stopifnot(!is.null(first_inds), !is.null(num_treats))
		p_sel <- ncol(d_inv_treat_sel)
		jacobian_mat <- matrix(as.numeric(NA), nrow = R, ncol = p_sel)

		sel_inds <- list()
		for (g in 1:R) {
			sel_inds[[g]] <- .cohort_block_inds(g, R, first_inds, num_treats)
			if (g > 1) {
				stopifnot(min(sel_inds[[g]]) > max(sel_inds[[g - 1]]))
				stopifnot(length(sel_inds[[g]]) <= length(sel_inds[[g - 1]]))
			}
		}
		stopifnot(all.equal(unlist(sel_inds), 1:num_treats))

		for (g in 1:R) {
			cons_g <- (sum(cohort_probs_overall) -
				cohort_probs_overall[g]) /
				sum(cohort_probs_overall)^2

			if (p_sel > 1) {
				jacobian_mat[g, ] <- cons_g *
					colMeans(d_inv_treat_sel[sel_inds[[g]], , drop = FALSE])
			} else {
				jacobian_mat[g, ] <- cons_g *
					mean(d_inv_treat_sel[sel_inds[[g]], , drop = FALSE])
			}

			for (g_double_prime in setdiff(1:R, g)) {
				cons_g_double_prime <- cohort_probs_overall[g_double_prime] /
					sum(cohort_probs_overall)^2

				jacobian_mat[g, ] <- jacobian_mat[g, ] -
					cons_g_double_prime *
						colMeans(d_inv_treat_sel[
							sel_inds[[g_double_prime]],
							,
							drop = FALSE
						])
			}
		}
		stopifnot(all(!is.na(jacobian_mat)))
		return(jacobian_mat)
	}

	# mode == "per_effect_masked" (reproduces .event_study_var2_fetwfe
	# L572-581): R x p_sel Jacobian for a specific event-time `e`'s V_e,
	# with cohort_probs_overall masked to zero outside V_e. Rows for
	# g not in V_e are left at zero.
	stopifnot(
		!is.null(V_e),
		!is.null(first_inds),
		!is.null(e),
		!is.null(d_inv_treat_sel)
	)
	masked <- numeric(R)
	masked[V_e] <- cohort_probs_overall[V_e]
	S_V <- sum(masked)
	if (S_V <= 0 || length(V_e) <= 1L) {
		return(matrix(0, nrow = R, ncol = ncol(d_inv_treat_sel)))
	}
	jacobian_e <- matrix(0, nrow = R, ncol = ncol(d_inv_treat_sel))
	for (g in V_e) {
		cons_g <- (S_V - masked[g]) / S_V^2
		jacobian_e[g, ] <- cons_g * d_inv_treat_sel[first_inds[g] + e, ]
		for (g_prime in setdiff(V_e, g)) {
			cons_off <- masked[g_prime] / S_V^2
			jacobian_e[g, ] <- jacobian_e[g, ] -
				cons_off * d_inv_treat_sel[first_inds[g_prime] + e, ]
		}
	}
	jacobian_e
}


# .assemble_joint_cov_var1
#' @title K x K regression-coefficient covariance block (Sigma_1)
#' @description
#' The K-effect generalization of the scalar `att_var_1` formula in
#' `getTeResults2()` / `getTeResultsOLS()`. Given a `p_sel x K` matrix `Psi`
#' whose column `k` is effect `k`'s contrast on the selected
#' treatment-effect support, returns the K x K block
#' `sig_eps_sq / (N*T) * t(Psi) %*% gram_inv %*% Psi` (model-based path) or,
#' under `se_type = "cluster"`, the zero-padded cluster-robust sandwich form
#' `t(Psi_full) %*% sandwich_full %*% Psi_full`. The K = 1 case reproduces the
#' scalar `att_var_1` exactly. Built by `.simultaneous_cis_impl()`; the joint
#' covariance is not persisted on the fit (see #192 Decision Log).
#' @param Psi Numeric matrix, `p_sel x K`. Column `k` is the contrast vector
#'   that picks effect `k` out of the selected treatment-effect support. The
#'   caller is responsible for shaping `Psi` in the right space (theta-space
#'   for FETWFE, beta-space for the OLS family).
#' @param gram_inv Numeric matrix; the Gram inverse restricted to the selected
#'   treatment-effect features (from `getGramInv()`). Used in the non-cluster
#'   path.
#' @param sig_eps_sq Numeric scalar; idiosyncratic error variance.
#' @param N,T Integers; units and time periods.
#' @param se_type Character; `"default"` / `"conservative"` use the model-based
#'   form, `"cluster"` uses the sandwich form.
#' @param sandwich_full Numeric matrix or `NULL`; the full cluster-robust
#'   sandwich on the selected support (required when `se_type = "cluster"`).
#' @param treat_block_mask Logical vector or `NULL`; marks which selected
#'   columns are treatment-effect features (required when
#'   `se_type = "cluster"`), used to zero-pad `Psi` to the sandwich's full
#'   feature space.
#' @return A numeric K x K matrix. Diagonal entries are floored at zero in the
#'   cluster path (issue #84 item 9 / issue #127), mirroring the scalar sites.
#' @keywords internal
#' @noRd
.assemble_joint_cov_var1 <- function(
	Psi,
	gram_inv,
	sig_eps_sq,
	N,
	T,
	se_type = "default",
	sandwich_full = NULL,
	treat_block_mask = NULL
) {
	K <- ncol(Psi)
	if (identical(se_type, "cluster")) {
		stopifnot(!is.null(sandwich_full))
		stopifnot(!is.null(treat_block_mask))
		stopifnot(sum(treat_block_mask) == nrow(Psi))
		Psi_full <- matrix(0, nrow = length(treat_block_mask), ncol = K)
		Psi_full[treat_block_mask, ] <- Psi
		out <- t(Psi_full) %*% sandwich_full %*% Psi_full
		# Floor diagonal entries at zero against floating-point drift; the
		# Liang-Zeger sandwich quadratic form is PSD in exact arithmetic
		# (issue #84 item 9 / issue #127). Off-diagonals can be either sign
		# and are left as-is.
		diag(out) <- pmax(diag(out), 0)
	} else {
		out <- sig_eps_sq * (t(Psi) %*% gram_inv %*% Psi) / (N * T)
	}
	out
}


# .assemble_joint_cov_var2
#' @title K x K cohort-probability covariance block (Sigma_2)
#' @description
#' The K-effect generalization of the scalar `att_var_2` (`getSecondVarTermOLS`
#' / `getSecondVarTermDataApp`). Block `(k, l)` of the K x K output is
#' `T/(N*T) * theta_sel' J_k' Sigma_pi_hat J_l theta_sel`, where `J_k` is the
#' per-effect Jacobian (R x p_sel for FETWFE; R x R for the OLS family) built by
#' `.build_jacobian()`. The K = 1 case reproduces the scalar `att_var_2`.
#'
#' Round-1 B1: the single-global-Jacobian sketch was empirically wrong for
#' `family = "event_study"` (it gave the wrong diagonal against the existing
#' `var_2(e)` scalars). Per-effect Jacobians (each masked to its own valid
#' cohort set) reproduce the existing scalars exactly. Per-family `J_list`
#' construction lives in `.simultaneous_cis_impl()`. Round-2 N3: the single
#' global `Sigma_pi_hat` is correct --- the zero-rows of `J_k` for cohorts not in
#' effect `k`'s valid set zero out the relevant `Sigma_pi_hat` entries
#' automatically, so a per-effect masked `Sigma_pi_hat` is unnecessary.
#' @param J_list Length-K list of per-effect Jacobian matrices (each R x p_sel
#'   for FETWFE, R x R for the OLS family).
#' @param theta_sel Numeric vector; the selected treatment-effect coefficient
#'   vector (theta-space for FETWFE, beta-space for the OLS family). Length
#'   matches `ncol(J_list[[k]])`.
#' @param Sigma_pi_hat Numeric matrix, R x R; the multinomial covariance of the
#'   cohort-count vector (from `.multinomial_cov(cohort_probs_overall[1:R])`).
#' @param N,T Integers; units and time periods.
#' @return A numeric K x K matrix. Diagonal entries are floored at zero (issue
#'   #127), mirroring the scalar sites.
#' @keywords internal
#' @noRd
.assemble_joint_cov_var2 <- function(
	J_list,
	theta_sel,
	Sigma_pi_hat,
	N,
	T
) {
	K <- length(J_list)
	out <- matrix(0, K, K)
	for (k in seq_len(K)) {
		for (l in seq_len(K)) {
			out[k, l] <- T *
				as.numeric(
					t(theta_sel) %*%
						t(J_list[[k]]) %*%
						Sigma_pi_hat %*%
						J_list[[l]] %*%
						theta_sel
				) /
				(N * T)
		}
	}
	# Floor diagonal (parallels the issue #127 floor at the scalar sites).
	diag(out) <- pmax(diag(out), 0)
	out
}


# getCohortATTsFinal
#' @description
#' Computes the **Cohort Average Treatment Effect on the Treated**
#' \deqn{\tau_{\text{ATT},g}\;=\;\frac1{T-g+1}\sum_{t=g}^{T}\tau_{\text{ATT}}(g,t)}
#' for every treated cohort \(g\in\{1,\dots,R\}\).\cr
#' In the fused-parameterisation used by the estimator, each
#' \(\tau_{\text{ATT}}(g,t)\) is a *row* of
#' \((D^{(2)}(\mathcal R))^{-1}\widehat\theta\).
#' Averaging those rows within a cohort produces a
#' **weight vector** \(\psi_g\) such that
#' \(\widehat{\tau}_{\text{ATT},g}=\psi_g^{\!\top}\widehat\theta\).
#' The function
#'
#' * constructs every \(\psi_g\) (via \code{getPsiRFused()})
#' * builds the Gram inverse
#'   \(\bigl((NT)^{-1}X_{\hat{\mathcal S}}^{\top}X_{\hat{\mathcal S}}\bigr)^{-1}\)
#'   for the *selected* treatment-effect columns
#' * returns point estimates, \(\sqrt{\widehat{\operatorname{Var}}}\) and
#'   Wald intervals
#'   \(\widehat{\tau}_{\text{ATT},g}\pm z_{1-\alpha/2}\,
#'     \widehat{\operatorname{SE}}(\widehat{\tau}_{\text{ATT},g})\).
#'
#' @param X_final Numeric matrix; the final design matrix, potentially
#'   transformed by `Omega_sqrt_inv` and the fusion transformation.
#' @param sel_feat_inds Integer vector OR `NULL`. Indices of all features
#'   selected by the penalized regression in the transformed space. Pass
#'   `NULL` for OLS callers (`etwfe()` / `twfeCovs()`) where no penalized
#'   selection occurred --- the unified function treats `NULL` as a sentinel
#'   for "use all features" and dispatches the OLS code path inside
#'   `getGramInv()` and `.assemble_cluster_robust_sandwich()`.
#' @param treat_inds Integer vector; indices in the original (untransformed)
#'   coefficient vector that correspond to the base treatment effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the block of `num_treats` treatment effect parameters.
#' @param sel_treat_inds_shifted Integer vector; indices of selected treatment
#'   effects within the `num_treats` block (shifted to start from 1). For OLS
#'   callers pass `seq_len(num_treats)`.
#' @param c_names Character vector; names of the `R` treated cohorts.
#' @param tes Numeric vector; estimated treatment effects in the original
#'   parameterization for all `num_treats` possible cohort-time combinations.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param R Integer; total number of treated cohorts.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param fused Logical; if `TRUE`, assumes fusion penalization was used and
#'   accumulates the `d_inv_treat_sel` block from the inverse fusion matrix
#'   for downstream variance-from-cohort-probabilities propagation. If
#'   `FALSE`, the unfused (OLS / BETWFE) path is taken: `psi_g` is computed by
#'   `getPsiRUnfused()` and no `d_inv_treat_sel` is returned.
#' @param calc_ses Logical; if `TRUE`, attempts to calculate standard errors.
#'   For OLS callers pass `TRUE` (the gram-matrix inverse is needed
#'   unconditionally for SEs); for bridge callers pass `q < 1`.
#' @param include_selected Logical; if `TRUE`, appends a trailing `selected =
#'   cohort_tes != 0` column to `cohort_te_df`. The column is meaningful only
#'   when a bridge penalty produced the `tes` (FETWFE / BETWFE), so OLS
#'   callers (`etwfe()` / `twfeCovs()`) pass `FALSE` to suppress the column.
#'   New in 1.9.27.
#' @param alpha Numeric scalar; significance level for confidence intervals
#'   (e.g., 0.05 for 95% CIs).
#' @param se_type Character; one of `"default"` (tight Gaussian variance under
#'   (Psi-IF), Theorem (c$'$)), `"conservative"` (Cauchy-Schwarz upper bound
#'   from Theorem (c) for non-(Psi-IF) propensity estimators), or `"cluster"`
#'   (experimental unit-clustered Liang-Zeger sandwich SE on the
#'   OLS-on-selected-support residuals). Default `"default"`. Note that
#'   `se_type` affects only the overall-ATT variance combination in
#'   `getTeResults2()` / `getTeResultsOLS()`; cohort-specific SEs computed
#'   here are governed by `se_type %in% c("default", "conservative")` (the
#'   model-based formula) vs `se_type = "cluster"` (the sandwich path).
#' @param y_final Numeric vector; the (fusion-then-)GLS-transformed response of
#'   length `N*T` (or `N*T + p` if a ridge augmentation was applied upstream).
#'   Required when `se_type = "cluster"` and `calc_ses = TRUE`; ignored
#'   otherwise.
#' @return A list containing:
#'   \item{cohort_te_df}{Data frame (with S3 class `c("catt_df", "data.frame")`)
#'     with columns `cohort`, `estimate`, `se`, `ci_low`, `ci_high`, and a
#'     `p_value` column (two-sided `2 * pnorm(-|estimate / se|)`, `NA` when
#'     `se` is zero or `NA`). When `include_selected = TRUE`, an additional
#'     trailing `selected` logical column is appended (`TRUE` when the
#'     bridge penalty left the cohort's `estimate` nonzero). The
#'     `catt_df` S3 class makes `[[` / `$` / `[` access on the
#'     pre-1.11.0 Title-Case column names `stop()` with a migration message
#'     pointing to the new name.}
#'   \item{cohort_tes}{Named numeric vector of estimated ATTs for each cohort.}
#'   \item{cohort_te_ses}{Named numeric vector of standard errors for cohort ATTs.}
#'   \item{psi_mat}{Matrix used in SE calculation for overall ATT.}
#'   \item{gram_inv}{(Potentially NA) Inverse of the Gram matrix for selected
#'     features, used in SE calculation.}
#'   \item{d_inv_treat_sel}{(If `fused=TRUE`) Relevant block of the inverse
#'     fusion matrix for selected treatment effects.}
#'   \item{calc_ses}{Logical, indicating if SEs were actually calculated.}
#'   \item{sandwich_full}{(`se_type = "cluster"` only.) The full cluster-robust
#'     sandwich variance on the selected support, reusable downstream by
#'     `getTeResults2()` / `getTeResultsOLS()` so the same matrix backs both
#'     cohort and overall ATT SEs.}
#'   \item{treat_block_mask}{(`se_type = "cluster"` only.) Logical vector of
#'     length `length(sel_feat_inds)` (or `ncol(X_final)` when
#'     `sel_feat_inds = NULL`) marking which selected columns correspond to
#'     treatment-effect features; used to zero-pad `psi_g` / `psi_att`
#'     against `sandwich_full`.}
#' @details The function first computes the Gram matrix inverse (`gram_inv`) if
#'   `calc_ses` is `TRUE`. Then, for each cohort `g`, it calculates the average
#'   of the relevant `tes`. If SEs are calculated, it uses `getPsiRFused` or
#'   `getPsiRUnfused` to get a `psi_g` vector, which is then used with
#'   `gram_inv` to find the standard error for that cohort's ATT.
#'
#' The finite-sample variance estimator is
#' \deqn{\widehat{\operatorname{Var}}(\widehat{\tau}_{\text{ATT},g})
#'       =\frac{\sigma^{2}}{NT}\;
#'         \psi_g^{\!\top}\,
#'         \widehat G^{-1}\psi_g ,}
#' where \(\widehat G^{-1}\) is the Gram inverse restricted to the selected
#' treatment-effect features.
#' All matrix slices come from the *inverse fusion matrix*
#' \(D^{(2)}(\mathcal R)^{-1}\) and therefore depend only on the known
#' cohort structure.
#'
#' When `calc_ses = TRUE`, the routine also accumulates a
#' block `d_inv_treat_sel` -
#' the sub-matrix of \(D^{(2)}(\mathcal R)^{-1}\) with
#' **all rows** but **only the selected columns**.
#' This block is later required by \code{getSecondVarTermDataApp()}
#' to propagate sampling noise in the cohort-probability estimates.
#' @seealso \code{getPsiRFused()}, \code{getTeResults2()}
#' @keywords internal
#' @noRd
getCohortATTsFinal <- function(
	X_final,
	sel_feat_inds,
	treat_inds,
	num_treats,
	first_inds,
	sel_treat_inds_shifted,
	c_names,
	tes,
	sig_eps_sq,
	R,
	N,
	T,
	fused,
	calc_ses,
	include_selected,
	alpha = 0.05,
	se_type = "default",
	y_final = NULL
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)

	stopifnot(max(sel_treat_inds_shifted) <= num_treats)
	stopifnot(min(sel_treat_inds_shifted) >= 1)
	stopifnot(length(tes) == num_treats)
	stopifnot(all(!is.na(tes)))

	stopifnot(nrow(X_final) == N * T)
	X_to_pass <- X_final

	# Translate the OLS-caller `sel_feat_inds = NULL` sentinel to the value
	# each downstream helper expects: `getGramInv()` uses NA as its "all
	# features" sentinel, while `.assemble_cluster_robust_sandwich()` uses
	# NULL. See #118 round-1 plan review for the empirical equivalence
	# verification.
	gram_sel_feat <- if (is.null(sel_feat_inds)) NA else sel_feat_inds
	gram_sel_treat <- if (is.null(sel_feat_inds)) NA else sel_treat_inds_shifted

	# Start by getting Gram matrix needed for standard errors
	if (calc_ses) {
		res <- getGramInv(
			N = N,
			T = T,
			X_final = X_to_pass,
			sel_feat_inds = gram_sel_feat,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_treat_inds_shifted = gram_sel_treat,
			calc_ses = calc_ses
		)

		gram_inv <- res$gram_inv
		calc_ses <- res$calc_ses
	} else {
		gram_inv <- NA
	}

	# Cluster-robust sandwich (computed once outside the cohort loop and
	# reused for the overall ATT in getTeResults2() / getTeResultsOLS()).
	sandwich_full <- NULL
	treat_block_mask <- NULL

	if (identical(se_type, "cluster") && calc_ses) {
		stopifnot(!is.null(y_final))
		stopifnot(length(y_final) >= N * T)
		res <- .assemble_cluster_robust_sandwich(
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			treat_inds = treat_inds,
			sel_feat_inds = sel_feat_inds
		)
		sandwich_full <- res$sandwich_full
		treat_block_mask <- res$treat_block_mask
		# Branch-aware mask-length check: in the OLS-caller path the mask
		# spans ncol(X_final); in the bridge path it spans length(sel_feat_inds).
		expected_mask_len <- if (is.null(sel_feat_inds)) {
			ncol(X_final)
		} else {
			length(sel_feat_inds)
		}
		stopifnot(length(treat_block_mask) == expected_mask_len)
		stopifnot(sum(treat_block_mask) == length(sel_treat_inds_shifted))
	}

	## We need the tau sub-matrix of D^{-1} ONLY if we are in fused workflow
	if (fused) {
		# Get the parts of D_inv that have to do with treatment effects
		d_inv_treat <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

		## we will progressively rbind() the selected-column rows of this matrix
		d_inv_treat_sel <- matrix(
			0,
			nrow = 0,
			ncol = length(sel_treat_inds_shifted)
		)
	}

	# First, each cohort
	cohort_tes <- rep(as.numeric(NA), R) # cohort average treatment effect point estimates
	cohort_te_ses <- rep(as.numeric(NA), R) # standard errors

	psi_mat <- matrix(0, length(sel_treat_inds_shifted), R)

	# loop over cohorts
	for (g in 1:R) {
		# Get indices corresponding to marginal treatment effects for rth cohort.
		# Endpoint-preservation form: downstream `getPsiRFused()` takes
		# `first_ind_g, last_ind_g` as separate positional args.
		inds_g <- .cohort_block_inds(g, R, first_inds, num_treats)
		first_ind_g <- inds_g[1]
		last_ind_g <- inds_g[length(inds_g)]

		stopifnot(last_ind_g >= first_ind_g)
		stopifnot(all(first_ind_g:last_ind_g %in% 1:num_treats))

		# average treatment effect for cohort g: simple mean of these treatment
		# effects
		cohort_tes[g] <- mean(tes[first_ind_g:last_ind_g])

		if (calc_ses) {
			# Calculate standard errors

			if (fused) {
				# build psi_g and extract block of D^{-1} for these rows/
				# selected columns
				res_g <- getPsiRFused(
					first_ind_g,
					last_ind_g,
					sel_treat_inds_shifted,
					d_inv_treat
				)

				psi_g <- res_g$psi_g

				stopifnot(
					nrow(res_g$d_inv_treat_sel) ==
						last_ind_g -
							first_ind_g +
							1
				)
				stopifnot(
					ncol(res_g$d_inv_treat_sel) ==
						length(sel_treat_inds_shifted)
				)

				stopifnot(is.matrix(res_g$d_inv_treat_sel))

				d_inv_treat_sel <- rbind(d_inv_treat_sel, res_g$d_inv_treat_sel)

				if (nrow(d_inv_treat_sel) != last_ind_g) {
					err_mes <- paste(
						"nrow(d_inv_treat_sel) == last_ind_g is not TRUE. ",
						"nrow(d_inv_treat_sel): ",
						nrow(d_inv_treat_sel),
						". num_treats: ",
						num_treats,
						". R: ",
						R,
						". first_inds: ",
						paste(first_inds, collapse = ", "),
						". g: ",
						g,
						". first_ind_g: ",
						first_ind_g,
						". last_ind_g: ",
						last_ind_g,
						". nrow(res_g$d_inv_treat_sel):",
						nrow(res_g$d_inv_treat_sel)
					)
					stop(err_mes)
				}

				rm(res_g)
			} else {
				psi_g <- getPsiRUnfused(
					first_ind_g,
					last_ind_g,
					sel_treat_inds_shifted
				)
			}

			stopifnot(length(psi_g) == length(sel_treat_inds_shifted))

			psi_mat[, g] <- psi_g
			# Get standard errors

			if (identical(se_type, "cluster")) {
				# Use length(treat_block_mask), not length(sel_feat_inds):
				# the OLS-caller path passes sel_feat_inds = NULL, and
				# `length(NULL) == 0` would silently yield a zero-length
				# psi_g_full and a broken cluster SE. `treat_block_mask`
				# always has the right length for this zero-padding —
				# `ncol(X_final)` in the OLS path, `length(sel_feat_inds)`
				# in the bridge path. Matches the legacy OLS form (see
				# git history for `getCohortATTsFinalOLS()` pre-#118).
				psi_g_full <- numeric(length(treat_block_mask))
				psi_g_full[treat_block_mask] <- psi_g
				# Issue #84 item 9: floor the cluster-sandwich quadratic
				# form at zero before sqrt(), defending against a
				# rounding-induced negative; see the matching guard in
				# `getTeResults2` (same file) for rationale. Issue #139
				# layers a two-tier (warning / error) diagnostic on top of
				# the floor via `.floor_cluster_quad()`; on well-conditioned
				# data the behavior is unchanged.
				cohort_te_ses[g] <- sqrt(.floor_cluster_quad(
					as.numeric(
						t(psi_g_full) %*%
							sandwich_full %*%
							psi_g_full
					),
					"getCohortATTsFinal/cohort_te_se"
				))
			} else {
				## Variance of the treatment effect for cohort g is
				## sigma^2 /(NT) x t(psi_g) %*% gram_inv %*% psi_g (see paper)
				cohort_te_ses[g] <- sqrt(
					sig_eps_sq *
						as.numeric(
							t(psi_g) %*%
								gram_inv %*%
								psi_g
						) /
						(N * T)
				)
			}
		}
	}

	if (fused & calc_ses) {
		if (nrow(d_inv_treat_sel) != num_treats) {
			err_mes <- paste(
				"nrow(d_inv_treat_sel) == num_treats is not TRUE. ",
				"nrow(d_inv_treat_sel): ",
				nrow(d_inv_treat_sel),
				". num_treats: ",
				num_treats,
				". R: ",
				R,
				". first_inds: ",
				paste(first_inds, collapse = ", "),
				"."
			)
			stop(err_mes)
		}
	}

	stopifnot(length(c_names) == R)
	stopifnot(length(cohort_tes) == R)

	if (calc_ses & all(!is.na(gram_inv))) {
		stopifnot(length(cohort_te_ses) == R)

		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			cohort_te_ses,
			cohort_tes - stats::qnorm(1 - alpha / 2) * cohort_te_ses,
			cohort_tes + stats::qnorm(1 - alpha / 2) * cohort_te_ses,
			.compute_p_values(cohort_tes, cohort_te_ses)
		)

		names(cohort_te_ses) <- c_names
		names(cohort_tes) <- c_names
	} else {
		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			rep(NA, R),
			rep(NA, R),
			rep(NA, R),
			rep(NA_real_, R)
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

	# The trailing `selected` column is meaningful only for bridge callers
	# (FETWFE / BETWFE) where some cohorts may be zeroed out by the penalty.
	# OLS callers (etwfe / twfeCovs) suppress it via include_selected = FALSE.
	# `$<-` appends both the column and its name implicitly, so the
	# colnames() block above lists only the 6 unconditional names.
	if (include_selected) {
		cohort_te_df$selected <- cohort_tes != 0
	}

	# Attach the `catt_df` S3 class so [[ / $ / [ accessors fire the
	# helpful-error layer (R/catt_df_class.R) when users hit the old
	# Title-Case column names from versions <= 1.10.0.
	class(cohort_te_df) <- c("catt_df", "data.frame")

	if (fused) {
		stopifnot(is.matrix(d_inv_treat_sel))
		ret <- list(
			cohort_te_df = cohort_te_df,
			cohort_tes = cohort_tes,
			cohort_te_ses = cohort_te_ses,
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			d_inv_treat_sel = d_inv_treat_sel,
			calc_ses = calc_ses,
			sandwich_full = sandwich_full,
			treat_block_mask = treat_block_mask
		)
	} else {
		ret <- list(
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
	return(ret)
}
