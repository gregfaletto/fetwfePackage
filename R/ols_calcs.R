# getCohortATTsFinalOLS
#' @title Calculate Cohort-Specific ATTs and Standard Errors for ETWFE
#' @description Computes the Average Treatment Effect on the Treated (ATT) for
#'   each cohort, along with their standard errors and confidence intervals if
#'   requested and feasible.
#' @param X_final Numeric matrix; the final design matrix, potentially
#'   transformed by `Omega_sqrt_inv` and the fusion transformation.
#' @param treat_inds Integer vector; indices in the original (untransformed)
#'   coefficient vector that correspond to the base treatment effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the block of `num_treats` treatment effect parameters.
#' @param c_names Character vector; names of the `R` treated cohorts.
#' @param tes Numeric vector; estimated treatment effects in the original
#'   parameterization for all `num_treats` possible cohort-time combinations.
#' @param sig_eps_sq Numeric scalar; variance of the idiosyncratic error term.
#' @param R Integer; total number of treated cohorts.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param p Integer; total number of parameters in the model.
#' @param alpha Numeric scalar; significance level for confidence intervals
#'   (e.g., 0.05 for 95% CIs).
#' @param se_type Character; one of `"default"` (Assumption-F1 model-based SE)
#'   or `"cluster"` (experimental unit-clustered Liang-Zeger sandwich SE on the
#'   OLS-on-selected-support residuals). Default `"default"`.
#' @param y_final Numeric vector; GLS-transformed response of length `N*T`
#'   (or `N*T + p` if a ridge augmentation was applied upstream). Required
#'   when `se_type = "cluster"`; ignored otherwise.
#' @return A list containing:
#'   \item{cohort_te_df}{Dataframe with cohort names, estimated ATTs, SEs,
#'     confidence interval bounds, and a `P_value` column (two-sided
#'     `2 * pnorm(-|estimate / se|)`, `NA` when `SE` is zero or `NA`).}
#'   \item{cohort_tes}{Named numeric vector of estimated ATTs for each cohort.}
#'   \item{cohort_te_ses}{Named numeric vector of standard errors for cohort ATTs.}
#'   \item{psi_mat}{Matrix used in SE calculation for overall ATT.}
#'   \item{gram_inv}{(Potentially NA) Inverse of the Gram matrix for selected
#'     features, used in SE calculation.}
#'   \item{d_inv_treat_sel}{(If `fused=TRUE`) Relevant block of the inverse
#'     fusion matrix for selected treatment effects.}
#'   \item{sandwich_full}{(`se_type = "cluster"` only.) The full cluster-robust
#'     sandwich variance on the selected support, reusable downstream by
#'     `getTeResultsOLS()` so the same matrix backs both cohort and overall
#'     ATT SEs.}
#'   \item{treat_block_mask}{(`se_type = "cluster"` only.) Logical vector of
#'     length `ncol(X_final)` marking which selected columns correspond to
#'     treatment-effect features; used to zero-pad `psi_r` / `psi_att` against
#'     `sandwich_full`.}
#' @details The function first computes the Gram matrix inverse (`gram_inv`).
#'   Then, for each cohort `r`, it calculates the average
#'   of the relevant `tes`. It uses `getPsiRFused` or
#'   `getPsiRUnfused` to get a `psi_r` vector, which is then used with
#'   `gram_inv` to find the standard error for that cohort's ATT.
#'   Under `se_type = "cluster"`, the function additionally re-solves OLS on
#'   the selected support (here the full `X_final`, since ETWFE/twfeCovs do
#'   not penalise) and forms a unit-clustered Liang-Zeger sandwich via
#'   `.compute_cluster_robust_sandwich()`; the cohort SE is then the
#'   quadratic form of the zero-padded `psi_r` against this sandwich.
#' @keywords internal
#' @noRd
getCohortATTsFinalOLS <- function(
	X_final,
	treat_inds,
	num_treats,
	first_inds,
	c_names,
	tes,
	sig_eps_sq,
	R,
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

	# Start by getting Gram matrix needed for standard errors
	res <- getGramInv(
		N = N,
		T = T,
		X_final = X_to_pass,
		treat_inds = treat_inds,
		num_treats = num_treats,
		calc_ses = TRUE
	)

	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses

	# Cluster-robust sandwich (computed once outside the cohort loop and
	# reused for the overall ATT in getTeResultsOLS).
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
			sel_feat_inds = NULL
		)
		sandwich_full <- res$sandwich_full
		treat_block_mask <- res$treat_block_mask
	}

	# First, each cohort
	cohort_tes <- rep(as.numeric(NA), R)
	cohort_te_ses <- rep(as.numeric(NA), R)

	psi_mat <- matrix(0, num_treats, R)

	for (r in 1:R) {
		# Get indices corresponding to rth treatment
		first_ind_r <- first_inds[r]
		if (r < R) {
			last_ind_r <- first_inds[r + 1] - 1
		} else {
			last_ind_r <- num_treats
		}

		stopifnot(last_ind_r >= first_ind_r)
		stopifnot(all(first_ind_r:last_ind_r %in% 1:num_treats))

		cohort_tes[r] <- mean(tes[first_ind_r:last_ind_r])

		psi_r <- getPsiRUnfused(
			first_ind_r,
			last_ind_r,
			sel_treat_inds_shifted = 1:num_treats
		)

		stopifnot(length(psi_r) == num_treats)

		psi_mat[, r] <- psi_r
		# Get standard errors

		if (calc_ses) {
			if (identical(se_type, "cluster")) {
				psi_r_full <- numeric(length(treat_block_mask))
				psi_r_full[treat_block_mask] <- psi_r
				cohort_te_ses[r] <- sqrt(as.numeric(
					t(psi_r_full) %*%
						sandwich_full %*%
						psi_r_full
				))
			} else {
				cohort_te_ses[r] <- sqrt(
					sig_eps_sq *
						as.numeric(t(psi_r) %*% gram_inv %*% psi_r) /
						(N * T)
				)
			}
		}
	}

	stopifnot(length(c_names) == R)
	stopifnot(length(cohort_tes) == R)

	if (all(!is.na(gram_inv))) {
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
		"Cohort",
		"Estimated TE",
		"SE",
		"ConfIntLow",
		"ConfIntHigh",
		"P_value"
	)

	stopifnot(length(tes) == nrow(psi_mat))

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
	return(ret)
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
#'   estimated probabilities of belonging to cohort `r` conditional on being
#'   treated. Length `R`. Sums to 1.
#' @param psi_mat Numeric matrix; matrix where column `r` is `psi_r` (from
#'   `getCohortATTsFinal`). Dimensions: `length(sel_treat_inds_shifted)` x `R`.
#' @param gram_inv Numeric matrix; inverse of the Gram matrix for selected
#'   treatment effect features.
#' @param tes Numeric vector; all `num_treats` estimated treatment effects
#'   (original parameterization).
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort P(W=r). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort.
#' @param calc_ses Logical; if `TRUE`, calculate standard errors.
#' @param indep_probs Logical; if `TRUE`, assumes `cohort_probs` (and
#'   `cohort_probs_overall`) were estimated from an independent sample, leading
#'   to a different SE formula (sum of variances) compared to when they are
#'   estimated from the same sample (conservative SE including a covariance term).
#' @param se_type Character; one of `"default"` (model-based SE) or `"cluster"`
#'   (experimental unit-clustered Liang-Zeger sandwich SE). Default
#'   `"default"`.
#' @param sandwich_full Numeric matrix (or `NULL`); the full cluster-robust
#'   sandwich variance on the selected support, as returned by
#'   `getCohortATTsFinalOLS()`. Required when `se_type = "cluster"` and
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
#' @details The overall ATT (`att_hat`) is `cohort_tes %*% cohort_probs`.
#'   If `calc_ses` is `TRUE`:
#'   - `att_var_1` (variance from `theta_hat` estimation) is computed using
#'     `psi_att = psi_mat %*% cohort_probs` and `gram_inv` (or, under
#'     `se_type = "cluster"`, using the zero-padded `psi_att_full` against
#'     `sandwich_full`).
#'   - `att_var_2` (variance from cohort probability estimation) is computed by
#'     calling `getSecondVarTermDataApp`. It is unchanged under
#'     `se_type = "cluster"`.
#'   - `att_te_se` is `sqrt(att_var_1 + att_var_2)` if `indep_probs` is `TRUE`,
#'     otherwise it's a conservative SE: `sqrt(att_var_1 + att_var_2 + 2*sqrt(att_var_1 * att_var_2))`.
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
	first_inds,
	calc_ses,
	indep_probs = FALSE,
	se_type = "default",
	sandwich_full = NULL,
	treat_block_mask = NULL
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

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

			att_var_1 <- as.numeric(
				t(psi_att_full) %*% sandwich_full %*% psi_att_full
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
			first_inds = first_inds,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R
		)

		if (indep_probs) {
			att_te_se <- sqrt(att_var_1 + att_var_2)
		} else {
			att_te_se <- sqrt(
				att_var_1 +
					att_var_2 +
					2 *
						sqrt(
							att_var_1 * att_var_2
						)
			)
		}

		att_te_se_no_prob <- sqrt(att_var_1)
	} else {
		att_te_se <- NA
		att_te_se_no_prob <- NA
	}

	return(list(
		att_hat = att_hat,
		att_te_se = att_te_se,
		att_te_se_no_prob = att_te_se_no_prob
	))
}

# getSecondVarTermOLS
#' @title Calculate Second Variance Term for ATT Standard Error for ETWFE
#' @description Computes the second component of the variance for the Average
#'   Treatment Effect on the Treated (ATT). This component accounts for the
#'   variability due to the estimation of cohort membership probabilities.
#' @param psi_mat Numeric matrix; a matrix where each column `r` is the `psi_r`
#'   vector used in calculating the ATT for cohort `r`. Dimensions:
#'   `length(sel_treat_inds_shifted)` x `R`.
#' @param tes Numeric vector; the estimated treatment effects for all
#'   `num_treats` possible cohort-time combinations.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort (P(W=r)). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the `num_treats` block.
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
	first_inds,
	num_treats,
	N,
	T,
	R
) {
	stopifnot(length(tes) == nrow(psi_mat))
	# Get Sigma_pi_hat, the (sample-estimated) covariance matrix for the
	# sample proportions (derived from the multinomial distribution)
	Sigma_pi_hat <- -outer(
		cohort_probs_overall[1:(R)],
		cohort_probs_overall[1:(R)]
	)
	diag(Sigma_pi_hat) <- cohort_probs_overall[1:R] *
		(1 - cohort_probs_overall[1:R])

	stopifnot(nrow(Sigma_pi_hat) == R)
	stopifnot(ncol(Sigma_pi_hat) == R)

	# Gather a list of the indices corresponding to the treatment coefficients
	# for each cohort. (Will be used to construct the Jacobian matrix.)
	sel_inds <- list()

	for (r in 1:R) {
		first_ind_r <- first_inds[r]

		if (r < R) {
			last_ind_r <- first_inds[r + 1] - 1
		} else {
			last_ind_r <- num_treats
		}
		stopifnot(last_ind_r >= first_ind_r)
		sel_inds[[r]] <- first_ind_r:last_ind_r
		if (r > 1) {
			stopifnot(min(sel_inds[[r]]) > max(sel_inds[[r - 1]]))
			stopifnot(length(sel_inds[[r]]) <= length(sel_inds[[r - 1]]))
		}
	}

	stopifnot(all.equal(unlist(sel_inds), 1:num_treats))

	# Construct Jacobian matrix corresponding to the mapping from the
	# individual cohort probabilities to the proportions calculated for the ATT.
	# See Proof of Theorem 6.1 for form.
	jacobian_mat <- matrix(
		as.numeric(NA),
		nrow = R,
		ncol = R
	)

	stopifnot(length(cohort_probs_overall) == R)
	stopifnot(sum(cohort_probs_overall) < 1 - 1e-6)

	for (r in 1:R) {
		# All terms in rth column have the same value except the diagonal term
		col_r_val <- -cohort_probs_overall[r] / sum(cohort_probs_overall)^2

		jacobian_mat[, r] <- rep(col_r_val, R)

		# Diagonal term
		cons_r <- (sum(cohort_probs_overall) -
			cohort_probs_overall[r]) /
			sum(cohort_probs_overall)^2

		jacobian_mat[r, r] <- cons_r
	}

	stopifnot(all(!is.na(jacobian_mat)))

	stopifnot(nrow(psi_mat) <= num_treats)
	stopifnot(ncol(psi_mat) == R)

	stopifnot(length(tes) == nrow(psi_mat))

	# Finally, calculate variance term for ATT from Sigma_pi_hat, Jacobian
	# matrix, psi_mat, and beta_hat. See proof of Theorem D.2 for details.

	att_var_2 <- T *
		as.numeric(
			t(tes) %*%
				psi_mat %*%
				t(jacobian_mat) %*%
				Sigma_pi_hat %*%
				jacobian_mat %*%
				t(psi_mat) %*%
				tes
		) /
		(N * T)

	return(att_var_2)
}


# getPsiRUnfused
#' @title Calculate Psi Vector for Cohort ATT (Unfused Case)
#' @description Computes the `psi_r` vector for a specific cohort `r` when no
#'   fusion penalization is applied to the treatment effects (or when calculating
#'   SEs as if it were an OLS on selected variables). This vector is used in
#'   standard error calculations for the cohort's Average Treatment Effect on
#'   the Treated (ATT). Specifically, `psi_r` places the constant
#'   `1 / k_full` (where `k_full = last_ind_r - first_ind_r + 1` is the
#'   cohort's full treatment-block size) at each position corresponding to a
#'   selected treatment effect in cohort `r`, and zero elsewhere — so that
#'   `t(psi_r) %*% theta_hat_treat_sel == cohort_tes[r]`.
#' @param first_ind_r Integer; the index of the first treatment effect for
#'   cohort `r` within the `num_treats` block of treatment effects.
#' @param last_ind_r Integer; the index of the last treatment effect for
#'   cohort `r` within the `num_treats` block.
#' @param sel_treat_inds_shifted Integer vector; indices of all selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @return A numeric vector `psi_r` of length
#'   `length(sel_treat_inds_shifted)`. Positions corresponding to selected
#'   treatment effects in cohort `r` carry `1 / k_full`; other positions are
#'   zero. Returns a zero vector if no coefficients in cohort `r`'s block are
#'   selected.
#' @details The function identifies which of `sel_treat_inds_shifted` fall
#'   within the cohort-r block `[first_ind_r, last_ind_r]`. For those
#'   identified indices, `psi_r` is set to `1 / k_full`, where
#'   `k_full = last_ind_r - first_ind_r + 1` is the cohort's full
#'   treatment-block size. This normalization (rather than dividing by the
#'   *selected* count `k_sel`) is what makes
#'   `t(psi_r) %*% theta_hat_treat_sel == cohort_tes[r]`, since the
#'   cohort point estimate `cohort_tes[r] = mean(tes[first_ind_r:last_ind_r])`
#'   averages over the full block (unselected entries are exact zeros
#'   post-bridge). If no coefficients in cohort r's block are selected,
#'   `psi_r` is a zero vector. For ETWFE / `twfeCovs()` callers
#'   (which pass `sel_treat_inds_shifted = 1:num_treats`, so `k_sel = k_full`),
#'   this is bit-identical to the previous divide-by-sum normalization.
#' @keywords internal
#' @noRd
getPsiRUnfused <- function(
	first_ind_r,
	last_ind_r,
	sel_treat_inds_shifted
) {
	k_full <- last_ind_r - first_ind_r + 1
	stopifnot(k_full >= 1)

	which_inds_ir <- sel_treat_inds_shifted %in% (first_ind_r:last_ind_r)

	psi_r <- rep(0, length(sel_treat_inds_shifted))

	if (sum(which_inds_ir) > 0) {
		inds_r <- which(which_inds_ir)

		stopifnot(is.integer(inds_r) | is.numeric(inds_r))
		stopifnot(identical(inds_r, as.integer(round(inds_r))))
		stopifnot(length(inds_r) >= 1)
		stopifnot(length(inds_r) == length(unique(inds_r)))
		stopifnot(length(inds_r) <= length(sel_treat_inds_shifted))
		stopifnot(all(inds_r %in% 1:length(sel_treat_inds_shifted)))

		psi_r[inds_r] <- 1 / k_full
	}

	return(psi_r)
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

#' @title Assemble the cluster-robust sandwich for the OLS-selected support
#' @description
#' Wraps the 4-step assembly ritual used by `getCohortATTsFinalOLS()`,
#' `getCohortATTsFinal()`, and the two `event_study` dispatchers
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
