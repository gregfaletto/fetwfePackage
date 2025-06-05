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
#' @return A list containing:
#'   \item{cohort_te_df}{Dataframe with cohort names, estimated ATTs, SEs, and
#'     confidence interval bounds.}
#'   \item{cohort_tes}{Named numeric vector of estimated ATTs for each cohort.}
#'   \item{cohort_te_ses}{Named numeric vector of standard errors for cohort ATTs.}
#'   \item{psi_mat}{Matrix used in SE calculation for overall ATT.}
#'   \item{gram_inv}{(Potentially NA) Inverse of the Gram matrix for selected
#'     features, used in SE calculation.}
#'   \item{d_inv_treat_sel}{(If `fused=TRUE`) Relevant block of the inverse
#'     fusion matrix for selected treatment effects.}
#' @details The function first computes the Gram matrix inverse (`gram_inv`).
#'   Then, for each cohort `r`, it calculates the average
#'   of the relevant `tes`. It uses `getPsiRFused` or
#'   `getPsiRUnfused` to get a `psi_r` vector, which is then used with
#'   `gram_inv` to find the standard error for that cohort's ATT.
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
	alpha = 0.05
) {
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
			sel_treat_inds_shifted = 1:num_treats,
			gram_inv = gram_inv
		)

		stopifnot(length(psi_r) == num_treats)

		psi_mat[, r] <- psi_r
		# Get standard errors

		if (calc_ses) {
			cohort_te_ses[r] <- sqrt(
				sig_eps_sq *
					as.numeric(t(psi_r) %*% gram_inv %*% psi_r) /
					(N * T)
			)
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
			cohort_tes + stats::qnorm(1 - alpha / 2) * cohort_te_ses
		)

		names(cohort_te_ses) <- c_names
		names(cohort_tes) <- c_names
	} else {
		cohort_te_df <- data.frame(
			c_names,
			cohort_tes,
			rep(NA, R),
			rep(NA, R),
			rep(NA, R)
		)
	}

	colnames(cohort_te_df) <- c(
		"Cohort",
		"Estimated TE",
		"SE",
		"ConfIntLow",
		"ConfIntHigh"
	)

	stopifnot(length(tes) == nrow(psi_mat))

	ret <- list(
		cohort_te_df = cohort_te_df,
		cohort_tes = cohort_tes,
		cohort_te_ses = cohort_te_ses,
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		calc_ses = calc_ses
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
#'     `psi_att = psi_mat %*% cohort_probs` and `gram_inv`.
#'   - `att_var_2` (variance from cohort probability estimation) is computed by
#'     calling `getSecondVarTermDataApp`.
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
	indep_probs = FALSE
) {
	att_hat <- as.numeric(cohort_tes %*% cohort_probs)

	if (calc_ses) {
		stopifnot(nrow(psi_mat) == length(tes))
		stopifnot(nrow(psi_mat) <= num_treats)
		stopifnot(ncol(psi_mat) == R)

		# Get ATT standard error
		# first variance term: convergence of theta
		psi_att <- psi_mat %*% cohort_probs

		att_var_1 <- sig_eps_sq *
			as.numeric(t(psi_att) %*% gram_inv %*% psi_att) /
			(N * T)

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
	stopifnot(sum(cohort_probs_overall) < 1 - 10e-6)

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
