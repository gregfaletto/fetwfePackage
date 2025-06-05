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
#'       \sum_{r=1}^{R}\widehat{\tau}_{\text{ATT},r}\,
#'                     \widehat{\pi}_{r\mid\tau},}
#' a weighted mean of cohort-specific ATTs with weights
#' \(\widehat{\pi}_{r\mid\tau}=N_r/N_\tau\).
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
#'   estimated probabilities of belonging to cohort `r` conditional on being
#'   treated. Length `R`. Sums to 1.
#' @param psi_mat Numeric matrix; matrix where column `r` is `psi_r` (from
#'   `getCohortATTsFinal`). Dimensions: `length(sel_treat_inds_shifted)` x `R`.
#' @param gram_inv Numeric matrix; inverse of the Gram matrix for selected
#'   treatment effect features.
#' @param sel_treat_inds_shifted Integer vector; indices of selected treatment
#'   effects within the `num_treats` block (shifted to start from 1).
#' @param tes Numeric vector; all `num_treats` estimated treatment effects
#'   (original parameterization).
#' @param d_inv_treat_sel Numeric matrix; block of the inverse fusion matrix for
#'   selected treatment effects.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort P(W=r). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in
#'   transformed (fused) space for selected treatment effects.
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
#'
#' `indep_probs = TRUE` implements the independent-sample
#' variance (`att_var_1 + att_var_2`);
#' `indep_probs = FALSE` returns the conservative bound
#' `att_var_1 + att_var_2 + 2sqrt(att_var_1*att_var_2)`.
#'
#' All matrices required for Term 2 are produced upstream:
#' * `psi_mat` from [getCohortATTsFinal()]
#' * `d_inv_treat_sel` from the same routine
#' * the Jacobian \(J\) is built in
#'   [getSecondVarTermDataApp()] using `d_inv_treat_sel`
#' @inheritParams getTeResults2
#' @seealso [getSecondVarTermDataApp()]
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
	tes,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	calc_ses,
	indep_probs = FALSE
) {
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

		att_var_1 <- sig_eps_sq *
			as.numeric(t(psi_att) %*% gram_inv %*% psi_att) /
			(N * T)

		stopifnot(!is.na(att_var_1))

		# Second variance term: convergence of cohort membership probabilities
		att_var_2 <- getSecondVarTermDataApp(
			# cohort_probs = cohort_probs,
			psi_mat = psi_mat,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			tes = tes,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = cohort_probs_overall,
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel,
			num_treats = num_treats,
			N = N,
			T = T,
			R = R,
			fused = TRUE
		)

		stopifnot(!is.na(att_var_2))

		# Combine the two variance terms
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

# checkFetwfeInputs
#' @title Check Inputs for the main `fetwfe` function
#' @description Validates the inputs provided to the main `fetwfe` function,
#'   ensuring they meet type, dimension, and content requirements. Stops
#'   execution with an error message if any check fails.
#' @param pdata Dataframe; the panel data set.
#' @param time_var Character; name of the time variable column.
#' @param unit_var Character; name of the unit variable column.
#' @param treatment Character; name of the treatment indicator column.
#' @param response Character; name of the response variable column.
#' @param covs Character vector; names of covariate columns. Default `c()`.
#' @param indep_counts Integer vector or NA; counts for independent cohort data.
#'   Default `NA`.
#' @param sig_eps_sq Numeric or NA; variance of idiosyncratic noise. Default `NA`.
#' @param sig_eps_c_sq Numeric or NA; variance of unit-level random effects.
#'   Default `NA`.
#' @param lambda.max Numeric or NA; maximum lambda for `gBridge`. Default `NA`.
#' @param lambda.min Numeric or NA; minimum lambda for `gBridge`. Default `NA`.
#' @param nlambda Integer; number of lambdas for `gBridge`. Default `100`.
#' @param q Numeric; Lq penalty exponent for `gBridge`. Default `0.5`.
#' @param verbose Logical; if TRUE, print progress. Default `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals. Default `0.05`.
#' @param add_ridge Logical; if TRUE, add small ridge penalty. Default `FALSE`.
#' @return Logical `indep_count_data_available`, which is `TRUE` if valid
#'   `indep_counts` were provided, `FALSE` otherwise.
#' @details This function performs a series of `stopifnot` checks on each
#'   parameter. For example:
#'   - `pdata` must be a dataframe with at least 4 rows.
#'   - `time_var`, `unit_var`, `treatment`, `response` must be single characters,
#'     present in `pdata`, and the corresponding columns must have the correct
#'     type (e.g., integer for time, character for unit, 0/1 integer for treatment).
#'   - `covs` if provided, must be characters, present in `pdata`, and columns
#'     must be numeric, integer, or factor.
#'   - `indep_counts` if provided, must be positive integers.
#'   - `sig_eps_sq`, `sig_eps_c_sq` if provided, must be non-negative numerics.
#'   - `lambda.max`, `lambda.min` if provided, must be valid numerics
#'     (`lambda.max > lambda.min >= 0`).
#'   - `q` must be in `(0, 2]`.
#'   - `alpha` must be in `(0, 1)`.
#'   Issues a warning if `alpha > 0.5`.
#' @keywords internal
#' @noRd
checkFetwfeInputs <- function(
	pdata,
	time_var,
	unit_var,
	treatment,
	response,
	covs = c(),
	indep_counts = NA,
	sig_eps_sq = NA,
	sig_eps_c_sq = NA,
	lambda.max = NA,
	lambda.min = NA,
	nlambda = 100,
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE
) {
	res <- checkEtwfeInputs(
		pdata = pdata,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		response = response,
		covs = covs,
		indep_counts = indep_counts,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge
	)

	indep_count_data_available <- res$indep_count_data_available
	pdata <- res$pdata

	rm(res)

	if (any(!is.na(lambda.max))) {
		stopifnot(is.numeric(lambda.max) | is.integer(lambda.max))
		stopifnot(length(lambda.max) == 1)
		stopifnot(lambda.max > 0)
	}

	if (any(!is.na(lambda.min))) {
		stopifnot(is.numeric(lambda.min) | is.integer(lambda.min))
		stopifnot(length(lambda.min) == 1)
		stopifnot(lambda.min >= 0)
		if (any(!is.na(lambda.max))) {
			stopifnot(lambda.max > lambda.min)
		}
	}

	stopifnot(is.numeric(q) | is.integer(q))
	stopifnot(length(q) == 1)
	stopifnot(q > 0)
	stopifnot(q <= 2)

	return(list(
		pdata = pdata,
		indep_count_data_available = indep_count_data_available
	))
}

#' Generate the full \eqn{D^{-1}} transformation matrix.
#'
#' @param first_inds A vector of indices corresponding to the first treatment effect for each
#' treated cohort.
#' @param T Total number of time periods.
#' @param R Number of treated cohorts.
#' @param d Number of covariates (time-invariant).
#' @param num_treats Total number of base treatment effect parameters.
#'
#' @return A matrix of dimension p x p where p = R + (T - 1) + d + d*R + d*(T - 1) + num_treats +
#' d*num_treats.
#' @noRd
genFullInvFusionTransformMat <- function(first_inds, T, R, d, num_treats) {
	# Load required package for block diagonal concatenation.
	if (!requireNamespace("Matrix", quietly = TRUE)) {
		stop("The 'Matrix' package is required but not installed.")
	}

	# Block 1: Cohort fixed effects block, size R x R.
	block1 <- genBackwardsInvFusionTransformMat(R)

	# Block 2: Time fixed effects block, size (T - 1) x (T - 1).
	block2 <- genBackwardsInvFusionTransformMat(T - 1)

	# Block 3: Covariate main effects, identity of dimension d.
	block3 <- if (d > 0) diag(d) else NULL

	# Block 4: Cohort-X interactions: I_d \otimes genBackwardsInvFusionTransformMat(R)
	block4 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(R)) else NULL

	# Block 5: Time-X interactions: I_d \otimes genBackwardsInvFusionTransformMat(T - 1)
	block5 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(T - 1)) else NULL

	# Block 6: Base treatment effects: genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	block6 <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

	# Block 7: Treatment-X interactions: I_d \otimes genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
	block7 <- if (d > 0)
		kronecker(
			diag(d),
			genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
		) else NULL

	# Combine blocks into a block-diagonal matrix.
	# Use Matrix::bdiag which returns a sparse matrix; convert to dense if needed.
	blocks <- list(block1, block2)
	if (!is.null(block3)) blocks <- c(blocks, list(block3))
	if (!is.null(block4)) blocks <- c(blocks, list(block4))
	if (!is.null(block5)) blocks <- c(blocks, list(block5))
	blocks <- c(blocks, list(block6))
	if (!is.null(block7)) blocks <- c(blocks, list(block7))

	full_D_inv <- as.matrix(Matrix::bdiag(blocks))

	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(nrow(full_D_inv) == p)
	stopifnot(ncol(full_D_inv) == p)

	return(full_D_inv)
}

# getPsiRFused
#' @title Calculate Psi Vector and D-inverse Block for Cohort ATT (Fused Case)
#' @description Computes the `psi_r` vector and the relevant block of the
#'   inverse fusion transformation matrix (`d_inv_treat_sel`) for a specific
#'   cohort `r`. These are used in standard error calculations for the cohort's
#'   Average Treatment Effect on the Treated (ATT) when fusion penalization
#'   has been applied.
#'
#' For a treated cohort \(r\) let
#' \(M_r=T-r+1\) be the number of post-treatment periods.
#' The CATT is the *simple average*
#' \(M_r^{-1}\sum_{t=r}^{T}\tau_{\text{ATT}}(r,t)\).
#' In the fused basis each \(\tau_{\text{ATT}}(r,t)\) is a row of
#' \(D^{(2)}(\mathcal R)^{-1}\widehat\theta\),
#' so the averaging weights are the **row-means** of the corresponding block
#' of the inverse fusion matrix.
#' This helper:
#'
#' * extracts those rows (`first_ind_r:last_ind_r`) from
#'   `d_inv_treat`,
#' * averages them column-wise to form \(\psi_r\),
#' * returns the block itself (`d_inv_treat_sel`) because later routines need
#'   it for probability-variance propagation.
#' @param first_ind_r Integer; the index of the first treatment effect parameter
#'   for cohort `r` within the original `num_treats` block (1-based).
#' @param last_ind_r Integer; the index of the last treatment effect parameter
#'   for cohort `r` within the original `num_treats` block (1-based).
#' @param sel_treat_inds_shifted Integer vector; indices (1-based) of the
#'   treatment effects that were selected by the model, relative to the start
#'   of the `num_treats` block. E.g., if original indices 5, 7 were selected
#'   from a block starting at index 1, this would be c(5, 7).
#' @param d_inv_treat Numeric matrix; the full inverse two-way fusion
#'   transformation matrix for all `num_treats` treatment effects. Dimensions:
#'   `num_treats` x `num_treats`.
#' @return A list containing:
#'   \item{psi_r}{Numeric vector. It's the column means of the sub-matrix of
#'     `d_inv_treat` corresponding to rows `first_ind_r:last_ind_r` and columns
#'     specified by `sel_treat_inds_shifted`. If `first_ind_r == last_ind_r`,
#'     it's just that specific row of `d_inv_treat` (subsetted by selected columns).}
#'   \item{d_inv_treat_sel}{Numeric matrix. The sub-matrix of `d_inv_treat`
#'     with rows `first_ind_r:last_ind_r` and columns corresponding to
#'     `sel_treat_inds_shifted`.}
#' @details `psi_r` effectively averages the rows of `d_inv_treat` (that correspond
#'   to cohort `r`'s treatment effects) for the columns that were actually
#'   selected by the model. `d_inv_treat_sel` is this specific block of the
#'   `d_inv_treat` matrix.
#'
#' * If only one transformed treatment coefficient was selected,
#'   both the mean and the returned block are forced to the correct
#'   1 x 1 or 1 x *`k`* shape so that higher-level code can
#'   `rbind()` the blocks without special cases.
#' @inheritParams getPsiRFused
#' @seealso [getCohortATTsFinal()]
#' @keywords internal
#' @noRd
getPsiRFused <- function(
	first_ind_r,
	last_ind_r,
	sel_treat_inds_shifted,
	d_inv_treat
) {
	stopifnot(length(sel_treat_inds_shifted) >= 0)
	stopifnot(last_ind_r >= first_ind_r)
	# Get psi vector: the part of D inverse that we need to look at is the
	# block corresponding to the treatment effect estimates, which is the
	# num_treats x num_treats matrix yielded by
	# genInvTwoWayFusionTransformMat(num_treats, first_inds).

	# Correct rows of matrix

	## psi_r := column-wise mean of those rows  (weights for average treatment
	## effect for cohort r)
	##
	## * If |S| > 1, result is a vector length |S|
	## * If |S| == 1, treat the scalar mean as length-1 vector

	if (last_ind_r > first_ind_r) {
		if (length(sel_treat_inds_shifted) > 1) {
			psi_r <- colMeans(d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			])
		} else {
			psi_r <- mean(d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			])
		}
		if (length(sel_treat_inds_shifted) == 1) {
			# Need to coerce this object to be a matrix with one column so it
			# works smoothly with rbind() later
			d_inv_treat_sel <- matrix(
				d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted],
				ncol = 1
			)
		} else {
			d_inv_treat_sel <- d_inv_treat[
				first_ind_r:last_ind_r,
				sel_treat_inds_shifted
			]
		}
	} else {
		psi_r <- d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted]
		# Since first_ind_r and last_ind_r are the same, need to coerce this
		# object to be a matrix with one row so that it works smoothly with
		# rbind() later

		## Block of D^{-1} used later for probability-variance term
		d_inv_treat_sel <- matrix(
			d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted],
			nrow = 1
		)
	}

	stopifnot(is.matrix(d_inv_treat_sel))

	return(list(psi_r = psi_r, d_inv_treat_sel = d_inv_treat_sel))
}


# getBetaBIC
#' @title Select Optimal Coefficients using BIC from gBridge Fit
#' @description From a `gBridge` fit object (which contains solutions for a
#'   path of lambda penalties), this function selects the optimal set of
#'   coefficients based on the Bayesian Information Criterion (BIC). It also
#'   returns the chosen lambda index and the size of the selected model.
#'   Coefficients are returned on their original scale.
#' @param fit A `gBridge` fit object, typically the output from `grpreg::gBridge()`.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the total number of predictor variables (excluding intercept)
#'   in the model matrix `X_mod`.
#' @param X_mod Numeric matrix; the design matrix (potentially transformed for
#'   FETWFE, and **not** yet GLStransformed or scaled/centered by `my_scale`) that was used to generate `y`.
#'   It's used here to calculate SSE on the original scale of `y`.
#' @param y Numeric vector; the original response variable (before GLS transform and centering)
#'   used to fit the model. Length `N*T`.
#' @param scale_center Numeric vector; the centering values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @param scale_scale Numeric vector; the scaling values used to scale `X_mod`
#'   before fitting `gBridge`. Length `p`.
#' @return A list containing:
#'   \item{theta_hat}{Numeric vector of length `p+1`. The selected coefficients
#'     (including intercept at `theta_hat[1]`) on their original data scale.}
#'   \item{lambda_star_ind}{Integer; the index of the lambda value in `fit$lambda`
#'     that resulted in the best BIC.}
#'   \item{lambda_star_model_size}{Integer; the number of non-zero coefficients
#'     (excluding intercept) in the selected model.}
#' @details The function iterates through each lambda in `fit$lambda`. For each:
#'   1. It extracts the intercept (`eta_s`) and slopes (`beta_s`) on the scaled data.
#'   2. It converts these coefficients back to the original data scale using
#'      `scale_center` and `scale_scale`.
#'   3. It calculates the Sum of Squared Errors (SSE) using `sse_bridge()` with
#'      the original-scale coefficients, original `y`, and `X_mod`.
#'   4. It computes the BIC value: `N*T*log(SSE/(N*T)) + s*log(N*T)`, where `s`
#'      is the number of non-zero coefficients (including intercept).
#'   The set of coefficients corresponding to the minimum BIC is chosen. If multiple
#'   lambdas yield the same minimum BIC, the one resulting in the smallest model
#'   size (fewest non-zero coefficients) is selected.
#'   The final returned `theta_hat` also has its slopes and intercept adjusted back to the original scale.
#' @keywords internal
#' @noRd
getBetaBIC <- function(fit, N, T, p, X_mod, y, scale_center, scale_scale) {
	stopifnot(length(y) == N * T)
	n_lambda <- ncol(fit$beta)
	BICs <- rep(as.numeric(NA), n_lambda)
	model_sizes <- rep(as.integer(NA), n_lambda)

	stopifnot(nrow(fit$beta) == p + 1)

	for (k in 1:n_lambda) {
		## --- extract coefficients on the scaled data -------------
		eta_s <- fit$beta[1, k] # intercept (scaled space)
		beta_s <- fit$beta[2:(p + 1), k] # slopes    (scaled space)

		## --- convert to original scale ---------------------------
		beta_hat_k <- beta_s / scale_scale
		eta_hat_k <- eta_s - sum(scale_center * beta_hat_k)

		# Residual sum of squares
		mse_hat <- sse_bridge(
			eta_hat_k,
			beta_hat_k,
			y = y,
			X_mod = X_mod,
			N = N,
			T = T
		)
		# Number of fitted coefficients
		s <- sum(fit$beta[, k] != 0)
		model_sizes[k] <- s

		stopifnot(is.na(BICs[k]))
		BICs[k] <- N * T * log(mse_hat) + s * log(N * T)
	}

	lambda_star_ind <- which(BICs == min(BICs))
	if (length(lambda_star_ind) == 1) {
		lambda_star_final_ind <- lambda_star_ind
		theta_hat <- fit$beta[, lambda_star_final_ind]
	} else {
		# Choose smallest model size among models with equal BIC
		model_sizes_star <- model_sizes[lambda_star_ind]
		min_model_size_ind <- which(model_sizes_star == min(model_sizes_star))
		lambda_star_final_ind <- lambda_star_ind[min_model_size_ind][1]
		stopifnot(length(lambda_star_final_ind) == 1)
		theta_hat <- fit$beta[, lambda_star_final_ind]
	}
	stopifnot(length(lambda_star_final_ind) == 1)
	stopifnot(length(theta_hat) == p + 1)
	stopifnot(all(!is.na(theta_hat)))

	#
	# Rescale coefficients back to the original scale.
	# The coefficient vector theta_hat is of length (p+1), with the first entry
	# as the intercept.
	# For predictors: original coefficient = beta_scaled / scale_j.
	# The intercept is adjusted as: intercept_original =
	# intercept_scaled - sum(center_j * (beta_scaled/scale_j)).
	#
	adjusted_theta_hat <- theta_hat
	if (length(scale_scale) != p) {
		stop("Length of scale_scale does not match number of predictors (p).")
	}
	for (j in 2:(p + 1)) {
		adjusted_theta_hat[j] <- theta_hat[j] / scale_scale[j - 1]
	}
	adjusted_theta_hat[1] <- theta_hat[1] -
		sum(scale_center * (theta_hat[2:(p + 1)] / scale_scale))

	return(list(
		theta_hat = adjusted_theta_hat,
		lambda_star_ind = lambda_star_final_ind,
		lambda_star_model_size = model_sizes[lambda_star_final_ind]
	))
}


#' Core Estimation Logic for Fused Extended Two-Way Fixed Effects
#'
#' @description
#' This function implements the core estimation steps of the FETWFE methodology.
#' It takes a pre-processed design matrix and response, applies transformations
#' for fusion penalties, handles variance components, performs bridge regression,
#' selects the optimal penalty via BIC, and calculates treatment effects and
#' their standard errors.
#'
#' @param X_ints The design matrix with all fixed effects, covariates, treatment
#'   dummies, and their interactions, as produced by `prepXints`.
#' @param y The centered response vector, as produced by `prepXints`.
#' @param in_sample_counts An integer vector named with cohort identifiers
#'   (including "Never_treated"), indicating the number of units in each cohort
#'   within the data used for estimation.
#' @param N The number of unique units.
#' @param T The number of unique time periods.
#' @param d The number of covariates.
#' @param p The total number of columns in `X_ints` (total parameters).
#' @param num_treats The total number of unique treatment effect parameters.
#' @param first_inds A numeric vector indicating the starting column index for
#'   each cohort's first treatment effect within the treatment effect block.
#' @param indep_counts (Optional) An integer vector of counts for how many units
#'   appear in the untreated cohort plus each of the other `R` cohorts, derived
#'   from an independent dataset. Used for asymptotically exact standard errors for
#'   the ATT. Default is `NA`.
#' @param sig_eps_sq (Optional) Numeric; the known variance of the observation-level
#'   IID noise. If `NA`, it will be estimated. Default is `NA`.
#' @param sig_eps_c_sq (Optional) Numeric; the known variance of the unit-level IID
#'   noise (random effects). If `NA`, it will be estimated. Default is `NA`.
#' @param lambda.max (Optional) Numeric; the maximum `lambda` penalty parameter for
#'   the bridge regression grid search. If `NA`, `grpreg` selects it. Default is `NA`.
#' @param lambda.min (Optional) Numeric; the minimum `lambda` penalty parameter.
#'   If `NA`, `grpreg` selects it. Default is `NA`.
#' @param nlambda (Optional) Integer; the number of `lambda` values in the grid.
#'   Default is 100.
#' @param q (Optional) Numeric; the power of the Lq penalty for fusion regularization
#'   (0 < q <= 2). `q=0.5` is default, `q=1` is lasso, `q=2` is ridge.
#'   Default is 0.5.
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals (e.g., 0.05 for
#'   95% CIs). Default is 0.05.
#' @param add_ridge (Optional) Logical; if `TRUE`, adds a small L2 penalty to
#'   the untransformed coefficients to stabilize estimation. Default is `FALSE`.
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters.
#'   \item **Coordinate Transformation:** Calls `transformXintImproved` to transform
#'     `X_ints` into `X_mod`. This transformation allows a standard bridge
#'     regression penalty on `X_mod` to achieve the desired fusion penalties
#'     on the original coefficients.
#'   \item **Variance Component Handling:**
#'     \itemize{
#'       \item If `sig_eps_sq` or `sig_eps_c_sq` are `NA`, `estOmegaSqrtInv` is
#'         called to estimate them from the data using a fixed-effects ridge
#'         regression.
#'       \item Constructs the covariance matrix `Omega` and its inverse square
#'         root `Omega_sqrt_inv`.
#'       \item Pre-multiplies `y` and `X_mod` by `sqrt(sig_eps_sq) * Omega_sqrt_inv`
#'         (via Kronecker product) to obtain `y_final` and `X_final`, effectively
#'         performing a GLS transformation.
#'     }
#'   \item **Optional Ridge Penalty:** If `add_ridge` is `TRUE`, `X_final_scaled`
#'     (scaled version of `X_final`) and `y_final` are augmented to add an L2
#'     penalty on the *original* (untransformed) coefficient scale. This involves
#'     using `genFullInvFusionTransformMat` to get the inverse of the overall
#'     fusion transformation matrix.
#'   \item **Cohort Probabilities:** Calculates cohort membership probabilities
#'     conditional on being treated, using `in_sample_counts` and `indep_counts`
#'     if available.
#'   \item **Bridge Regression:** Fits a bridge regression model using
#'     `grpreg::gBridge` on `X_final_scaled` and `y_final` with the specified `q`
#'     and lambda sequence.
#'   \item **Coefficient Selection (BIC):** Calls `getBetaBIC` to select the
#'     optimal `lambda` using BIC and retrieve the corresponding estimated
#'     coefficients (`theta_hat` in the transformed space).
#'   \item **Handle Zero-Feature Case:** If BIC selects a model with zero features,
#'     treatment effects are set to zero.
#'   \item **Coefficient Untransformation:** Calls `untransformCoefImproved` to
#'     transform `theta_hat` back to the original coefficient space, yielding
#'     `beta_hat`. If `add_ridge` was true, `beta_hat` is scaled.
#'   \item **Treatment Effect Calculation:**
#'     \itemize{
#'       \item Extracts cohort-specific average treatment effects (CATTs) from
#'         `beta_hat`.
#'       \item Calls `getCohortATTsFinal` to calculate CATT point estimates,
#'         standard errors (if `q < 1`), and confidence intervals. This involves
#'         computing the Gram matrix and related quantities.
#'     }
#'   \item **Overall ATT Calculation:** Calls `getTeResults2` to calculate the
#'     overall average treatment effect on the treated (ATT) and its standard
#'     error, using both in-sample probabilities and independent probabilities
#'     if `indep_counts` were provided.
#' }
#' The standard errors for CATTs are asymptotically exact. For ATT, if
#' `indep_counts` are provided, the SE is asymptotically exact; otherwise, it's
#' asymptotically conservative (if `q < 1`).
#'
#' @return A list containing detailed estimation results:
#'   \item{in_sample_att_hat}{Estimated overall ATT using in-sample cohort probabilities.}
#'   \item{in_sample_att_se}{Standard error for `in_sample_att_hat`.}
#'   \item{in_sample_att_se_no_prob}{SE for `in_sample_att_hat` ignoring variability from estimating cohort probabilities.}
#'   \item{indep_att_hat}{Estimated overall ATT using `indep_counts` cohort probabilities (NA if `indep_counts` not provided).}
#'   \item{indep_att_se}{Standard error for `indep_att_hat` (NA if not applicable).}
#'   \item{catt_hats}{A named vector of estimated CATTs for each cohort.}
#'   \item{catt_ses}{A named vector of SEs for `catt_hats` (NA if `q >= 1`).}
#'   \item{catt_df}{A data.frame summarizing CATTs, SEs, and confidence intervals.}
#'   \item{theta_hat}{The vector of estimated coefficients in the *transformed* (fused) space, including the intercept as the first element.}
#'   \item{beta_hat}{The vector of estimated coefficients in the *original* space (after untransforming `theta_hat`, excluding intercept).}
#'   \item{treat_inds}{Indices in `beta_hat` corresponding to base treatment effects.}
#'   \item{treat_int_inds}{Indices in `beta_hat` corresponding to treatment-covariate interactions.}
#'   \item{cohort_probs}{Estimated cohort probabilities conditional on being treated, from `in_sample_counts`.}
#'   \item{indep_cohort_probs}{Estimated cohort probabilities from `indep_counts` (NA if not provided).}
#'   \item{sig_eps_sq}{The (possibly estimated) variance of observation-level noise.}
#'   \item{sig_eps_c_sq}{The (possibly estimated) variance of unit-level random effects.}
#'   \item{lambda.max}{The maximum lambda value used in `grpreg`.}
#'   \item{lambda.max_model_size}{Model size for `lambda.max`.}
#'   \item{lambda.min}{The minimum lambda value used in `grpreg`.}
#'   \item{lambda.min_model_size}{Model size for `lambda.min`.}
#'   \item{lambda_star}{The lambda value selected by BIC.}
#'   \item{lambda_star_model_size}{Model size for `lambda_star`.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after fusion transformation and GLS weighting.}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, R, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
fetwfe_core <- function(
	X_ints,
	y,
	in_sample_counts,
	N,
	T,
	d,
	p,
	num_treats,
	first_inds,
	indep_counts = NA,
	sig_eps_sq = NA,
	sig_eps_c_sq = NA,
	lambda.max = NA,
	lambda.min = NA,
	nlambda = 100,
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE
) {
	ret <- check_etwfe_core_inputs(
		in_sample_counts = in_sample_counts,
		N = N,
		T = T,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		indep_counts = indep_counts,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge
	)

	R <- ret$R
	c_names <- ret$c_names
	indep_count_data_available <- ret$indep_count_data_available

	rm(ret)

	if (any(!is.na(lambda.max))) {
		stopifnot(is.numeric(lambda.max) | is.integer(lambda.max))
		stopifnot(length(lambda.max) == 1)
		stopifnot(lambda.max >= 0)
	}

	if (any(!is.na(lambda.min))) {
		stopifnot(is.numeric(lambda.min) | is.integer(lambda.min))
		stopifnot(length(lambda.min) == 1)
		stopifnot(lambda.min >= 0)
		if (any(!is.na(lambda.max))) {
			stopifnot(lambda.max >= lambda.min)
		}
	}

	stopifnot(is.numeric(q) | is.integer(q))
	stopifnot(length(q) == 1)
	stopifnot(q > 0)
	stopifnot(q <= 2)

	#
	#
	# Step 1: change coordinates of data so that regular bridge regression
	# penalty applied to transformed dataframe results in FETWFE penalty applied
	# to original data
	#
	#

	if (verbose) {
		message("Transforming matrix...")
	}

	# Transform matrix (change of coordinates so that fitting regular bridge
	# regression results in FETWFE fusion penalties)
	X_mod <- transformXintImproved(
		X_int = X_ints,
		N = N,
		T = T,
		R = R,
		d = d,
		num_treats = num_treats,
		first_inds = first_inds
	)

	res <- prep_for_etwfe_regresion(
		verbose = verbose,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p,
		num_treats = num_treats,
		add_ridge = add_ridge,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts
	)

	X_final_scaled <- res$X_final_scaled
	y_final <- res$y_final
	scale_center <- res$scale_center
	scale_scale <- res$scale_scale
	cohort_probs <- res$cohort_probs
	cohort_probs_overall <- res$cohort_probs_overall
	indep_cohort_probs <- res$indep_cohort_probs
	indep_cohort_probs_overall <- res$indep_cohort_probs_overall
	X_final <- res$X_final
	lambda_ridge <- res$lambda_ridge
	sig_eps_sq <- res$sig_eps_sq
	sig_eps_c_sq <- res$sig_eps_c_sq

	rm(res)

	#
	#
	# Step 4: estimate bridge regression and extract fitted coefficients
	#
	#

	# Estimate bridge regression
	if (verbose) {
		message("Estimating bridge regression...")
		t0 <- Sys.time()
	}

	if (!is.na(lambda.max) & !is.na(lambda.min)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.max = lambda.max,
			lambda.min = lambda.min,
			nlambda = nlambda
		)
	} else if (!is.na(lambda.max)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.max = lambda.max,
			nlambda = nlambda
		)
	} else if (!is.na(lambda.min)) {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			lambda.min = lambda.min,
			nlambda = nlambda
		)
	} else {
		fit <- grpreg::gBridge(
			X = X_final_scaled,
			y = y_final,
			gamma = q,
			nlambda = nlambda
		)
	}

	if (verbose) {
		message("Done! Time for estimation:")
		message(Sys.time() - t0)
	}

	# For diagnostics later, store largest and smallest lambda, as well as
	# corresponding smallest and largest model sizes, to return.
	lambda.max <- max(fit$lambda)
	lambda.max_model_size <- sum(fit$beta[, ncol(fit$beta)] != 0)

	lambda.min <- min(fit$lambda)
	lambda.min_model_size <- sum(fit$beta[, 1] != 0)

	# Select a single set of fitted coefficients by using BIC to choose among
	# the penalties that were fitted
	res <- getBetaBIC(
		fit = fit,
		N = N,
		T = T,
		p = p,
		X_mod = X_mod,
		y = y,
		scale_center = scale_center,
		scale_scale = scale_scale
	)

	theta_hat <- res$theta_hat # This includes intercept
	lambda_star_ind <- res$lambda_star_ind
	lambda_star_model_size <- res$lambda_star_model_size

	lambda_star <- fit$lambda[lambda_star_ind]

	# c_names <- names(in_sample_counts)[2:(R + 1)] # Moved definition up
	stopifnot(length(c_names) == R)

	# Indices corresponding to base treatment effects
	treat_inds <- getTreatInds(R = R, T = T, d = d, num_treats = num_treats)

	if (d > 0) {
		stopifnot(max(treat_inds) + 1 <= p)
		stopifnot(
			max(treat_inds) == R + T - 1 + d + R * d + (T - 1) * d + num_treats
		)

		treat_int_inds <- (max(treat_inds) + 1):p

		stopifnot(length(treat_int_inds) == num_treats * d)
	} else {
		stopifnot(max(treat_inds) <= p)
		stopifnot(max(treat_inds) == R + T - 1 + num_treats)

		treat_int_inds <- c()
	}

	# Handle edge case where no features are selected (model_size includes intercept)
	if (lambda_star_model_size <= 1 && all(theta_hat[2:(p + 1)] == 0)) {
		# only intercept might be non-zero
		if (verbose) {
			message(
				"No features selected (or only intercept); all treatment effects estimated to be 0."
			)
		}

		if (q < 1) {
			ret_se <- 0
		} else {
			ret_se <- NA
		}

		catt_df_to_ret <- data.frame(
			Cohort = c_names,
			`Estimated TE` = rep(0, R),
			SE = rep(ret_se, R),
			ConfIntLow = rep(ret_se, R),
			ConfIntHigh = rep(ret_se, R),
			check.names = FALSE
		)

		return(list(
			in_sample_att_hat = 0,
			in_sample_att_se = ret_se,
			in_sample_att_se_no_prob = ret_se,
			indep_att_hat = 0,
			indep_att_se = ret_se,
			catt_hats = setNames(rep(0, R), c_names),
			catt_ses = setNames(rep(ret_se, R), c_names),
			catt_df = catt_df_to_ret,
			theta_hat = theta_hat, # Includes intercept
			beta_hat = rep(0, p), # Slopes are all zero
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			lambda.max = lambda.max,
			lambda.max_model_size = lambda.max_model_size,
			lambda.min = lambda.min,
			lambda.min_model_size = lambda.min_model_size,
			lambda_star = lambda_star,
			lambda_star_model_size = lambda_star_model_size,
			X_ints = X_ints,
			y = y,
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			R = R,
			d = d,
			p = p,
			calc_ses = q < 1
		))
	}

	# estimated coefficients (slopes in transformed space)
	theta_hat_slopes <- theta_hat[2:(p + 1)]

	# Indices of selected features in transformed feature space (among slopes)
	sel_feat_inds <- which(theta_hat_slopes != 0)

	sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

	stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
	stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_treat_inds) <= length(treat_inds))
	stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

	# 1. Get theta_hat_slopes[treat_inds] -> these are the transformed treatment coefficients
	# 2. Find which of these are non-zero: `which(theta_hat_slopes[treat_inds] != 0)` -> these are
	# indices *within* the treat_inds block. Let's call this `sel_treat_inds_relative_to_block`.
	# `sel_treat_inds` itself contains the global indices of selected treatment features.
	# So `theta_hat_slopes[sel_treat_inds]` are the non-zero transformed treatment coefs.

	theta_hat_treat_block_transformed = theta_hat_slopes[treat_inds]
	sel_treat_inds_shifted <- which(theta_hat_treat_block_transformed != 0) # these are 1 to num_treats

	stopifnot(all(sel_treat_inds_shifted >= 1))
	stopifnot(all(sel_treat_inds_shifted <= num_treats))

	# Handle edge case where no treatment features selected
	if (length(sel_treat_inds_shifted) == 0) {
		if (verbose) {
			message(
				"No treatment features selected; all treatment effects estimated to be 0."
			)
		}

		if (q < 1) {
			ret_se <- 0
		} else {
			ret_se <- NA
		}

		catt_df_to_ret <- data.frame(
			Cohort = c_names,
			`Estimated TE` = rep(0, R),
			SE = rep(ret_se, R),
			ConfIntLow = rep(ret_se, R),
			ConfIntHigh = rep(ret_se, R),
			check.names = FALSE
		)

		# Need to untransform theta_hat_slopes to get beta_hat for consistency
		beta_hat_early_exit <- untransformCoefImproved(
			beta_hat_mod = theta_hat_slopes, # Pass slopes only
			first_inds = first_inds,
			T = T,
			R = R,
			p = p,
			d = d,
			num_treats = num_treats
		)
		if (add_ridge) {
			lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
			beta_hat_early_exit <- beta_hat_early_exit * (1 + lambda_ridge)
		}

		return(list(
			in_sample_att_hat = 0,
			in_sample_att_se = ret_se,
			in_sample_att_se_no_prob = ret_se,
			indep_att_hat = 0,
			indep_att_se = ret_se,
			catt_hats = setNames(rep(0, R), c_names),
			catt_ses = setNames(rep(ret_se, R), c_names),
			catt_df = catt_df_to_ret,
			theta_hat = theta_hat, # Full theta_hat with intercept
			beta_hat = beta_hat_early_exit, # Untransformed slopes
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			lambda.max = lambda.max,
			lambda.max_model_size = lambda.max_model_size,
			lambda.min = lambda.min,
			lambda.min_model_size = lambda.min_model_size,
			lambda_star = lambda_star,
			lambda_star_model_size = lambda_star_model_size,
			X_ints = X_ints,
			y = y,
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			R = R,
			d = d,
			p = p,
			calc_ses = q < 1
		))
	}

	#
	#
	# Step 5: transform estimated coefficients back to original feature
	# space
	#
	#

	beta_hat <- untransformCoefImproved(
		beta_hat_mod = theta_hat_slopes, # Pass slopes only
		first_inds = first_inds,
		T = T,
		R = R,
		p = p,
		d = d,
		num_treats = num_treats
	)

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
		beta_hat <- beta_hat * (1 + lambda_ridge)
	}

	# Get actual estimated treatment effects (in original, untransformed space)
	tes <- beta_hat[treat_inds]

	stopifnot(length(tes) == num_treats)
	# Checks based on transformed coefficients (theta_hat_slopes)
	stopifnot(all(theta_hat_slopes[treat_inds][sel_treat_inds_shifted] != 0))
	stopifnot(all(
		theta_hat_slopes[treat_inds][setdiff(
			1:num_treats,
			sel_treat_inds_shifted
		)] ==
			0
	))

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= num_treats)

	stopifnot(length(sel_feat_inds) > 0) # sel_feat_inds are indices in theta_hat_slopes
	stopifnot(length(sel_treat_inds_shifted) > 0) # sel_treat_inds_shifted are indices within the treat_inds block of theta_hat_slopes

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = sel_feat_inds, # Indices of non-zero elements in theta_hat_slopes
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted, # Indices (1 to num_treats) of non-zero transformed treat. coefs.
		c_names = c_names,
		tes = tes, # Untransformed treatment effect estimates (beta_hat[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		fused = TRUE,
		calc_ses = q < 1,
		p = p, # Total number of original parameters (columns in X_ints)
		alpha = alpha
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	d_inv_treat_sel <- res$d_inv_treat_sel
	calc_ses <- res$calc_ses

	rm(res)

	if (calc_ses) {
		stopifnot(nrow(d_inv_treat_sel) == num_treats)
		stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
	}

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	# theta_hat_treat_sel needs to be the selected non-zero *transformed* treatment coefficients
	# sel_treat_inds contains global indices of selected transformed features that are treatment effects
	theta_hat_treat_sel_for_att <- theta_hat_slopes[sel_treat_inds]

	in_sample_te_results <- getTeResults2(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		R = R,
		num_treats = num_treats,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_r | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		sel_treat_inds_shifted = sel_treat_inds_shifted,
		tes = tes, # Untransformed treatment effect estimates beta_hat[treat_inds]
		d_inv_treat_sel = d_inv_treat_sel,
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_r (unconditional on treated)
		first_inds = first_inds,
		theta_hat_treat_sel = theta_hat_treat_sel_for_att, # Selected non-zero transformed treat coefs
		calc_ses = calc_ses,
		indep_probs = FALSE
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob

	if ((q < 1) & calc_ses) {
		stopifnot(!is.na(in_sample_att_se))
	}

	if (indep_count_data_available) {
		indep_te_results <- getTeResults2(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			R = R,
			num_treats = num_treats,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_r | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			tes = tes,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_r (unconditional)
			first_inds = first_inds,
			theta_hat_treat_sel = theta_hat_treat_sel_for_att,
			calc_ses = calc_ses,
			indep_probs = TRUE
		)
		indep_att_hat <- indep_te_results$att_hat
		indep_att_se <- indep_te_results$att_te_se
	} else {
		indep_att_hat <- NA
		indep_att_se <- NA
	}

	return(list(
		in_sample_att_hat = in_sample_att_hat,
		in_sample_att_se = in_sample_att_se,
		in_sample_att_se_no_prob = in_sample_att_se_no_prob,
		indep_att_hat = indep_att_hat,
		indep_att_se = indep_att_se,
		catt_hats = cohort_tes, # Already named if applicable from getCohortATTsFinal
		catt_ses = cohort_te_ses, # Already named if applicable
		catt_df = cohort_te_df,
		theta_hat = theta_hat, # Full theta_hat (with intercept)
		beta_hat = beta_hat, # Untransformed slopes
		treat_inds = treat_inds,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda.max = lambda.max,
		lambda.max_model_size = lambda.max_model_size,
		lambda.min = lambda.min,
		lambda.min_model_size = lambda.min_model_size,
		lambda_star = lambda_star,
		lambda_star_model_size = lambda_star_model_size,
		X_ints = X_ints,
		y = y,
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p,
		calc_ses = calc_ses
	))
}


#-------------------------------------------------------------------------------
# Transformation Matrices and Coefficient Transformations
#-------------------------------------------------------------------------------

#' Transform a Full Design Matrix by  \( \boldsymbol D_N^{-1} \)
#'
#' @description
#' Takes the *raw* stacked panel design matrix
#' \eqn{\tilde{\boldsymbol Z}\in\mathbb R^{NT\times p}}
#' and post-multiplies it by the block-diagonal inverse fusion matrix
#' \(\boldsymbol D_N^{-1}\) from Lemma 3:
#' \deqn{
#'   \boldsymbol D_N^{-1}
#'   = \operatorname{diag}\!\bigl(
#'       (D^{(1)}(R))^{-1},\;
#'       (D^{(1)}(T-1))^{-1},\;
#'       I_{d},\;
#'       I_{d}\otimes(D^{(1)}(R))^{-1},\;
#'       I_{d}\otimes(D^{(1)}(T-1))^{-1},\;
#'       (D^{(2)}(\mathcal R))^{-1},\;
#'       I_{d}\otimes(D^{(2)}(\mathcal R))^{-1}
#'     \bigr).
#' }
#' The result is a matrix ready for vanilla bridge regression (Lasso,
#' elastic-net, \eqn{\ell_q} etc.), where the penalty on the transformed
#' coefficients reproduces the complex fusion penalty on the original ones.
#'
#' @details
#' The columns of `X_int` **must** appear in the order:
#' \enumerate{
#'   \item Cohort fixed-effects (length \eqn{R})
#'   \item Time fixed-effects (length \eqn{T-1})
#'   \item Main covariates (length \eqn{d})
#'   \item Covariate \(\times\) cohort interactions (length \eqn{dR})
#'   \item Covariate \(\times\) time interactions (length \eqn{d(T-1)})
#'   \item Base treatment effects (length `num_treats`)
#'   \item Covariate \(\times\) treatment interactions (length `d*num_treats`)
#' }
#' Each block is transformed exactly by the corresponding diagonal block of
#' \(\boldsymbol D_N^{-1}\) using helper functions
#' `genBackwardsInvFusionTransformMat()` and
#' `genInvTwoWayFusionTransformMat()`.
#'
#' Side-effect `stopifnot()` guards verify both the expected column order and
#' the absence of `NA`s after transformation.
#'
#' @param X_int Numeric matrix \eqn{NT\times p}.  Original design matrix.
#' @param N,T,R Integers. Panel dimensions and number of treated cohorts.
#' @param d Integer. Number of time-invariant covariates.
#' @param num_treats Integer. Total number of base treatment-effect dummies
#'   \eqn{\mathfrak W}.
#' @param first_inds Optional integer vector of length \eqn{R}.
#'   Starting column indices of the first treatment dummy of each cohort inside
#'   the treatment block.  If `NA` (default) they are computed by
#'   `getFirstInds()`.
#'
#' @return
#' A numeric matrix `X_mod` with the **same dimensions** as `X_int` but whose
#' columns are the transformed regressors
#' \(\bigl[\tilde{\boldsymbol Z}\,\boldsymbol D_N^{-1}\bigr]_{NT\times p}\).
#'
#' @seealso
#' * `genBackwardsInvFusionTransformMat()`
#' * `genInvTwoWayFusionTransformMat()`
#'
#' @examples
#' set.seed(1)
#' R <- 2; T <- 5; d <- 1; N <- 10
#' num_treats  <- getNumTreats(R, T)
#' p <- R + (T-1) + d + d*R + d*(T-1) + num_treats + d*num_treats
#' X_int <- matrix(rnorm(N*T*p), N*T, p)
#' X_mod <- transformXintImproved(
#'   X_int, N=N, T=T, R=R, d=d, num_treats=num_treats
#' )
#' # The two matrices have identical dimensions:
#' dim(X_mod)  # NT x p
#' @keywords internal
#' @noRd
transformXintImproved <- function(
	X_int,
	N,
	T,
	R,
	d,
	num_treats,
	first_inds = NA
) {
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)
	stopifnot(p == ncol(X_int))
	X_mod <- matrix(as.numeric(NA), nrow = N * T, ncol = p)
	stopifnot(nrow(X_int) == N * T)

	# Transform cohort fixed effects
	X_mod[, 1:R] <- X_int[, 1:R] %*% genBackwardsInvFusionTransformMat(R)

	# Transform time fixed effects
	X_mod[, (R + 1):(R + T - 1)] <- X_int[, (R + 1):(R + T - 1)] %*%
		genBackwardsInvFusionTransformMat(T - 1)

	# Copy X (the main covariate block; may be empty when d==0)
	if (d > 0) {
		stopifnot(all(is.na(X_mod[, (R + T - 1 + 1):(R + T - 1 + d)])))
		X_mod[, (R + T - 1 + 1):(R + T - 1 + d)] <- X_int[,
			(R + T - 1 + 1):(R + T - 1 + d)
		]

		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d)])))
		stopifnot(all(is.na(X_mod[, (R + T - 1 + d + 1):p])))
	}

	# For cohort effects interacted with X: we have d*R columns to deal with.
	# For each individual feature, this will be handled using
	# genTransformedMatFusion.
	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(R = R, T = T)
	}

	if (d > 0) {
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on R times
			feat_1 <- R + T - 1 + d + j
			feat_R <- R + T - 1 + d + (R - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == R)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))

			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(R)
		}
		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R * d)])))
		stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- R + T - 1 + d + R * d + j
			feat_T_minus_1 <- R + T - 1 + d + R * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))
			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(T - 1)
		}
		stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R * d + (T - 1) * d)])))
		stopifnot(all(is.na(X_mod[,
			(R + T - 1 + d + R * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects. For each cohort, will penalize base term, then
	# fuse remaining terms toward it. Also, for each cohort, will penalize base
	# treatment effect of this cohort to base of previous cohort. New function
	# genTransformedMatTwoWayFusion does this.

	feat_inds <- (R + T - 1 + d + R * d + (T - 1) * d + 1):(R +
		T -
		1 +
		d +
		R * d +
		(T - 1) * d +
		num_treats)

	# Now ready to generate the appropriate transformed matrix
	stopifnot(all(is.na(X_mod[, feat_inds])))

	X_mod[, feat_inds] <- X_int[, feat_inds] %*%
		genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

	if (d > 0) {
		# Lastly, penalize interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with genTransformedMatTwoWayFusion, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + R + T - 1 + d + R * d + (T - 1) * d + num_treats

			# Now ready to generate the appropriate transformed matrix
			stopifnot(all(is.na(X_mod[, inds_j])))

			X_mod[, inds_j] <- X_int[, inds_j] %*%
				genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

			stopifnot(all(!is.na(X_mod[, inds_j])))
		}
	}

	stopifnot(all(!is.na(X_mod)))
	stopifnot(ncol(X_mod) == p)
	stopifnot(nrow(X_mod) == N * T)

	return(X_mod)
}

#' Back-transform Bridge-Regression Coefficients
#' \( \widehat{\boldsymbol\beta}
#'   = \boldsymbol D_N^{-1}\,\widehat{\boldsymbol\theta}\)
#'
#' @description
#' After fitting a bridge (Lasso, Elastic-Net, \eqn{\ell_q}) regression on the
#' transformed design matrix
#' \eqn{\widetilde{\boldsymbol Z}\,\boldsymbol D_N^{-1}}
#' (see `transformXintImproved()`), the solver returns parameter estimates
#' \eqn{\widehat{\boldsymbol\theta}}.
#' This helper multiplies that vector by the block-diagonal matrix
#' \eqn{\boldsymbol D_N^{-1}} from the paper,
#' thereby recovering the original-scale coefficients
#' \eqn{\widehat{\boldsymbol\beta}} for the FETWFE model.
#'
#' @details
#' The matrix
#' \deqn{
#' \boldsymbol D_N^{-1}
#'  = \operatorname{diag}\Bigl(
#'        (D^{(1)}(R))^{-1},\;
#'        (D^{(1)}(T\!-\!1))^{-1},\;
#'        I_d,\;
#'        I_d\!\otimes\!(D^{(1)}(R))^{-1},\;
#'        I_d\!\otimes\!(D^{(1)}(T\!-\!1))^{-1},\;
#'        (D^{(2)}(\mathcal R))^{-1},\;
#'        I_d\!\otimes\!(D^{(2)}(\mathcal R))^{-1}
#'     \Bigr)
#' }
#' corresponds to seven consecutive blocks in the column ordering used by
#' `transformXintImproved()`:
#' \enumerate{
#'   \item Cohort fixed-effect coefficients (length \eqn{R})
#'   \item Time fixed-effects (length \eqn{T-1})
#'   \item Main covariates (length \eqn{d})
#'   \item Covariate x cohort interactions (\eqn{dR})
#'   \item Covariate x time interactions \eqn{d(T-1)}
#'   \item Base treatment effects (\eqn{\mathfrak W =} `num_treats`)
#'   \item Covariate x treatment interactions (\eqn{d\mathfrak W})
#' }
#'
#' For each block the function premultiplies the slice of
#' \code{beta_hat_mod} with the *same* inverse-fusion matrix that was used as a
#' post-multiplier when building the transformed design matrix:
#'
#' | block | helper used | mathematical symbol |
#' |-------|-------------|---------------------|
#' | cohort FE | `genBackwardsInvFusionTransformMat(R)` | \((D^{(1)}(R))^{-1}\) |
#' | time FE | `genBackwardsInvFusionTransformMat(T-1)` | \((D^{(1)}(T-1))^{-1}\) |
#' | covariate blocks | identity copy | \(I_d\) |
#' | covariate x cohort | same helper, repeated for each feature | \(I_d\otimes(D^{(1)}(R))^{-1}\) |
#' | covariate x time | idem with \(T-1\) | \(I_d\otimes(D^{(1)}(T-1))^{-1}\) |
#' | base treatment | `genInvTwoWayFusionTransformMat()` | \((D^{(2)}(\mathcal R))^{-1}\) |
#' | covariate x treatment | same two-way helper for each feature | \(I_d\otimes(D^{(2)}(\mathcal R))^{-1}\) |
#'
#' @param beta_hat_mod Numeric vector (length \eqn{p}).
#'   Estimated coefficients returned by a penalised fit on the transformed
#'   design matrix.
#' @param T Integer. Total time periods.
#' @param R Integer. Number of treated cohorts.
#' @param p Integer. Total number of columns in the **original** design matrix
#'   \eqn{p = p_N}.
#' @param d Integer. Number of time-invariant covariates.
#' @param num_treats Integer. Count of base treatment-effect dummies
#'   \eqn{\mathfrak W}.
#' @param first_inds Optional integer vector of length \eqn{R}.
#'   Starting indices (1-based, inside the treatment-dummy block) of the first
#'   effect for each cohort.  If `NA` they are reconstructed with
#'   **`getFirstInds()`**.
#'
#' @return Numeric vector \code{beta_hat} (length \eqn{p}) equal to
#'   \eqn{\boldsymbol D_N^{-1}\,\widehat{\boldsymbol\theta}}.
#'
#' @references
#' Wooldridge, J. M. (2021). Two-way fixed effects, the two-way mundlak
#' regression, and difference-in-differences estimators.
#' \emph{Available at SSRN 3906345}.
#' \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3906345}.
#' Tibshirani & Taylor (2011), "The Solution Path of the Generalized
#' Lasso".
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @seealso
#' * `transformXintImproved()` - forward transformation of the design matrix.
#' * `genBackwardsInvFusionTransformMat()`,
#'   `genInvTwoWayFusionTransformMat()` - block-wise inverse matrices.
#'
#' @examples
#' ## toy example: one covariate, two treated cohorts, T = 5
#' R <- 2; T <- 5; d <- 1
#' num_treats <- getNumTreats(R, T)
#' p <- R + (T-1) + d + d*R + d*(T-1) + num_treats + d*num_treats
#'
#' ## pretend we ran a penalised regression in the transformed space
#' beta_hat_mod <- rnorm(p)
#'
#' ## back-transform:
#' beta_hat <- untransformCoefImproved(
#'   beta_hat_mod, T=T, R=R, p=p, d=d,
#'   num_treats=num_treats
#' )
#' length(beta_hat)  # == p
#' @keywords internal
#' @noRd
untransformCoefImproved <- function(
	beta_hat_mod,
	T,
	R,
	p,
	d,
	num_treats,
	first_inds = NA
) {
	stopifnot(length(beta_hat_mod) == p)
	beta_hat <- rep(as.numeric(NA), p)

	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(R = R, T = T)
	}

	# First handle R cohort fixed effects effects
	beta_hat[1:R] <- genBackwardsInvFusionTransformMat(R) %*% beta_hat_mod[1:R]

	stopifnot(all(!is.na(beta_hat[1:R])))
	stopifnot(all(is.na(beta_hat[(R + 1):p])))

	# Next, T - 1 time fixed effects
	beta_hat[(R + 1):(R + T - 1)] <- genBackwardsInvFusionTransformMat(
		T - 1
	) %*%
		beta_hat_mod[(R + 1):(R + T - 1)]

	stopifnot(all(!is.na(beta_hat[1:(R + T - 1)])))
	stopifnot(all(is.na(beta_hat[(R + T):p])))

	# Coefficients for X (if any)
	if (d > 0) {
		beta_hat[(R + T):(R + T - 1 + d)] <- beta_hat_mod[
			(R + T):(R + T - 1 + d)
		]

		stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d)])))
		stopifnot(all(is.na(beta_hat[(R + T + d):p])))

		# Next, coefficients for cohort effects interacted with X. For each individual
		# feature, this will be handled using untransformVecFusion.
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on R times
			feat_1 <- R + T - 1 + d + j
			feat_R <- R + T - 1 + d + (R - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == R)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(R) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d + R * d)])))
		stopifnot(all(is.na(beta_hat[(R + T - 1 + d + R * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- R + T - 1 + d + R * d + j
			feat_T_minus_1 <- R + T - 1 + d + R * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(
				T - 1
			) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(
			!is.na(beta_hat[1:(R + T - 1 + d + R * d + (T - 1) * d)])
		))
		stopifnot(all(is.na(beta_hat[
			(R + T - 1 + d + R * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects.

	feat_inds <- (R + T - 1 + d + R * d + (T - 1) * d + 1):(R +
		T -
		1 +
		d +
		R * d +
		(T - 1) * d +
		num_treats)

	stopifnot(all(is.na(beta_hat[feat_inds])))

	beta_hat[feat_inds] <- genInvTwoWayFusionTransformMat(
		num_treats,
		first_inds,
		R
	) %*%
		beta_hat_mod[feat_inds]

	if (d > 0) {
		stopifnot(all(
			!is.na(beta_hat[
				1:(R + T - 1 + d + R * d + (T - 1) * d + num_treats)
			])
		))
		stopifnot(all(is.na(beta_hat[
			(R + T - 1 + d + R * d + (T - 1) * d + num_treats + 1):p
		])))
		# Lastly, interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with untransformTwoWayFusionCoefs, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + R + T - 1 + d + R * d + (T - 1) * d + num_treats

			# Now ready to untransform the estimated coefficients
			stopifnot(all(is.na(beta_hat[inds_j])))

			beta_hat[inds_j] <- genInvTwoWayFusionTransformMat(
				num_treats,
				first_inds,
				R
			) %*%
				beta_hat_mod[inds_j]

			stopifnot(all(!is.na(beta_hat[inds_j])))
		}
	}

	stopifnot(all(!is.na(beta_hat)))

	return(beta_hat)
}

# genBackwardsFusionTransformMat
#' @title Generate Backward Fusion Transformation Matrix
#' @description Creates a square transformation matrix `D` of size `n_vars` x
#'   `n_vars`. When pre-multiplied by a coefficient vector `beta`, `D %*% beta`
#'   yields a transformed vector `theta` where `theta_i = beta_i - beta_{i+1}`
#'   for `i < n_vars`, and `theta_{n_vars} = beta_{n_vars}`. This is used to
#'   penalize coefficients towards the *next* coefficient in sequence, and the
#'   last coefficient directly.
#' @param n_vars Integer; the number of variables (coefficients) in the block
#'   to be transformed. This will be the dimension of the output matrix.
#' @return A numeric matrix of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix `D` has 1s on the main diagonal. For each row
#'   `i` (from 1 to `n_vars - 1`), it has a -1 at column `i+1`. All other
#'   elements are 0.
#' @examples
#'   genBackwardsFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1   -1    0
#'   # [2,]    0    1   -1
#'   # [3,]    0    0    1
#' @keywords internal
#' @noRd
genBackwardsFusionTransformMat <- function(n_vars) {
	# Generates D matrix in relation theta = D beta, where D beta is what
	# we want to penalize (for a single set of coefficients where we want to
	# penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D <- matrix(0, n_vars, n_vars)

	for (i in 1:n_vars) {
		for (j in 1:n_vars) {
			if (i == j) {
				D[i, j] <- 1
			}
			if (j == i + 1) {
				D[i, j] <- -1
			}
		}
	}

	return(D)
}

#' Inverse "Backwards-Difference" Transformation Matrix  \( \bigl(D^{(1)}(t)\bigr)^{-1} \)
#'
#' @description
#' Generates the \eqn{t\times t} upper-triangular matrix of 1s whose inverse is
#' the first-difference operator
#' \eqn{D^{(1)}(t)}
#' defined in Eq. (14) of the paper:
#' \deqn{
#'   D^{(1)}(t) =
#'   \begin{bmatrix}
#'     1 & -1 &        &        & 0\\
#'       &  1 & -1     &        &  \\
#'       &    & \ddots & \ddots &  \\
#'       &    &        & 1 & -1 \\
#'     0 &    &        &   &  1
#'   \end{bmatrix}.
#' }
#' Multiplying a block of coefficients by this inverse converts finite
#' differences back to cumulative sums, which is what the bridge-penalty
#' expects.
#'
#' @details
#' The matrix has ones on and strictly above the main diagonal
#' (\eqn{D^{-1}_{ij}=1_{\{i\le j\}}}).
#' It is its own Cholesky factor, \eqn{D^{-1} = U = U^\top}.
#'
#' After constructing `D_inv`, the routine calls
#' `genBackwardsFusionTransformMat(n_vars)`, forms both
#' `D %*% D_inv` and `D_inv %*% D`, and checks that their
#' maximum element-wise deviation from the identity matrix is below
#' `tol = 1e-12`.  Failing the check raises an error, protecting downstream
#' computations from a silent algebraic bug.
#'
#' @param n_vars Integer. Dimension \eqn{t}.
#'
#' @return A \eqn{n\_vars \times n\_vars} numeric matrix equal to
#'   \(\bigl(D^{(1)}(n\_vars)\bigr)^{-1}\).
#' @keywords internal
#' @noRd
genBackwardsInvFusionTransformMat <- function(n_vars) {
	# Generates inverse of D matrix in relation theta = D beta, where D beta is
	# what we want to penalize (for a single set of coefficients where we want
	# to penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D_inv <- matrix(0, n_vars, n_vars)

	diag(D_inv) <- 1

	D_inv[upper.tri(D_inv)] <- 1

	stopifnot(nrow(D_inv) == n_vars)
	stopifnot(ncol(D_inv) == n_vars)

	## -- self-test: D_inv %*% D  ==  I
	D <- genBackwardsFusionTransformMat(n_vars)
	tol <- 1e-12
	if (!all.equal(D_inv %*% D, diag(n_vars), tolerance = tol)) {
		stop(
			"genBackwardsInvFusionTransformMat(): self-test failed - ",
			"result is not the matrix inverse of D."
		)
	}

	return(D_inv)
}

# genInvFusionTransformMat (TODO: unused for now)
#' @title Generate Inverse of Forward Fusion Transformation Matrix
#' @description Creates the inverse of a "forward" fusion transformation matrix.
#'   A forward fusion matrix `D_forward` (not explicitly generated here) would
#'   transform `beta` to `theta` such that `theta_1 = beta_1` and
#'   `theta_i = beta_i - beta_{i-1}` for `i > 1`. This function returns
#'   `D_forward_inv`. If `theta = D_forward %*% beta`, then
#'   `beta = D_forward_inv %*% theta`. This is used when the first coefficient
#'   in a sequence is penalized directly, and subsequent coefficients are
#'   penalized towards the *previous* coefficient.
#' @param n_vars Integer; the number of variables (coefficients), determining
#'   the dimension of the output matrix.
#' @return A numeric matrix `D_forward_inv` of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix is a lower triangular matrix with all elements
#'   on and below the main diagonal equal to 1, and all elements above the main
#'   diagonal equal to 0.
#' @examples
#'   genInvFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1    0    0
#'   # [2,]    1    1    0
#'   # [3,]    1    1    1
#' @keywords internal
#' @noRd
genInvFusionTransformMat <- function(n_vars) {
	# Generates inverse of D matrix in relation theta = D beta, where D beta is
	# what we want to penalize (for a single set of coefficients where we want
	# to penalize first coefficient directly and penalize remaining coefficints
	# towards the previous coefficient)
	D_inv <- matrix(0, n_vars, n_vars)

	diag(D_inv) <- 1

	D_inv[lower.tri(D_inv)] <- 1

	return(D_inv)
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
#'         (\(\widehat\pi_r(1-\widehat\pi_r)\) on the diagonal,
#'         \(-\widehat\pi_r\widehat\pi_s\) off-diagonal).
#'   \item \(J\) is the Jacobian of the weighting function
#'         \(f_r(\pi)=\pi_r/\sum_{k}\pi_k\) evaluated at
#'         \(\widehat\pi\).  In matrix form
#'         \eqn{J_{rs}=
#'           \begin{cases}
#'             (1-\widehat\pi_r)/S^2,& r=s,\\[4pt]
#'             -\,\widehat\pi_r    /S^2,& r\neq s,
#'           \end{cases}} with \(S=\sum_k\widehat\pi_k\).
#'   \item Each column of `d_inv_treat_sel` is a selected **transformed**
#'         treatment coefficient; row-averaging over the rows belonging to
#'         cohort \(s\) produces the part of \(J\) that multiplies
#'         \(\theta_{\text{sel}}\).
#' }
#' @param psi_mat Numeric matrix; a matrix where each column `r` is the `psi_r`
#'   vector used in calculating the ATT for cohort `r`. Dimensions:
#'   `length(sel_treat_inds_shifted)` x `R`.
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#' @param tes Numeric vector; the estimated treatment effects for all
#'   `num_treats` possible cohort-time combinations.
#' @param cohort_probs_overall Numeric vector; estimated marginal probabilities
#'   of belonging to each treated cohort (P(W=r)). Length `R`.
#' @param first_inds Integer vector; indices of the first treatment effect for
#'   each cohort within the `num_treats` block.
#' @param theta_hat_treat_sel Numeric vector; estimated coefficients in the
#'   transformed (fused) space, corresponding only to the selected treatment
#'   effects.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param R Integer; total number of treated cohorts.
#' @param fused Logical; if `TRUE`, assumes fusion penalization was used,
#'   affecting how standard errors and related matrices are computed.
#' @param d_inv_treat_sel Numeric matrix; the relevant block of the inverse
#'   two-way fusion transformation matrix corresponding to selected treatment
#'   effects. Dimensions: `num_treats` (or fewer if selection occurs) x
#'   `length(sel_treat_inds_shifted)`. Does not need to be provided if fused
#'   is FALSE.
#' @return A numeric scalar representing the second variance component for the
#'   ATT.
#' @details This function calculates `Sigma_pi_hat`, the covariance matrix of
#'   the cohort assignment indicators, and a Jacobian matrix. These are then
#'   combined with `theta_hat_treat_sel` to compute the variance term as
#'   `T * t(theta_hat_treat_sel) %*% t(jacobian_mat) %*% Sigma_pi_hat %*% jacobian_mat %*%
#'   theta_hat_treat_sel / (N * T)`. The construction of the Jacobian involves averaging parts of
#'   `d_inv_treat_sel` corresponding to different cohorts.
#'
#' The function is currentky written for the fused workflow only
#' (`fused = TRUE` guard).  It first reconstructs \(J\) from
#' `d_inv_treat_sel` and the cohort probabilities, then plugs everything into
#' the quadratic form above and finally rescales by \(T/(N T)=1/N\).
#' @inheritParams getSecondVarTermDataApp
#' @seealso [getTeResults2()]
#' @keywords internal
#' @noRd
getSecondVarTermDataApp <- function(
	# cohort_probs,
	psi_mat,
	sel_treat_inds_shifted,
	tes,
	cohort_probs_overall,
	first_inds,
	theta_hat_treat_sel,
	num_treats,
	N,
	T,
	R,
	fused = TRUE,
	d_inv_treat_sel = NA
) {
	if (fused) {
		stopifnot(all(!is.na(d_inv_treat_sel)))
		stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
	}
	stopifnot(length(theta_hat_treat_sel) == length(sel_treat_inds_shifted))

	# Get Sigma_pi_hat, the (sample-estimated) covariance matrix for the
	# sample proportions (derived from the multinomial distribution)
	Sigma_pi_hat <- -outer(
		cohort_probs_overall[1:(R)],
		cohort_probs_overall[1:(R)]
	)
	diag(Sigma_pi_hat) <- cohort_probs_overall[1:(R)] *
		(1 - cohort_probs_overall[1:(R)])

	stopifnot(nrow(Sigma_pi_hat) == R)
	stopifnot(ncol(Sigma_pi_hat) == R)

	# Jacobian
	##
	## J_{rs} =  (S-pi_r)/S^2      if r = s
	##           -pi_r   /S^2      if r != s
	##
	## where S := sum(cohort_probs_overall)
	##
	## Each column block of d_inv_treat_sel belongs to a selected theta coordinate;
	## averaging the rows of cohorts gives the vector needed to multiply theta_sel.
	##
	jacobian_mat <- matrix(
		as.numeric(NA),
		nrow = R,
		ncol = length(sel_treat_inds_shifted)
	)

	if (fused) {
		# Gather a list of the indices corresponding to the treatment coefficients
		# for each cohort
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
				stopifnot(length(sel_inds[[r]]) < length(sel_inds[[r - 1]]))
			}
		}
		stopifnot(all.equal(unlist(sel_inds), 1:num_treats))

		for (r in 1:R) {
			## diagonal contribution
			cons_r <- (sum(cohort_probs_overall) -
				cohort_probs_overall[r]) /
				sum(cohort_probs_overall)^2

			if (length(sel_treat_inds_shifted) > 1) {
				jacobian_mat[r, ] <- cons_r *
					colMeans(d_inv_treat_sel[sel_inds[[r]], , drop = FALSE])
			} else {
				jacobian_mat[r, ] <- cons_r *
					mean(d_inv_treat_sel[sel_inds[[r]], , drop = FALSE])
			}

			## off-diagonal: subtract sum_{s!=r} pi_r / S^2  x  block-mean of cohort s
			for (r_double_prime in setdiff(1:R, r)) {
				cons_r_double_prime <- cohort_probs_overall[r] /
					sum(cohort_probs_overall)^2

				jacobian_mat[r, ] <- jacobian_mat[r, ] -
					cons_r_double_prime *
						colMeans(d_inv_treat_sel[
							sel_inds[[r_double_prime]],
							,
							drop = FALSE
						])
			}
		}
	}

	stopifnot(all(!is.na(jacobian_mat)))

	## variance term: theta_sel' J' sum_pi J theta_sel / N
	if (fused) {
		att_var_2 <- T *
			as.numeric(
				t(theta_hat_treat_sel) %*%
					t(jacobian_mat) %*%
					Sigma_pi_hat %*%
					jacobian_mat %*%
					theta_hat_treat_sel
			) /
			(N * T)
	}

	if (!fused) {
		stop("this function is not implemented as of now.")
	}

	return(att_var_2)
}

#' Inverse Two-Way-Fusion Transformation Matrix  \( \bigl(D^{(2)}(\mathcal R)\bigr)^{-1} \)
#'
#' @description
#' Constructs the **inverse** of the block-lower-triangular matrix
#' \eqn{D^{(2)}(\mathcal R)} that appears in Lemma \eqn{3} of the paper.
#' Each treated cohort \(r\) has a "first" post-treatment coefficient
#' \(\tau_{r,0}\) and a string of subsequent coefficients
#' \(\tau_{r,1},\dots,\tau_{r,T-r}\).
#' The fusion penalty (1) pulls every \(\tau_{r,k}\;(k>0)\) toward its
#' predecessor \(\tau_{r,k-1}\) *within* the same cohort **and**
#' (2) pulls the first effect of cohort \(r\) toward the first effect of cohort
#' \(r-1\).
#' Multiplying the raw design sub-matrix by the output of
#' `genInvTwoWayFusionTransformMat()` therefore changes coordinates from
#' \eqn{\boldsymbol\beta} to
#' \eqn{\boldsymbol\theta = D^{(2)}(\mathcal R)\,\boldsymbol\beta},
#' so that an \eqn{\ell_q} penalty on \eqn{\theta} is exactly the desired
#' two-level fusion penalty on \eqn{\beta}.
#'
#' @details
#' * The returned matrix is **upper-triangular with 1s** on and above the main
#'   diagonal and zeros elsewhere, except for extra 1s that implement the
#'   cross-cohort link on the first effect of each cohort (see paper, Eq. (18)).
#' * Its inverse contains only \{-1,0,1\} and recreates Eq. (17) of the paper.
#' * Determinant is 1, so the transformation is volume-preserving.
#'
#' @param n_vars Integer. Total number of base treatment-effect coefficients
#'   ( \eqn{\mathfrak W}  in the paper).
#' @param first_inds Integer vector of length \eqn{R}.
#'   `first_inds[r]` is the **1-based** column index of \(\tau_{r,0}\)
#'   inside the block of the \eqn{n\_vars} treatment columns.
#' @param R Integer. Number of treated cohorts.
#'
#' @return A numeric matrix of size \eqn{n\_vars \times n\_vars} that is
#'   exactly \(\bigl(D^{(2)}(\mathcal R)\bigr)^{-1}\).
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' R  <- 3;  T <- 6
#' num_treats <- getNumTreats(R, T)
#' first <- getFirstInds(R, T)
#' Dinv <- genInvTwoWayFusionTransformMat(num_treats, first, R)
#' # verify Dinv %*% solve(Dinv) = I
#' all.equal(Dinv %*% solve(Dinv), diag(num_treats))
#' @keywords internal
#' @noRd
genInvTwoWayFusionTransformMat <- function(n_vars, first_inds, R) {
	stopifnot(length(n_vars) == 1, n_vars >= 0)
	stopifnot(is.numeric(R), length(R) == 1, R >= 0)
	if (R > 0) {
		stopifnot(is.numeric(first_inds), length(first_inds) == R)
		# Add checks for consistency of first_inds if n_vars > 0
		if (n_vars > 0) {
			stopifnot(first_inds[1] == 1)
			if (R > 1) {
				stopifnot(all(diff(first_inds) > 0)) # Must be strictly increasing
			}
			# Check that the last effect of the last cohort aligns with n_vars
			# M_R = n_vars - first_inds[R] + 1. This M_R must be > 0 if first_inds[R] <= n_vars.
			stopifnot(first_inds[R] <= n_vars)
		} else {
			# n_vars == 0
			stopifnot(R == 0) # If n_vars is 0, R must be 0. first_inds should be empty.
		}
	} else {
		# R == 0
		stopifnot(length(first_inds) == 0, n_vars == 0)
	}

	D_inv <- matrix(0, nrow = n_vars, ncol = n_vars)

	if (n_vars == 0) {
		# Handles R=0 correctly
		return(D_inv)
	}

	# Part 1: Set up the column structure for \tilde{U} blocks and
	# the first column of each diagonal block.
	# For each cohort `c` (from 1 to R), the column in D_inv corresponding to its
	# first treatment effect (i.e., absolute column index `first_inds[c]`)
	# should have 1s from its own row (`first_inds[c]`) down to `n_vars`.
	if (R > 0) {
		for (cohort_idx in 1:R) {
			col_to_fill <- first_inds[cohort_idx]
			D_inv[col_to_fill:n_vars, col_to_fill] <- 1
		}
	}

	# Part 2: Form the correct diagonal blocks.
	# Each diagonal block is (D^{(1)}(M_c)^{-1})^T, which is a fully
	# lower triangular matrix of 1s (all elements on and below diagonal are 1).
	# This will correctly overwrite the 1s placed in Part 1 for these diagonal parts.
	if (R > 0) {
		for (cohort_idx in 1:R) {
			block_start_idx <- first_inds[cohort_idx]

			if (cohort_idx < R) {
				block_end_idx <- first_inds[cohort_idx + 1] - 1
			} else {
				# This is the last cohort
				block_end_idx <- n_vars
			}

			# Ensure block indices are valid (e.g. cohort has at least one effect)
			if (block_start_idx > block_end_idx) {
				# This case implies cohort_idx has zero effects.
				# This should ideally not happen if first_inds and n_vars are consistent
				# with cohorts having at least one effect.
				# If it can happen, 'continue' or 'next' might be appropriate.
				# For now, assume M_x >= 1 for all cohorts included in R.
				next
			}

			# Fill this diagonal block D_inv[block_start_idx:block_end_idx, block_start_idx:block_end_idx]
			for (row_abs in block_start_idx:block_end_idx) {
				# Absolute row index in D_inv
				# For the current row_abs within its block, set columns from
				# the start of the block (block_start_idx) up to the current row_abs to 1.
				D_inv[row_abs, block_start_idx:row_abs] <- 1
			}
		}
	}

	# The original code had a stop for R < 2.
	# This revised logic handles R=0 and R=1 correctly without a special stop.
	# If R=1, Part 1 sets D_inv[1:n_vars, 1] <- 1.
	# Part 2 (with cohort_idx=1, block_start_idx=1, block_end_idx=n_vars) then correctly
	# forms the full (D^{(1)}(n_vars)^{-1})^T matrix.

	return(D_inv)
}
