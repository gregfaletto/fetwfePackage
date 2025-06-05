#' @title Prepare Transformed Design Matrix and Response for ETWFE Regression
#'
#' @description
#' Generates all matrix- and vector-level inputs required by the bridge/OLS
#' fitting step inside \code{\link{etwfe_core}}.
#' Its responsibilities include: estimating or assembling the covariance
#' matrix `\Omega`; performing the GLS whitening transformation;
#' (optionally) appending ridge-penalty rows; computing cohort probability
#' weights; and returning both raw and scaled versions of the final design
#' matrix.
#'
#' @param verbose Logical.  If \code{TRUE}, prints timing information for
#'   each major sub-task.  Default \code{FALSE}.
#' @param sig_eps_sq,sig_eps_c_sq Numeric or \code{NA}.  Observation-level
#'   and unit-level variance components.  If either is \code{NA},
#'   they are estimated via \code{estOmegaSqrtInv()}.
#' @param y Numeric vector of length \(N\times T\).  Centered response.
#' @param X_ints Numeric matrix.  Full design matrix in the \emph{original}
#'   parameterisation.
#' @param X_mod Numeric matrix.  Same dimensions as \code{X_ints}.  In the
#'   ETWFE path this is usually a transformed version, but for OLS it is
#'   often identical to \code{X_ints}.
#' @param N,T,R,d,p,num_treats Integers giving key problem dimensions:
#'   number of units, time periods, treated cohorts, covariates,
#'   total parameters, and base treatment-effect parameters, respectively.
#' @param add_ridge Logical.  Whether to append rows that implement a small
#'   L2 penalty on the \emph{untransformed} coefficients.  Default \code{FALSE}.
#' @param first_inds Integer vector (length \(R\)).  Column indices of the
#'   first base treatment-effect for each cohort within the block of
#'   \code{num_treats} columns.
#' @param in_sample_counts Integer vector of cohort sizes in the estimation
#'   sample (see \code{check_etwfe_core_inputs}).
#' @param indep_count_data_available Logical.  Indicates whether
#'   \code{indep_counts} contains valid cohort counts from an independent
#'   sample; affects the probability calculations below.
#' @param indep_counts Optional integer vector.  Cohort counts from an
#'   independent sample; only used when
#'   \code{indep_count_data_available = TRUE}.  Default \code{NA}.
#' @param is_fetwfe Logical.  If \code{TRUE}, a fusion transformation matrix
#'   has been applied upstream (this matters for how the optional ridge rows
#'   are constructed).  Default \code{TRUE}.
#' @param is_twfe_covs Logical.  If \code{TRUE}, columns will be removed and
#'   consolidated as required for `twfeCovs()`.  Default \code{FALSE}.
#'
#' @details
#' The routine carries out the following steps in order:
#' \enumerate{
#'   \item \strong{Variance-Component Handling}\newline
#'     Calls \code{estOmegaSqrtInv()} when either variance component is
#'     unknown; otherwise uses supplied values.
#'   \item \strong{GLS Transformation}\newline
#'     Forms \(\Omega = \sigma_\varepsilon^2 I_T + \sigma_{\varepsilon c}^2 J_T\)
#'     and multiplies both \code{y} and \code{X_mod} on the left by
#'     \(\sqrt{\sigma_\varepsilon^2}\,\Omega^{-1/2} \otimes I_N\).
#'   \item \strong{Optional Ridge Augmentation}\newline
#'     If \code{add_ridge = TRUE}, appends \(p\) rows to the scaled design
#'     matrix and response, thereby imposing a tiny L2 penalty equal to
#'     \code{lambda_ridge}.  When \code{is_fetwfe = TRUE} the rows are first
#'     premultiplied by the inverse fusion transformation so the penalty
#'     applies on the original coefficient scale.
#'   \item \strong{Cohort Probability Estimation}\newline
#'     Computes \eqn{\hat\pi_r = n_r / \sum_{s=1}^R n_s} from the in-sample
#'     counts; when an independent split is available the same is done for
#'     \code{indep_counts}.
#' }
#'
#' @return A list with everything \code{etwfe_core} needs next:
#'   \describe{
#'     \item{\code{X_final_scaled}}{Design matrix after GLS transform
#'       and column-wise scaling (and ridge rows if requested).}
#'     \item{\code{X_final}}{Design matrix after GLS transform but before scaling.}
#'     \item{\code{y_final}}{Response after GLS transform (and augmentation).}
#'     \item{\code{scale_center}, \code{scale_scale}}{Vectors used to undo
#'       column scaling.}
#'     \item{\code{cohort_probs}}{Vector of \(\hat\pi_r \mid \text{treated}\)
#'       from the estimation sample.}
#'     \item{\code{cohort_probs_overall}}{Vector of unconditional
#'       probabilities \(P(W=r)\).}
#'     \item{\code{indep_cohort_probs}, \code{indep_cohort_probs_overall}}{Same
#'       probabilities if an independent split was provided; otherwise \code{NA}.}
#'     \item{\code{sig_eps_sq}, \code{sig_eps_c_sq}}{Possibly estimated
#'       variance components carried forward.}
#'     \item{\code{lambda_ridge}}{Numeric scalar.  Value of the ridge penalty
#'       used (or \code{NA} if none).}
#'   }
#' @keywords internal
#' @noRd
prep_for_etwfe_regresion <- function(
	verbose,
	sig_eps_sq,
	sig_eps_c_sq,
	y,
	X_ints,
	X_mod,
	N,
	T,
	R,
	d,
	p,
	num_treats,
	add_ridge,
	first_inds,
	in_sample_counts,
	indep_count_data_available,
	indep_counts = NA,
	is_fetwfe = TRUE,
	is_twfe_covs = FALSE
) {
	if (verbose) {
		message("Getting omega sqrt inverse estimate...")
		t0 <- Sys.time()
	}

	if (is.na(sig_eps_sq) | is.na(sig_eps_c_sq)) {
		# Get omega_sqrt_inv matrix to multiply y and X_mod by on the left
		omega_res <- estOmegaSqrtInv(
			y,
			X_ints,
			N = N,
			T = T,
			p = p
		)

		sig_eps_sq <- omega_res$sig_eps_sq
		sig_eps_c_sq <- omega_res$sig_eps_c_sq

		rm(omega_res)

		if (verbose) {
			message("Done! Time to estimate noise variances:")
			message(Sys.time() - t0)
			t0 <- Sys.time()
		}
	}

	stopifnot(!is.na(sig_eps_sq) & !is.na(sig_eps_c_sq))

	Omega <- diag(rep(sig_eps_sq, T)) + matrix(sig_eps_c_sq, T, T)

	Omega_sqrt_inv <- expm::sqrtm(solve(Omega))

	if (verbose) {
		message("Time to get sqrt inverse matrix:")
		message(Sys.time() - t0)
	}

	y_final <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% y
	X_final <- kronecker(diag(N), sqrt(sig_eps_sq) * Omega_sqrt_inv) %*% X_mod

	stopifnot(ncol(X_final) == p)

	#
	#
	# twfeCovs modifications
	#
	#

	if (is_twfe_covs) {
		stopifnot(nrow(X_final) == N * T)
		# stopifnot(nrow(X_final_scaled) == N * T)

		# Drop columns corresponding to treatment effect interactions
		X_final <- X_final[, 1:(R + T - 1 + d * (1 + R + T - 1) + num_treats)]
		# X_final_scaled <- X_final_scaled[, 1:(R + T - 1 + d * (1 + R + T - 1) + num_treats)]

		# Drop columns corresponding to interactions between X and time, cohorts
		first_treat_ind <- R + T - 1 + d * (1 + R + T - 1)
		treat_inds <- first_treat_ind:(first_treat_ind + num_treats - 1)

		stopifnot(length(treat_inds) == num_treats)

		X_final <- X_final[, c(1:(R + T - 1 + d), treat_inds)]
		# X_final_scaled <- X_final_scaled[, c(1:(R + T - 1 + d), treat_inds)]

		stopifnot(nrow(X_final) == N * T)

		# Collapse together all columns corresponding to the same cohort
		# stopifnot(nrow(X_final_scaled) == N * T)
		treat_inds_mat <- matrix(as.numeric(NA), nrow = N * T, ncol = R)

		for (r in 1:R) {
			first_ind_r <- first_inds[r]
			if (r == R) {
				last_ind_r <- num_treats
			} else {
				last_ind_r <- first_inds[r + 1] - 1
			}

			cols_r <- R + T - 1 + d + first_ind_r:last_ind_r

			treat_inds_mat[, r] <- rowSums(X_final[, cols_r, drop = FALSE])
		}

		stopifnot(all(!is.na(treat_inds_mat)))

		X_final <- cbind(X_final[, 1:(R + T - 1 + d)], treat_inds_mat)
		# X_final_scaled <- cbind(X_final_scaled[, 1:(R + T - 1 + d)], treat_inds_mat)

		p_short <- R + T - 1 + d + R

		stopifnot(ncol(X_final) == p_short)
		# stopifnot(ncol(X_final_scaled) == p_short)

		treat_inds_short <- (R + T - 1 + d + 1):p_short

		stopifnot(length(treat_inds_short) == R)

		#
		#
		# Wrap up, store needed values
		#
		#

		p <- p_short
		treat_inds <- treat_inds_short
		num_treats <- R
	}

	#
	#
	# Optional: if using ridge regularization on untransformed coefficients,
	# add those rows now
	#
	#

	X_final_scaled <- my_scale(X_final)
	stopifnot(ncol(X_final_scaled) == p)
	scale_center <- attr(X_final_scaled, "scaled:center")
	scale_scale <- attr(X_final_scaled, "scaled:scale")

	if (add_ridge) {
		# Initialize identity matrix
		mat_to_multiply <- diag(p)

		if (is_fetwfe) {
			# First need to get D^{-1}:
			D_inverse <- genFullInvFusionTransformMat(
				first_inds = first_inds,
				T = T,
				R = R,
				d = d,
				num_treats = num_treats
			)

			stopifnot(ncol(D_inverse) == p)
			stopifnot(nrow(D_inverse) == p)

			mat_to_multiply <- D_inverse
		}

		# Now add rows to X_final_scaled
		lambda_ridge <- 0.00001 *
			(sig_eps_sq + sig_eps_c_sq) *
			sqrt(p / (N * T))

		X_final_scaled <- rbind(
			X_final_scaled,
			sqrt(lambda_ridge) * mat_to_multiply
		)
		y_final <- c(y_final, rep(0, p))

		stopifnot(length(y_final) == N * T + p)
		stopifnot(nrow(X_final_scaled) == N * T + p)
	} else {
		lambda_ridge <- as.numeric(NA)
	}

	#
	#
	# Step 2: get cohort-specific sample proportions (estimated treatment
	# probabilities)
	#
	#

	cohort_probs <- in_sample_counts[2:(R + 1)] /
		sum(in_sample_counts[2:(R + 1)])

	stopifnot(all(!is.na(cohort_probs)))
	stopifnot(all(cohort_probs >= 0))
	stopifnot(all(cohort_probs <= 1))
	stopifnot(length(cohort_probs) == R)
	stopifnot(abs(sum(cohort_probs) - 1) < 10^(-6))

	cohort_probs_overall <- in_sample_counts[2:(R + 1)] / N

	stopifnot(
		abs(1 - sum(cohort_probs_overall) - in_sample_counts[1] / N) < 10^(-6)
	)

	if (indep_count_data_available) {
		indep_cohort_probs <- indep_counts[2:(R + 1)] /
			sum(indep_counts[2:(R + 1)])

		stopifnot(all(!is.na(indep_cohort_probs)))
		stopifnot(all(indep_cohort_probs >= 0))
		stopifnot(all(indep_cohort_probs <= 1))
		stopifnot(length(indep_cohort_probs) == R)
		stopifnot(abs(sum(indep_cohort_probs) - 1) < 10^(-6))

		indep_cohort_probs_overall <- indep_counts[2:(R + 1)] / N

		stopifnot(
			abs(
				1 -
					sum(
						indep_cohort_probs_overall
					) -
					indep_counts[1] / N
			) <
				10^(-6)
		)
	} else {
		indep_cohort_probs <- NA
		indep_cohort_probs_overall <- NA
	}

	return(list(
		X_final_scaled = X_final_scaled,
		X_final = X_final,
		y_final = y_final,
		scale_center = scale_center,
		scale_scale = scale_scale,
		cohort_probs = cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs = indep_cohort_probs,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda_ridge = lambda_ridge
	))
}


# getGramInv
#' @title Compute Inverse of Gram Matrix for Selected Features
#' @description Calculates the inverse of the Gram matrix formed by the selected
#'   features from the (potentially transformed) design matrix `X_final`. This
#'   is a key component in calculating standard errors for the estimated
#'   coefficients.
#' @param N Integer; total number of units.
#' @param T Integer; total number of time periods.
#' @param X_final Numeric matrix; the final design matrix (e.g., after GLS and
#'   fusion transformations). Dimensions: `N*T` x `p_model`, where `p_model` is
#'   the total number of columns in this matrix.
#' @param sel_feat_inds Integer vector; indices of the features selected by the
#'   penalized regression, corresponding to columns in `X_final`.
#' @param treat_inds Integer vector; original indices (before selection) of the
#'   base treatment effect parameters. Used to subset `sel_feat_inds` to get
#'   only selected *treatment* features for the final Gram matrix.
#' @param num_treats Integer; total number of base treatment effect parameters.
#' @param sel_treat_inds_shifted Integer vector; indices of the selected
#'   treatment effects within the `num_treats` block, shifted to start from 1.
#'   Used to check dimensions.
#' @param calc_ses Logical; if `FALSE`, the function may return `NA` for
#'   `gram_inv`.
#' @return A list containing:
#'   \item{gram_inv}{The inverse of the Gram matrix corresponding to the
#'     *selected treatment effect features*. Returns `NA` if `calc_ses` is
#'     `FALSE` or if the Gram matrix is found to be singular.}
#'   \item{calc_ses}{Logical, potentially modified to `FALSE` if the Gram
#'     matrix is singular.}
#' @details
#'   1. Subsets `X_final` to include only columns specified by `sel_feat_inds`.
#'   2. Centers these selected columns.
#'   3. Computes the Gram matrix: `(1/(N*T)) * t(X_sel_centered) %*% X_sel_centered`.
#'   4. Checks if the minimum eigenvalue of the Gram matrix is too small (close to zero).
#'      If so, it issues a warning and sets `calc_ses` to `FALSE`, returning `NA`
#'      for `gram_inv`.
#'   5. Otherwise, it computes the inverse of the Gram matrix.
#'   6. It then subsets this inverse Gram matrix to retain only the rows/columns
#'      that correspond to the *selected treatment effects* (identified via
#'      `sel_feat_inds` and `treat_inds`).
#' @keywords internal
#' @noRd
getGramInv <- function(
	N,
	T,
	X_final,
	treat_inds,
	num_treats,
	calc_ses,
	sel_feat_inds = NA,
	sel_treat_inds_shifted = NA
) {
	stopifnot(nrow(X_final) == N * T)
	if (any(!is.na(sel_feat_inds))) {
		X_sel <- X_final[, sel_feat_inds, drop = FALSE]
		p_sel <- length(sel_feat_inds)
	} else {
		X_sel <- X_final
		p_sel <- ncol(X_final)
	}

	stopifnot(length(treat_inds) == num_treats)

	# Centering X_sel even when estimating via OLS is harmless
	X_sel_centered <- scale(X_sel, center = TRUE, scale = FALSE)

	gram <- 1 / (N * T) * (t(X_sel_centered) %*% X_sel_centered)

	stopifnot(nrow(gram) == p_sel)
	stopifnot(ncol(gram) == p_sel)

	min_gram_eigen <- min(
		eigen(gram, symmetric = TRUE, only.values = TRUE)$values
	)

	if (min_gram_eigen < 10^(-16)) {
		warning(
			"Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated."
		)
		return(list(gram_inv = NA, calc_ses = FALSE))
	}

	gram_inv <- solve(gram)

	if (any(!is.na(sel_feat_inds))) {
		# Get only the parts of gram_inv that have to do with treatment effects
		sel_treat_inds <- sel_feat_inds %in% treat_inds
		stopifnot(sum(sel_treat_inds) <= length(sel_feat_inds))
		stopifnot(length(sel_treat_inds) == length(sel_feat_inds))
	} else {
		sel_treat_inds <- rep(FALSE, p_sel)
		sel_treat_inds[treat_inds] <- TRUE
	}

	stopifnot(is.logical(sel_treat_inds))
	stopifnot(sum(sel_treat_inds) <= num_treats)
	stopifnot(length(sel_treat_inds) == p_sel)
	stopifnot(nrow(gram_inv) == ncol(gram_inv))

	stopifnot(all(!is.na(gram_inv)))

	gram_inv <- gram_inv[sel_treat_inds, sel_treat_inds]

	stopifnot(nrow(gram_inv) <= num_treats)
	stopifnot(nrow(gram_inv) == ncol(gram_inv))

	if (any(!is.na(sel_treat_inds_shifted))) {
		stopifnot(nrow(gram_inv) == length(sel_treat_inds_shifted))
		if (any(!is.na(sel_feat_inds))) {
			stopifnot(nrow(gram_inv) <= length(sel_feat_inds))
		}
	}

	return(list(gram_inv = gram_inv, calc_ses = calc_ses))
}


# estOmegaSqrtInv
#' @title Estimate Noise Variance Components and Omega Matrix
#' @description Estimates the idiosyncratic error variance (`sig_eps_sq`) and
#'   the unit-level random effect variance (`sig_eps_c_sq`) using the method
#'   described by Pesaran (2015, Section 26.5.1) with ridge regression.
#'   This function is called when these variances are not provided by the user.
#' @param y Numeric vector; the observed response variable, length `N*T`.
#' @param X_ints Numeric matrix; the design matrix including all fixed effects,
#'   covariates, treatment dummies, and interactions. `N*T` rows.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @param p Integer; the number of columns in `X_ints`.
#' @return A list containing two named elements:
#'   \item{sig_eps_sq}{Estimated variance of the idiosyncratic error term.}
#'   \item{sig_eps_c_sq}{Estimated variance of the unit-level random effect.}
#' @details The function first demeans `y` and `X_ints` within each unit (fixed
#'   effects transformation). Then, it fits a ridge regression (`alpha=0` in
#'   `glmnet::cv.glmnet`) to the demeaned data. Residuals from this model are
#'   used to estimate `sigma_hat_sq` (which corresponds to `sig_eps_sq`).
#'   The estimated unit-specific intercepts (`alpha_hat`) are then used to
#'   estimate `sigma_c_sq_hat` (corresponding to `sig_eps_c_sq`).
#' @references Pesaran, M. H. (2015). Time Series and Panel Data Econometrics.
#'   Oxford University Press.
#' @keywords internal
#' @noRd
estOmegaSqrtInv <- function(y, X_ints, N, T, p) {
	if (N * (T - 1) - p <= 0) {
		stop("Not enough units available to estimate the noise variance.")
	}
	stopifnot(N > 1)

	# Estimate standard deviations
	y_fe <- y
	X_ints_fe <- X_ints
	for (i in 1:N) {
		inds_i <- ((i - 1) * T + 1):(i * T)
		stopifnot(length(inds_i) == T)
		y_fe[inds_i] <- y[inds_i] - mean(y[inds_i])
		X_ints[inds_i, ] <- X_ints[inds_i, ] - colMeans(X_ints[inds_i, ])
	}

	lin_mod_fe <- glmnet::cv.glmnet(x = X_ints_fe, y = y_fe, alpha = 0)

	# Get residuals
	y_hat_fe <- stats::predict(lin_mod_fe, s = "lambda.min", newx = X_ints_fe)

	# Get coefficients
	beta_hat_fe <- stats::coef(lin_mod_fe, s = "lambda.min")[2:(p + 1)]

	resids <- y_fe - y_hat_fe

	tau <- rep(1, T)
	M <- diag(rep(1, T)) - outer(tau, tau) / T

	sigma_hat_sq <- 0
	alpha_hat <- rep(as.numeric(NA), N)

	for (i in 1:N) {
		inds_i <- ((i - 1) * T + 1):(i * T)
		sigma_hat_sq <- sigma_hat_sq +
			as.numeric(resids[inds_i] %*% M %*% resids[inds_i])
		alpha_hat[i] <- mean(y_hat_fe[inds_i]) -
			colMeans(X_ints_fe[inds_i, ] %*% beta_hat_fe)
	}

	stopifnot(all(!is.na(alpha_hat)))

	sigma_hat_sq <- sigma_hat_sq / (N * (T - 1) - p)

	sigma_c_sq_hat <- sum((alpha_hat - mean(alpha_hat))^2) / (N - 1)

	return(list(sig_eps_sq = sigma_hat_sq, sig_eps_c_sq = sigma_c_sq_hat))
}

#' @title Validate Inputs for the ETWFE Core Estimator
#'
#' @description
#' Performs a complete set of argument, dimension, and consistency checks
#' for \code{\link{etwfe_core}}.
#' The function halts with informative error messages whenever a necessary
#' condition is violated (e.g., negative variances, inconsistent cohort
#' counts, or ill-posed significance levels).
#' When all checks pass it returns a compact list of derived quantities
#' used downstream by the core estimator.
#'
#' @param in_sample_counts Integer named vector. Length `R+1`.
#'   The first element must be the number of never-treated units; the
#'   remaining `R` elements give the number of units in each treated
#'   cohort.  Names must be unique and correspond to cohort identifiers.
#' @param N Integer. Total number of unique units in the (filtered) data.
#' @param T Integer. Total number of time periods.
#' @param sig_eps_sq Numeric scalar or \code{NA}.  Prespecified observation-level
#'   variance component (see Section 2 of Faletto 2025).  If \code{NA},
#'   it will later be estimated.
#' @param sig_eps_c_sq Numeric scalar or \code{NA}.  Prespecified unit-level
#'   random-effect variance component.  If \code{NA}, it will later be estimated.
#' @param indep_counts Optional integer vector.  Cohort counts from an
#'   independent sample used to obtain asymptotically exact SEs for the ATT.
#'   Must have the same length and ordering as \code{in_sample_counts}.
#'   Default is \code{NA}.
#' @param verbose Logical.  Forwarded verbosity flag.  Default \code{FALSE}.
#' @param alpha Numeric in \eqn{(0,1)}.  Significance level requested for
#'   confidence intervals.  Default \code{0.05}.
#' @param add_ridge Logical.  Whether a small ridge penalty will later be
#'   applied inside \code{etwfe_core}.  Used only to customize warning
#'   messages.  Default \code{FALSE}.
#'
#' @details
#' Key validations include:
#' \itemize{
#'   \item Ensuring the total of \code{in_sample_counts} equals \code{N}.
#'   \item Checking that at least one never-treated unit is present.
#'   \item Verifying that cohort names are unique and counts non-negative.
#'   \item Confirming that \eqn{1 \le R \le T-1}.
#'   \item Basic sanity checks for all numeric scalars (\code{alpha} inside \eqn{(0,1)}, non-negative variances, etc.).
#'   \item All structural requirements on \code{indep_counts} when supplied.
#' }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{R}}{Integer. The number of treated cohorts `=length(in_sample_counts)-1`).}
#'     \item{\code{c_names}}{Character vector of cohort names (length `R`).}
#'     \item{\code{indep_count_data_available}}{Logical. \code{TRUE} if
#'       valid \code{indep_counts} were supplied, \code{FALSE} otherwise.}
#'   }
#' @keywords internal
#' @noRd
check_etwfe_core_inputs <- function(
	in_sample_counts,
	N,
	T,
	sig_eps_sq,
	sig_eps_c_sq,
	indep_counts,
	verbose,
	alpha,
	add_ridge
) {
	R <- length(in_sample_counts) - 1

	c_names <- names(in_sample_counts)[2:(R + 1)]

	stopifnot(N >= 2) # bare minimum, 2 units at 2 times

	stopifnot(T >= 2) # bare minimum, 2 units at 2 times

	if (any(!is.na(sig_eps_sq))) {
		stopifnot(is.numeric(sig_eps_sq) | is.integer(sig_eps_sq))
		stopifnot(length(sig_eps_sq) == 1)
		stopifnot(sig_eps_sq >= 0)
	}

	if (any(!is.na(sig_eps_c_sq))) {
		stopifnot(is.numeric(sig_eps_c_sq) | is.integer(sig_eps_c_sq))
		stopifnot(length(sig_eps_c_sq) == 1)
		stopifnot(sig_eps_c_sq >= 0)
	}

	stopifnot(sum(in_sample_counts) == N)
	stopifnot(all(in_sample_counts >= 0))
	if (in_sample_counts[1] == 0) {
		stop(
			"No never-treated units detected in data to fit model; estimating treatment effects is not possible"
		)
	}
	if (length(names(in_sample_counts)) != length(in_sample_counts)) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}

	if (
		length(names(in_sample_counts)) !=
			length(unique(names(in_sample_counts)))
	) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)"
		)
	}

	stopifnot(R >= 1)
	stopifnot(R <= T - 1)

	indep_count_data_available <- FALSE
	if (any(!is.na(indep_counts))) {
		if (sum(indep_counts) != N) {
			stop(
				"Number of units in independent cohort count data does not equal number of units in data to be used to fit model."
			)
		}
		if (length(indep_counts) != length(in_sample_counts)) {
			stop(
				"Number of counts in independent counts does not match number of cohorts in data to be used to fit model."
			)
		}
		if (any(indep_counts <= 0)) {
			stop(
				"At least one cohort in the independent count data has 0 members"
			)
		}
		indep_count_data_available <- TRUE
	}

	stopifnot(is.logical(verbose))
	stopifnot(length(verbose) == 1)

	stopifnot(is.numeric(alpha))
	stopifnot(length(alpha) == 1)
	stopifnot(alpha > 0)
	stopifnot(alpha < 1)

	stopifnot(is.logical(add_ridge))
	stopifnot(length(add_ridge) == 1)

	return(list(
		R = R,
		c_names = c_names,
		indep_count_data_available = indep_count_data_available
	))
}

#' Build the **entire** inverse-fusion matrix \(D_N^{-1}\)
#'
#' Constructs the \eqn{p\times p} block–diagonal matrix that sends the sparse
#' coefficient vector \(\theta\) used in the bridge loss back to the original
#' coefficient scale \(\beta\), and is also used to multiply \(Z\) on the right
#' to transform the features into a transformed feature space:
#' \deqn{\beta \;=\;D_N^{-1}\,\theta.}
#' The block layout follows the factorisation proved in the
#' paper – repeated here using the helper generators already available in the
#' package.
#'
#' \preformatted{
#'       D_N^{-1} = diag(
#'         (D^{(1)}(R))^{-1},                          # 1. cohort FEs
#'         (D^{(1)}(T-1))^{-1},                        # 2. time  FEs
#'         I_d,                                        # 3. X main effects
#'         I_d ⊗ (D^{(1)}(R))^{-1},                   # 4. cohort × X
#'         I_d ⊗ (D^{(1)}(T-1))^{-1},                 # 5. time   × X
#'         (D^{(2)}(𝓡))^{-1},                         # 6. base τ_{r,t}
#'         I_d ⊗ (D^{(2)}(𝓡))^{-1} )                  # 7. τ_{r,t} × X
#' }
#'
#' @section Block dimensions:
#' \itemize{
#'   \item Cohort FEs: \(R\times R\)
#'   \item Time-period FEs: \((T-1)\times(T-1)\)
#'   \item Identitites: \(d\), \(dR\), \(d(T-1)\)
#'   \item Treatment blocks: \( \mathfrak W \times \mathfrak W\) with
#'         \(\mathfrak W = \texttt{num_treats}\)
#' }
#' All non–identity blocks contain only 0/1 entries, so the determinant of
#' the whole matrix is 1 (volume–preserving transform).
#'
#' @param first_inds Integer vector (length \code{R}).
#'   `first_inds[r]` is the 1-based column index of the first base
#'   treatment-effect parameter \(\tau_{r,0}\) for cohort \code{r}.
#' @param T Integer. Total number of time periods \(\ge 3\).
#' @param R Integer. Number of treated cohorts (\(\ge 1\)).
#'   The function stops if you accidentally pass \code{R = 0}.
#' @param d Integer. Number of time-invariant covariates (can be 0).
#' @param num_treats Integer.  Total number of base treatment-effect
#'   coefficients \(\mathfrak W = T R - R(R+1)/2\).
#'
#' @return A dense base-R matrix of size
#'   \eqn{p \times p} with
#'   \eqn{p = R + (T-1) + d + dR + d(T-1) + \mathfrak W + d\mathfrak W}.
#'
#' @examples
#' R  <- 3; T <- 6; d <- 2
#' nt <- getNumTreats(R, T)
#' Dinv <- genFullInvFusionTransformMat(getFirstInds(R,T), T, R, d, nt)
#' dim(Dinv)   # should be p × p
#'
#' @keywords internal
#' @noRd
genFullInvFusionTransformMat <- function(first_inds, T, R, d, num_treats) {
	##———— Safety checks ————————————————————————————————————————————
	stopifnot(is.numeric(R), length(R) == 1L, R >= 1L)
	stopifnot(is.numeric(T), length(T) == 1L, T >= 3L, R <= T - 1)
	stopifnot(is.numeric(d), length(d) == 1L, d >= 0L)
	stopifnot(length(first_inds) == R)

	##———— 1. Cohort fixed-effects block:  (D^{(1)}(R))^{-1} ————————
	block1 <- genBackwardsInvFusionTransformMat(R)

	##———— 2. Time fixed-effects block:    (D^{(1)}(T-1))^{-1} ———————
	block2 <- genBackwardsInvFusionTransformMat(T - 1)

	##———— 3. Covariate main effects:      I_d ————————————————
	block3 <- if (d > 0) diag(d) else NULL

	##———— 4. Cohort × X interactions:     I_d ⊗ (D^{(1)}(R))^{-1} ——
	block4 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(R)) else NULL

	##———— 5. Time × X interactions:       I_d ⊗ (D^{(1)}(T-1))^{-1} —
	block5 <- if (d > 0)
		kronecker(diag(d), genBackwardsInvFusionTransformMat(T - 1)) else NULL

	##———— 6. Base treatment effects:      (D^{(2)}(𝓡))^{-1} ————————
	block6 <- genInvTwoWayFusionTransformMat(
		n_vars = num_treats,
		first_inds = first_inds,
		R = R
	)

	##———— 7. Treatment × X interactions:  I_d ⊗ (D^{(2)}(𝓡))^{-1} ——
	block7 <- if (d > 0)
		kronecker(
			diag(d),
			genInvTwoWayFusionTransformMat(
				n_vars = num_treats,
				first_inds = first_inds,
				R = R
			)
		) else NULL

	## Gather present blocks in the same order as the theoretical expression
	blocks <- list(block1, block2, block3, block4, block5, block6, block7)
	blocks <- Filter(Negate(is.null), blocks) # drop NULLs for d = 0

	##———— Assemble block-diagonal matrix ————————————————
	## Matrix::bdiag() returns a sparse dgCMatrix.  We convert to base-R matrix
	## here because downstream (ridge-row augmentation) works with dense objects
	## via rbind().
	full_D_inv <- as.matrix(Matrix::bdiag(blocks))

	##———— Dimension cross-check ————————————————————————————
	p <- getP(R = R, T = T, d = d, num_treats = num_treats)
	stopifnot(nrow(full_D_inv) == p, ncol(full_D_inv) == p)

	return(full_D_inv)
}
