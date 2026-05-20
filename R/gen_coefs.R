# Coefficient generators and treatment-effect truth extraction for the
# simulation pipeline: `genCoefs()`, its core (`genCoefsCore()`), and
# `getTes()` / `getActualCohortTes()`. Moved from R/gen_funcs.R in
# 1.9.25.

#' Generate Coefficient Vector for Data Generation
#'
#' This function generates a coefficient vector \code{beta} for simulation studies of the fused
#' extended two-way fixed effects estimator. It returns an S3 object of class
#' \code{"FETWFE_coefs"} containing \code{beta} along with simulation parameters \code{R},
#' \code{T}, and \code{d}. See the simulation studies section of Faletto (2025) for details.
#'
#' @param R Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to
#'   \code{R + 1}).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms
#'   corresponding to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1). The probability that any given entry in the initial sparse
#'   coefficient vector \code{theta} is nonzero.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each
#'   nonzero entry is set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a
#'   positive value).
#' @param seed (Optional) Integer. Seed for reproducibility.
#'
#' @return An object of class \code{"FETWFE_coefs"}, which is a list containing:
#' \describe{
#'   \item{beta}{A numeric vector representing the full coefficient vector after the inverse fusion
#'      transform.}
#'   \item{theta}{A numeric vector representing the coefficient vector in the transformed feature
#'		space. \code{theta} is a sparse vector, which aligns with an assumption that deviations from the
#'		restrictions encoded in the FETWFE model are sparse. \code{beta} is derived from
#'		\code{theta}.}
#'   \item{R}{The provided number of treated cohorts.}
#'   \item{T}{The provided number of time periods.}
#'   \item{d}{The provided number of covariates.}
#'   \item{seed}{The provided seed.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = R + (T - 1) + d + dR + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = R + (T - 1) + d + dR + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times R - \frac{R(R+1)}{2}}{num_treats = T * R - R(R+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries
#'   occurring with probability \code{density}. Nonzero entries are set to \code{eff_size} or
#'   \code{-eff_size} (with a 60\% chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion
#'   transform to \code{theta} using internal routines (e.g.,
#'   \code{genBackwardsInvFusionTransformMat()} and \code{genInvTwoWayFusionTransformMat()}).
#' }
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#' }
#'
#' @export
genCoefs <- function(R, T, d, density, eff_size, seed = NULL) {
	# Check that T is a numeric scalar and at least 3.
	if (!is.numeric(T) || length(T) != 1 || T < 3) {
		stop("T must be a numeric value greater than or equal to 3")
	}

	# Check that R is a numeric scalar and at least 2.
	if (!is.numeric(R) || length(R) != 1 || R < 2) {
		stop(
			"R must be a numeric value greater than or equal to 2 (currently there is only support for data sets with staggered adoptions, so at least two treated cohorts)"
		)
	}

	# Check that R does not exceed T - 1.
	if (R > T - 1) {
		stop("R must be less than or equal to T - 1")
	}

	# Check that d is a numeric scalar and is non-negative.
	if (!is.numeric(d) || length(d) != 1 || d < 0) {
		stop("d must be a non-negative numeric value")
	}

	# Check that density is a numeric scalar strictly between 0 and 1.
	if (
		!is.numeric(density) ||
			length(density) != 1 ||
			density <= 0 ||
			density >= 1
	) {
		stop("density must be numeric and strictly between 0 and 1")
	}

	# Check that eff_size is numeric.
	if (!is.numeric(eff_size) || length(eff_size) != 1) {
		stop("eff_size must be a numeric value")
	}

	stopifnot(R >= 2)
	stopifnot(T >= 3)
	stopifnot(R <= T - 1)

	core_obj <- genCoefsCore(
		R = R,
		T = T,
		d = d,
		density = density,
		eff_size = eff_size,
		seed = seed
	)
	if (is.null(core_obj$beta)) {
		stop(
			"Internal error: genCoefsCore() did not return expected components."
		)
	}

	# Create an S3 object of class "FETWFE_coefs"
	obj <- list(
		beta = core_obj$beta,
		theta = core_obj$theta,
		R = R,
		T = T,
		d = d,
		seed = seed
	)
	class(obj) <- "FETWFE_coefs"
	return(obj)
}


#' Compute True Treatment Effects
#'
#' @description
#' This function extracts the true treatment effects from a full coefficient vector
#' as generated by \code{genCoefs()}. It calculates the overall average treatment effect on the
#' treated (ATT) as the equal-weighted average of the cohort-specific treatment effects, and also
#' returns the individual treatment effects for each treated cohort.
#'
#' @param coefs_obj An object of class \code{"FETWFE_coefs"} containing the coefficient vector
#' and simulation parameters.
#'
#' @return An object of class \code{"FETWFE_tes"}, which is a list with the
#' following elements:
#' \describe{
#'   \item{att_true}{A numeric value representing the overall average treatment
#'         effect on the treated. It is computed as the (equal-weighted) mean of
#'         the cohort-specific treatment effects.}
#'   \item{actual_cohort_tes}{A numeric vector of length \code{R} containing the
#'         true cohort-specific treatment effects, calculated by averaging the
#'         coefficients corresponding to the treatment dummies for each cohort.}
#'   \item{R, T, d, seed}{The generating parameters carried over from
#'         \code{coefs_obj} so that \code{print()} and \code{summary()} on the
#'         returned object are self-describing.}
#' }
#' Use \code{print()} or \code{summary()} on the returned object for a
#' formatted display.
#'
#' @details
#' The function internally uses auxiliary routines \code{getNumTreats()}, \code{getP()},
#' \code{getFirstInds()}, \code{getTreatInds()}, and \code{getActualCohortTes()} to determine the
#' correct indices of treatment effect coefficients in \code{beta}. The overall treatment effect
#' is computed as the simple average of these cohort-specific effects.
#'
#' @examples
#' \dontrun{
#' # Generate coefficients
#' coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#' # Compute the true treatment effects:
#' te_results <- getTes(coefs)
#'
#' # Overall average treatment effect on the treated:
#' print(te_results$att_true)
#'
#' # Cohort-specific treatment effects:
#' print(te_results$actual_cohort_tes)
#'
#' # Or use the new print method for a self-describing display:
#' print(te_results)
#' }
#'
#' @export
getTes <- function(coefs_obj) {
	if (!inherits(coefs_obj, "FETWFE_coefs")) {
		stop("coefs_obj must be an object of class 'FETWFE_coefs'")
	}

	# Unpack components from the coefs object
	beta <- coefs_obj$beta
	R <- coefs_obj$R
	T <- coefs_obj$T
	d <- coefs_obj$d

	num_treats <- getNumTreats(R = R, T = T)

	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(length(beta) == p)

	first_inds <- getFirstInds(R = R, T = T)
	treat_inds <- getTreatInds(R = R, T = T, d = d, num_treats = num_treats)

	actual_cohort_tes <- getActualCohortTes(
		R = R,
		first_inds = first_inds,
		treat_inds = treat_inds,
		coefs = beta,
		num_treats = num_treats
	)

	att_true <- as.numeric(mean(actual_cohort_tes))

	# Cohort adoption times in the simulator's convention: cohort r adopts
	# at calendar time r + 1 (cohort 0 = never-treated, by convention
	# encoded in the panel's `time` integer values 1..T). Stored so
	# downstream tooling (e.g., `tidy.FETWFE_tes`) can label rows with the
	# same scheme that `tidy.<estimator>` uses on a fitted panel.
	cohort_times <- as.integer(seq_len(R) + 1L)

	out <- list(
		att_true = att_true,
		actual_cohort_tes = actual_cohort_tes,
		R = R,
		T = T,
		d = d,
		seed = coefs_obj$seed,
		cohort_times = cohort_times
	)
	class(out) <- "FETWFE_tes"
	return(out)
}


#' Generate Coefficient Vector for Data Generation
#'
#' This function generates a coefficient vector \code{beta} along with a sparse auxiliary vector
#' \code{theta} for simulation studies of the fused extended two-way fixed effects estimator. The
#' returned \code{beta} is formatted to align with the design matrix created by
#' \code{genRandomData()}, and is a valid input for the \code{beta} argument of that function. The
#' vector \code{theta} is sparse, with nonzero entries occurring with probability \code{density} and
#' scaled by \code{eff_size}. See the simulation studies section of Faletto (2025) for details.
#'
#' @param R Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to
#' \code{R + 1}).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms
#' corresponding to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1). The probability that any given entry in the initial sparse
#' coefficient vector \code{theta} is nonzero.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each
#' nonzero entry is set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a
#' positive value).
#' @param seed (Optional) Integer. Seed for reproducibility.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{beta}}{A numeric vector representing the full coefficient vector after the inverse
#'   fusion transform.}
#'   \item{theta}{A numeric vector representing the coefficient vector in the transformed feature
#'		space. \code{theta} is a sparse vector, which aligns with an assumption that deviations from the
#'		restrictions encoded in the FETWFE model are sparse. \code{beta} is derived from
#'		\code{theta}.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = R + (T - 1) + d + dR + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = R + (T - 1) + d + dR + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times R - \frac{R(R+1)}{2}}{num_treats = T * R - R(R+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries
#' occurring
#'   with probability \code{density}. Nonzero entries are set to \code{eff_size} or \code{-eff_size}
#'   (with a 60\% chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion
#'   transform to \code{theta} using internal routines (e.g.,
#'   \code{genBackwardsInvFusionTransformMat()} and \code{genInvTwoWayFusionTransformMat()}).
#' }
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' \dontrun{
#'   # Set parameters for the coefficient generation
#'   R <- 3         # Number of treated cohorts
#'   T <- 6         # Total number of time periods
#'   d <- 2         # Number of covariates
#'   density <- 0.1 # Probability that an entry in the initial vector is nonzero
#'   eff_size <- 1.5  # Scaling factor for nonzero coefficients
#'   seed <- 789    # Seed for reproducibility
#'
#'   # Generate coefficients using genCoefsCore()
#'   coefs_core <- genCoefsCore(R = R, T = T, d = d, density = density,
#'   eff_size = eff_size, seed = seed)
#'   beta <- coefs_core$beta
#'   theta <- coefs_core$theta
#'
#'   # For diagnostic purposes, compute the expected length of beta.
#'   # The length p is defined internally as:
#'   #   p = R + (T - 1) + d + d*R + d*(T - 1) + num_treats + num_treats*d,
#'   # where num_treats = T * R - (R*(R+1))/2.
#'   num_treats <- T * R - (R * (R + 1)) / 2
#'   p_expected <- R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
#'
#'   cat("Length of beta:", length(beta), "\nExpected length:", p_expected, "\n")
#' }
#'
#' @export
genCoefsCore <- function(R, T, d, density, eff_size, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	# Check that T is a numeric scalar and at least 3.
	if (!is.numeric(T) || length(T) != 1 || T < 3) {
		stop("T must be a numeric value greater than or equal to 3")
	}

	# Check that R is a numeric scalar and at least 2.
	if (!is.numeric(R) || length(R) != 1 || R < 2) {
		stop(
			"R must be a numeric value greater than or equal to 2 (currently there is only support for data sets with staggered adoptions, so at least two treated cohorts)"
		)
	}

	# Check that R does not exceed T - 1.
	if (R > T - 1) {
		stop("R must be less than or equal to T - 1")
	}

	# Check that d is a numeric scalar and is non-negative.
	if (!is.numeric(d) || length(d) != 1 || d < 0) {
		stop("d must be a non-negative numeric value")
	}

	# Check that density is a numeric scalar strictly between 0 and 1.
	if (
		!is.numeric(density) ||
			length(density) != 1 ||
			density <= 0 ||
			density >= 1
	) {
		stop("density must be numeric and strictly between 0 and 1")
	}

	# Check that eff_size is numeric.
	if (!is.numeric(eff_size) || length(eff_size) != 1) {
		stop("eff_size must be a numeric value")
	}

	stopifnot(R >= 2)
	stopifnot(T >= 3)
	stopifnot(R <= T - 1)

	num_treats <- getNumTreats(R = R, T = T)

	p <- getP(R = R, T = T, d = d, num_treats = num_treats)

	theta <- rep(0, p)

	# Make sure at least one feature is selected
	pass_condition <- FALSE
	while (!pass_condition) {
		theta_inds <- which(as.logical(rbinom(n = p, size = 1, prob = density)))
		pass_condition <- length(theta_inds) > 0
	}

	num_coefs <- length(theta_inds)
	# Generate signs of coefficients in transformed space, and bias away from
	# 0.5 (as described in paper)
	signs <- rfunc(num_coefs, prob = 0.6)

	theta[theta_inds] <- eff_size * signs

	# Now we have coefficients that are sparse in the appropriate feature space.
	# The last step is to transform them to the original feature space. Since
	# theta = D %*% beta, beta = solve(D) %*% theta.
	beta <- rep(as.numeric(NA), p)

	beta[1:R] <- genBackwardsInvFusionTransformMat(R) %*% theta[1:R]

	stopifnot(all(is.na(beta[(R + 1):(R + T - 1)])))
	beta[(R + 1):(R + T - 1)] <- genBackwardsInvFusionTransformMat(T - 1) %*%
		theta[(R + 1):(R + T - 1)]

	if (d > 0) {
		# Coefficients corresponding to X don't need to be transformed
		stopifnot(all(is.na(beta[(R + T - 1 + 1):(R + T - 1 + d)])))
		beta[(R + T - 1 + 1):(R + T - 1 + d)] <- theta[
			(R + T - 1 + 1):(R + T - 1 + d)
		]

		# Cohort-X interactions (one cohort at a time, with all interactions for
		# X. So R blocks of size d.)

		for (j in 1:d) {
			first_ind_j <- R + T - 1 + d + j
			last_ind_j <- R + T - 1 + d + (R - 1) * d + j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)

			stopifnot(length(inds_j) == R)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- genBackwardsInvFusionTransformMat(R) %*%
				theta[inds_j]
		}

		stopifnot(all(!is.na(beta[1:(R + T - 1 + d + R * d)])))
		stopifnot(all(is.na(beta[(R + T - 1 + d + R * d + 1):p])))

		# Time-X interactions
		for (j in 1:d) {
			first_ind_j <- R + T - 1 + d + R * d + j
			last_ind_j <- R + T - 1 + d + R * d + (T - 2) * d + j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)
			stopifnot(length(inds_j) == T - 1)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- genBackwardsInvFusionTransformMat(T - 1) %*%
				theta[inds_j]
		}

		stopifnot(all(!is.na(beta[1:(R + T - 1 + d + R * d + (T - 1) * d)])))
		stopifnot(all(is.na(beta[(R + T - 1 + d + R * d + (T - 1) * d + 1):p])))
	}

	# Base treatment effects: need to identify indices of first treatment
	# effect for each cohort
	first_inds <- getFirstInds(R = R, T = T)

	treat_inds <- getTreatInds(R = R, T = T, d = d, num_treats = num_treats)

	stopifnot(all(is.na(beta[treat_inds])))

	beta[treat_inds] <- genInvTwoWayFusionTransformMat(
		num_treats,
		first_inds,
		R
	) %*%
		theta[treat_inds]

	stopifnot(all(
		!is.na(beta[1:(R + T - 1 + d + R * d + (T - 1) * d + num_treats)])
	))

	if (d > 0) {
		stopifnot(all(is.na(beta[
			(R + T - 1 + d + R * d + (T - 1) * d + num_treats + 1):p
		])))

		# Treatment effect-X interactions
		for (j in 1:d) {
			first_ind_j <- R + T - 1 + d + R * d + (T - 1) * d + num_treats + j
			last_ind_j <- R +
				T -
				1 +
				d +
				R * d +
				(T - 1) * d +
				num_treats +
				(num_treats - 1) * d +
				j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)

			stopifnot(length(inds_j) == num_treats)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- genInvTwoWayFusionTransformMat(
				num_treats,
				first_inds,
				R
			) %*%
				theta[inds_j]
		}
	}

	stopifnot(all(!is.na(beta)))

	# Confirm beta satisfies input requirements of genRandomData() (make
	# up values for N, sig_eps_sq, and sig_eps_c_sq that meet requirements)

	testGenRandomDataInputs(
		beta = beta,
		R = R,
		T = T,
		d = d,
		N = R + 1,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	return(list(beta = beta, theta = theta))
}


#' Calculate True Cohort Average Treatment Effects from Coefficients
#'
#' Given a full coefficient vector and information about treatment effect indices,
#' this function calculates the true average treatment effect for each cohort by
#' averaging the relevant treatment effect coefficients.
#'
#' @param R Integer. Number of treated cohorts.
#' @param first_inds Integer vector. `first_inds[r]` is the index (within the
#'   block of treatment effect coefficients) of the first treatment effect for cohort `r`.
#' @param treat_inds Integer vector. Indices in the full `coefs` vector that
#'   correspond to the block of all treatment effect coefficients.
#' @param coefs Numeric vector. The full true coefficient vector \eqn{\beta}.
#' @param num_treats Integer. Total number of treatment effect parameters.
#'
#' @return A numeric vector of length R, where the r-th element is the true
#'   average treatment effect for cohort r.
#' @keywords internal
#' @noRd
getActualCohortTes <- function(R, first_inds, treat_inds, coefs, num_treats) {
	actual_cohort_tes <- rep(as.numeric(NA), R)

	for (r in 1:R) {
		inds_r <- .cohort_block_inds(r, R, first_inds, num_treats)
		actual_cohort_tes[r] <- mean(coefs[treat_inds][inds_r])
	}

	return(actual_cohort_tes)
}
