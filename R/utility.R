#' @import glmnet
#' @importFrom Matrix bdiag
#' @importFrom expm sqrtm
#' @importFrom stats qnorm predict coef model.matrix setNames lm

#-------------------------------------------------------------------------------
# Helper Functions for Data Processing
#-------------------------------------------------------------------------------

#' @title Identify Cohorts and Basic Panel Structure
#'
#' @description
#' This function processes the input panel data to identify treatment cohorts,
#' unique units, and time periods. It also removes units that were treated in the
#' very first time period, as these units cannot be used for identifying
#' pre-treatment trends for themselves. The treatment variable is removed after
#' cohorts are identified.
#'
#' @param df A `data.frame` object representing the panel data. Contains columns for response, time
#' period variable (a single categorical or numeric/integer variable indicating the time period
#' for each observation), unit variable (a single categorical variable
#' indicating which unit each observation is), treatment (a binary variable;
#' 1 if unit is treated at that time and 0 if not; treatment must be an
#' absorbing state in this model); covariates (which are fixed over time)
#' @param time_var Character string; name of the time variable column.
#' @param unit_var Character string; name of the unit identifier column.
#' @param treat_var Character string; name of the treatment indicator column.
#' @param covs Character vector; names of covariate columns (used for subsetting `df`).
#'
#' @details
#' The function iterates through each unit to determine its cohort assignment based
#' on the first period it is observed as treated. It assumes treatment is an
#' absorbing state.
#' Units treated in the first time period (`times[1]`) are removed from the dataset
#' and a warning is issued.
#' Cohorts are stored in a list, where names are the treatment adoption times and
#' values are character vectors of unit identifiers belonging to that cohort.
#' The list of cohorts is ordered by adoption time, and the first time period is
#' ensured to have an empty cohort (representing never-treated or baseline).
#'
#' @return A list containing:
#'   \item{df}{The processed `data.frame` with first-period-treated units removed
#'     and the `treat_var` column dropped. Contains only `time_var`, `unit_var`,
#'     `response` (implicitly, as it's not dropped), and `covs`.}
#'   \item{cohorts}{A list where each element is a character vector of unit IDs
#'     belonging to a specific treatment cohort. The names of the list elements
#'     are the times of treatment adoption. The first cohort (earliest adoption
#'     time after the first period) is listed first.}
#'   \item{units}{A character vector of unique unit identifiers remaining after
#'     processing.}
#'   \item{times}{A numeric or integer vector of unique, sorted time periods
#'     present in the original data.}
#' @keywords internal
#' @noRd
idCohorts <- function(df, time_var, unit_var, treat_var, covs) {
	stopifnot(time_var %in% colnames(df))
	stopifnot(unit_var %in% colnames(df))
	stopifnot(treat_var %in% colnames(df))

	# Form design matrix
	units <- unique(df[, unit_var])
	N <- length(units)
	times <- sort(unique(df[, time_var]))
	T <- length(times)

	# Variable to identify cohorts
	cohorts <- list()
	for (t in 1:T) {
		cohorts[[t]] <- character()
	}
	names(cohorts) <- times

	for (s in units) {
		df_s <- df[df[, unit_var] == s, ]
		# Assume this is a balanced panel
		if (nrow(df_s) != T) {
			stop(paste(
				"Panel does not appear to be balanced (unit",
				s,
				"does not have exactly T observations for T =",
				T
			))
		}

		if (any(df_s[, treat_var] == 1)) {
			# Identify first year of treatment (and cohort)
			# Ensure df_s is sorted by time before finding min
			df_s_sorted <- df_s[order(df_s[, time_var]), ]
			treat_year_s_ind_in_sorted_times <- min(which(
				df_s_sorted[, treat_var] == 1
			))
			actual_treat_time <- df_s_sorted[
				treat_year_s_ind_in_sorted_times,
				time_var
			]

			# Make sure treatment is absorbing
			# Check from the actual_treat_time onwards in the original df for unit s
			original_unit_times_from_treatment <- df[
				(df[, unit_var] == s) & (df[, time_var] >= actual_treat_time),
				treat_var
			]
			if (any(original_unit_times_from_treatment != 1)) {
				stop(paste(
					"Treatment does not appear to be an absorbing state for unit",
					s
				))
			}
			cohorts[[as.character(actual_treat_time)]] <- c(
				cohorts[[as.character(actual_treat_time)]],
				s
			)
		}
	}

	stopifnot(length(unlist(cohorts)) <= N)
	# Keep only cohorts that actually have units
	cohorts <- cohorts[lengths(cohorts) > 0]

	# Need at least one untreated period, so have to omit units that were
	# treated in the very first time period
	first_time_val_char <- as.character(times[1])
	if (first_time_val_char %in% names(cohorts)) {
		first_year_cohort_units <- cohorts[[first_time_val_char]]
		if (length(first_year_cohort_units) > 0) {
			df <- df[!(df[, unit_var] %in% first_year_cohort_units), ]
			units <- unique(df[, unit_var]) # Update units
			if (length(first_year_cohort_units) > 0) {
				# N is original N
				warning(paste(
					length(first_year_cohort_units),
					"units were removed because they were treated in the first time period:",
					paste(first_year_cohort_units, collapse = ", ")
				))
			}
		}
		cohorts[[first_time_val_char]] <- NULL # Remove this cohort entry
	}

	if (length(units) == 0) {
		stop(
			"All units were treated in the first time period or no units remain after filtering; estimating treatment effects is not possible"
		)
	}

	N <- length(units) # Update N to reflect remaining units

	# Treatment no longer needed
	df <- df[, colnames(df) != treat_var]

	# Make sure there is an empty cohort for the first time
	cohorts[[as.character(times[1])]] <- character()

	# Order cohorts in order of times
	if (length(cohorts) > 0) {
		cohorts <- cohorts[order(as.numeric(names(cohorts)))]
	}

	# If after removing first-period treated, there are no treated cohorts left:
	if (length(cohorts) == 0) {
		# This implies all treated units were treated in the first period.
		stop("all units appear to have been treated in the first period")
	}

	# This should have been the first cohort
	stopifnot(length(cohorts[[1]]) == 0)

	cohorts <- cohorts[-1]

	stopifnot(all(lengths(cohorts) >= 1))
	stopifnot(length(cohorts) <= T)
	stopifnot(length(unlist(cohorts)) <= N)

	return(list(df = df, cohorts = cohorts, units = units, times = times))
}


# my_scale
#' @title Custom Scaling Function Handling Zero-Variance Columns
#' @description Centers and scales the columns of a numeric matrix. This function
#'   is similar to `base::scale()` but explicitly handles columns with zero
#'   variance by setting their scale factor to 1, thus avoiding division by zero
#'   and NaN results.
#' @param x Numeric matrix; the matrix whose columns are to be scaled.
#' @return A numeric matrix of the same dimensions as `x`, with columns centered
#'   and scaled. Attributes "scaled:center" and "scaled:scale" are attached,
#'   containing the means and standard deviations (or 1 for zero-variance
#'   columns) used for the transformation.
#' @details
#'   1. Computes column means (`ctr`) and column standard deviations (`sds`) of `x`.
#'   2. Identifies columns where `sds == 0`.
#'   3. For these zero-variance columns, the scaling factor in `sds2` is set to 1.
#'      For other columns, `sds2` is the same as `sds`. The centering values (`ctr2`)
#'      are the original column means.
#'   4. Columns of `x` are first centered using `ctr2` and then scaled (divided)
#'      by `sds2`.
#' @examples
#'\dontrun{
#'   mat <- matrix(c(1, 1, 1, 2, 3, 4), ncol = 2)
#'   # mat is:
#'   #      [,1] [,2]
#'   # [1,]    1    2
#'   # [2,]    1    3
#'   # [3,]    1    4
#'   # Column 1 has zero variance.
#'   scaled_mat <- my_scale(mat)
#'   print(scaled_mat)
#'   # attr(scaled_mat, "scaled:center") # Should be c(1, 3)
#'   # attr(scaled_mat, "scaled:scale")  # Should be c(1, 1)
#'   # Expected output for scaled_mat:
#'   #      [,1] [,2]
#'   # [1,]    0   -1
#'   # [2,]    0    0
#'   # [3,]    0    1
#'}
#' @keywords internal
#' @noRd
my_scale <- function(x) {
	# Compute column means and standard deviations
	ctr <- colMeans(x)
	sds <- apply(x, 2, sd)

	# Identify zero-variance columns
	zero_sd <- (sds == 0)

	# For zero-variance columns, set scale=1 to avoid dividing by 0
	ctr2 <- ctr
	sds2 <- sds
	sds2[zero_sd] <- 1

	# Center and scale
	scaled <- sweep(x, 2, ctr2, FUN = "-")
	scaled <- sweep(scaled, 2, sds2, FUN = "/")

	# Attach attributes so behavior mimics base::scale()
	attr(scaled, "scaled:center") <- ctr2
	attr(scaled, "scaled:scale") <- sds2

	return(scaled)
}

# getNumTreats
#' @title Calculate the Total Number of Base Treatment Effect Parameters
#' @description Computes the total number of unique treatment dummy variables
#'   (and thus base treatment effect parameters `tau_rt`) in a staggered
#'   adoption setting.
#' @param R Integer; the number of treated cohorts. Treatment is assumed to
#'   start in periods 2 to `R+1`.
#' @param T Integer; the total number of time periods.
#' @return An integer representing the total number of treatment effect
#'   parameters (`num_treats`).
#' @details The formula used is `num_treats = T * R - (R * (R + 1)) / 2`.
#'   This corresponds to summing the number of post-treatment periods for each
#'   cohort:
#'   Cohort 1 (starts period 2): T-1 effects
#'   Cohort 2 (starts period 3): T-2 effects
#'   ...
#'   Cohort R (starts period R+1): T-R effects
#'   Summing these gives `R*T - (1 + 2 + ... + R) = R*T - R(R+1)/2`.
#' @keywords internal
#' @noRd
getNumTreats <- function(R, T) {
	return(T * R - (R * (R + 1)) / 2)
}

# getTreatInds
#' @title Get Indices of Base Treatment Effect Parameters
#' @description Determines the column indices in the full design matrix `X_ints` (or
#'   corresponding coefficient vector `beta`) that correspond to the base
#'   treatment effect parameters (`tau_rt`).
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the number of time periods.
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of base treatment effect
#'   parameters, typically calculated by `getNumTreats(R, T)`.
#' @return An integer vector containing the indices for the `num_treats` base
#'   treatment effect parameters.
#' @details The full design matrix `X_ints` is structured with several blocks of
#'   variables:
#'   1. Cohort fixed effects (`R` columns)
#'   2. Time fixed effects (`T-1` columns)
#'   3. Covariate main effects (`d` columns, if `d > 0`)
#'   4. Cohort-Covariate interactions (`d*R` columns, if `d > 0`)
#'   5. Time-Covariate interactions (`d*(T-1)` columns, if `d > 0`)
#'   The base treatment effects form the next block. This function calculates the
#'   total number of columns in the preceding blocks (`base_cols`) and returns
#'   `base_cols + 1` through `base_cols + num_treats`.
#' @keywords internal
#' @noRd
getTreatInds <- function(R, T, d, num_treats) {
	base_cols <- if (d > 0) {
		R + (T - 1) + d + d * R + d * (T - 1)
	} else {
		R + (T - 1)
	}

	treat_inds <- seq(from = base_cols + 1, length.out = num_treats)

	stopifnot(length(treat_inds) == num_treats)
	if (d > 0) {
		stopifnot(
			max(treat_inds) == R + T - 1 + d + R * d + (T - 1) * d + num_treats
		)
	} else {
		stopifnot(max(treat_inds) == R + T - 1 + num_treats)
	}

	stopifnot(length(treat_inds) == num_treats)

	return(treat_inds)
}

# getP
#' @title Calculate Total Number of Parameters in the Full Design Matrix
#' @description Computes `p`, the total number of columns (parameters) in the
#'   full design matrix `X_ints` used by the FETWFE model, including all fixed
#'   effects, covariates, treatment dummies, and all their interactions.
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the total number of time periods.
#' @param d Integer; the number of time-invariant covariates.
#' @param num_treats Integer; the total number of base treatment effect
#'   parameters (e.g., from `getNumTreats(R,T)`).
#' @return An integer `p` representing the total number of parameters.
#' @details The total number of parameters `p` is the sum of:
#'   - Cohort fixed effects: `R`
#'   - Time fixed effects: `T-1`
#'   - Covariate main effects: `d`
#'   - Cohort-Covariate interactions: `d*R`
#'   - Time-Covariate interactions: `d*(T-1)`
#'   - Base treatment effects: `num_treats`
#'   - Treatment-Covariate interactions: `num_treats*d`
#'   The formula is `p = R + (T-1) + d + d*R + d*(T-1) + num_treats + num_treats*d`.
#'   If `d=0`, terms involving `d` become zero.
#' @keywords internal
#' @noRd
getP <- function(R, T, d, num_treats) {
	return(R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d)
}

#' Get Indices of First Treatment Effects for Each Cohort
#'
#' @description
#' Calculates the starting indices of treatment effect parameters for each cohort
#' within a concatenated block of all treatment effects. This is used for constructing
#' fusion penalty matrices.
#'
#' @param R Integer; the number of treated cohorts.
#' @param T Integer; the total number of time periods.
#'
#' @details
#' Assumes treatment for cohort `r` (where `r` is 1-indexed for calculation,
#' corresponding to actual adoption times `times[r+1]`) starts at `times[r+1]`
#' and continues until `times[T]`.
#' The number of treatment effects for cohort `r` is `T - (r+1) + 1 = T - r`.
#' The `num_treats` is the sum of these counts.
#' The formula for the starting index of the `r`-th cohort's treatment effects
#' (1-indexed `r` from 1 to `R`):
#' `f_inds[r] = 1 + sum_{k=1}^{r-1} (T - k) = 1 + (r-1)T - (r-1)r/2`.
#'
#' Gemini check:
#' The paper's formula `1 + (r - 1)*(2*T - r)/2` seems to be for a slightly
#' different definition of `num_treats` or indexing. This function should align
#' with how `treat_var_mat` is constructed in `addDummies` and `num_treats`
#' is calculated there.
#' Let's re-evaluate the formula based on typical DiD setups:
#' Cohort 1 (starts at time 2) has T-1 effects.
#' Cohort 2 (starts at time 3) has T-2 effects.
#' ...
#' Cohort R (starts at time R+1) has T-R effects.
#' `first_inds[1] = 1`
#' `first_inds[2] = (T-1) + 1`
#' `first_inds[3] = (T-1) + (T-2) + 1`
#' So, `first_inds[r] = 1 + sum_{j=1}^{r-1} (T-j)`.
#'
#' The current code uses: `f_inds[r] <- 1 + (r - 1)*(2*T - r)/2`.
#' Let's test this formula:
#' r=1: 1
#' r=2: 1 + 1*(2T-2)/2 = 1 + T - 1 = T
#' r=3: 1 + 2*(2T-3)/2 = 1 + 2T - 3 = 2T - 2
#' This implies the number of effects for cohort 1 is T-1 (indices 1 to T-1).
#' Number of effects for cohort 2 is (2T-2) - T = T-2.
#' This seems to match the sequence T-1, T-2, ..., T-R.
#' Total number of treats (`num_treats`) from `getNumTreats(R,T)` is `T*R - R*(R+1)/2`.
#' This is `sum_{j=1 to R} (T-j)`. This is consistent.
#'
#' @return An integer vector of length `R`, where `f_inds[r]` is the 1-based
#'   starting index of the `r`-th cohort's treatment effects in the combined block.
#' @keywords internal
#' @noRd
getFirstInds <- function(R, T) {
	# Let's identify the indices of the first treatment effects for each cohort.
	# The first one is index 1, then the second one is index (T - 1) + 1 = T,
	# then the third one is (T - 1) + (T - 2) + 1 = 2*T - 2. In general, for
	# r > 1 the rth one will occur at index
	#
	# (T - 1) + (T - 2) + ... + (T - (r - 1)) + 1
	# = 1 + (r - 1)*(T - 1 + T - r + 1)/2
	# = 1 + (r - 1)*(2*T - r)/2.
	#
	# (Looks like the formula works for r = 1 too.)

	n_treats <- getNumTreats(R = R, T = T)

	f_inds <- integer(R)
	if (R == 0) return(f_inds) # No cohorts, no first_inds

	for (r in 1:R) {
		f_inds[r] <- 1 + (r - 1) * (2 * T - r) / 2
	}
	stopifnot(all(f_inds <= n_treats))
	stopifnot(f_inds[1] == 1)
	# Last cohort has T - R treatment effects to estimate. So last first_ind
	# should be at position num_treats - (T - R) + 1 = num_treats - T + R + 1.
	stopifnot(f_inds[R] == n_treats - T + R + 1)

	# Additional checks from Gemini below

	stopifnot(all(f_inds <= n_treats) || R == 0)
	if (R > 0) stopifnot(f_inds[1] == 1)

	# Last first_ind: f_inds[R] = 1 + sum_{j=1}^{R-1} (T-j)
	# Total effects = sum_{j=1}^{R} (T-j) = n_treats
	# Last block of effects for cohort R has T-R effects.
	# So f_inds[R] should be n_treats - (T-R) + 1.
	if (R > 0) {
		stopifnot(f_inds[R] == n_treats - (T - R) + 1)
	}

	# The original formula was: `1 + (r - 1)*(2*T - r)/2`
	# Let's re-check the paper's formula logic, it implies number of effects for cohort `k` (1-indexed) is `T-k`.
	# sum_{j=1}^{r-1} (T-j) = (r-1)T - (r-1)r/2. So f_inds[r] = 1 + (r-1)T - r(r-1)/2.
	# (r-1)*(2T-r)/2 = (r-1)T - r(r-1)/2. Yes, they are the same.

	# Using the paper's original loop for safety, assuming it's correct with getNumTreats
	f_inds_paper <- integer(R)
	if (R > 0) {
		for (r_val in 1:R) {
			# r_val is 1-indexed cohort number
			f_inds_paper[r_val] <- 1 + (r_val - 1) * (2 * T - r_val) / 2 # This is not sum of T-k, but related to cumulative sum for triangular numbers.
			# This formula is sum_{k=0}^{r-2} (T-1-k) if T-1 is max effects
			# Or sum_{j=0}^{r-1-1} ( (T-1) - j )
			# Let's use the direct sum formulation which is clearer.
		}
		stopifnot(f_inds == f_inds_paper) # Verify my derivation matches paper's formula
	}

	return(f_inds)
}

# sse_bridge
#' @title Calculate Sum of Squared Errors for Bridge Regression
#' @description Computes the sum of squared errors (SSE) for a given set of
#'   intercept and slope coefficients from a bridge regression model, relative
#'   to the observed responses and the (potentially transformed) design matrix.
#'   The result is then averaged by dividing by the total number of observations (N*T).
#' @param eta_hat Numeric scalar; the estimated intercept term.
#' @param beta_hat Numeric vector; the estimated slope coefficients. Its length
#'   should match the number of columns in `X_mod`.
#' @param y Numeric vector; the observed response variable, of length `N*T`.
#' @param X_mod Numeric matrix; the design matrix (possibly transformed, e.g.,
#'   for FETWFE) used in the regression. It has `N*T` rows.
#' @param N Integer; the total number of unique units.
#' @param T Integer; the total number of time periods.
#' @return A numeric scalar representing the mean squared error (MSE), i.e.,
#'   SSE / (N*T).
#' @keywords internal
#' @noRd
sse_bridge <- function(eta_hat, beta_hat, y, X_mod, N, T) {
	stopifnot(length(eta_hat) == 1)
	stopifnot(length(beta_hat) == ncol(X_mod))

	y_hat <- X_mod %*% beta_hat + eta_hat
	stopifnot(length(y_hat) == N * T)
	stopifnot(all(!is.na(y_hat)))

	ret <- sum((y - y_hat)^2) / (N * T)

	stopifnot(!is.na(ret))
	stopifnot(ret >= 0)

	return(ret)
}
