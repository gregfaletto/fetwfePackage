#' @import glmnet
#' @importFrom Matrix bdiag
#' @importFrom stats qnorm pnorm predict coef model.matrix setNames lm

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
idCohorts <- function(df, time_var, unit_var, treat_var) {
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

	# Collect per-unit violations across BOTH balance and absorbing-state
	# checks, then stop() once at the end with a grouped listing. Previously
	# each check stop()-ed on the first offending unit, forcing the user
	# into one error round-trip per malformed unit. See issue #64.
	balance_violations <- character(0)
	absorbing_violations <- character(0)

	for (s in units) {
		df_s <- df[df[, unit_var] == s, ]

		# Balance check (gates the absorbing-state check on the same unit;
		# we can't trust the per-period treatment vector when observations
		# are missing or when a time period is duplicated). The
		# setequal() clause catches the (1, 2, 2)-style duplicate-time
		# case that nrow == T silently passed pre-#75.
		if (nrow(df_s) != T || !setequal(df_s[, time_var], times)) {
			balance_violations <- c(
				balance_violations,
				paste0(
					"unit ",
					s,
					" has ",
					nrow(df_s),
					" observations (",
					length(unique(df_s[, time_var])),
					" distinct time periods)"
				)
			)
			next
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

			# Absorbing-state check.
			# Check from the actual_treat_time onwards in the original df for unit s
			original_unit_times_from_treatment <- df[
				(df[, unit_var] == s) & (df[, time_var] >= actual_treat_time),
				treat_var
			]
			if (any(original_unit_times_from_treatment != 1)) {
				absorbing_violations <- c(absorbing_violations, as.character(s))
				# Do NOT append to cohorts: only units that pass BOTH the
				# balance and absorbing-state checks contribute to cohort
				# membership. (load-bearing: without `next`, this unit would
				# leak into the cohorts list before stop() fires below.)
				next
			}
			cohorts[[as.character(actual_treat_time)]] <- c(
				cohorts[[as.character(actual_treat_time)]],
				s
			)
		}
	}

	if (
		length(balance_violations) > 0 ||
			length(absorbing_violations) > 0
	) {
		# unit_var is required to be character (enforced at the public
		# entry points; see `.collect_etwfe_input_violations()` later in
		# this file), so sort() is lexicographic. For unit names like
		# "unit10" / "unit2", the ordering is "unit10" < "unit2" <
		# "unit3"; that's the same lex order as v1.9.3's cohort-sort
		# caveat (#53).
		msgs <- character(0)
		if (length(balance_violations) > 0) {
			msgs <- c(
				msgs,
				paste0(
					"Panel does not appear to be balanced (expected T = ",
					T,
					" observations per unit). Violations:\n  ",
					paste(
						.truncate_violations(sort(balance_violations)),
						collapse = "\n  "
					)
				)
			)
		}
		if (length(absorbing_violations) > 0) {
			msgs <- c(
				msgs,
				paste0(
					"Treatment does not appear to be an absorbing state. ",
					"Violations (units whose treatment was 1 and later 0):\n  ",
					paste(
						.truncate_violations(sort(absorbing_violations)),
						collapse = "\n  "
					)
				)
			)
		}
		stop(paste(msgs, collapse = "\n"))
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

#' @title Truncate a violations listing for an error message
#' @description Used by `idCohorts()` to cap a long list of per-unit
#'   violations at a small representative sample, so a user with hundreds
#'   of malformed units does not get a multi-page error message. Returns
#'   the first `max_show` entries verbatim followed by `"... and N more"`
#'   when the input exceeds `max_show`; returns the input unchanged
#'   otherwise.
#' @param xs Character vector of violation descriptions.
#' @param max_show Integer scalar; maximum number of entries to show
#'   verbatim before truncating. Defaults to 20.
#' @return Character vector of length `min(length(xs), max_show + 1)`.
#' @keywords internal
#' @noRd
.truncate_violations <- function(xs, max_show = 20L) {
	if (length(xs) <= max_show) {
		return(xs)
	}
	c(
		xs[seq_len(max_show)],
		paste0("... and ", length(xs) - max_show, " more")
	)
}


# .collect_etwfe_input_violations
#' @title Collect input-validation violations for the etwfe-family entry points
#' @description Internal helper that runs the per-argument checks of
#'   `checkEtwfeInputs()` without `stop()`-ing on the first failure. Each
#'   malformed argument contributes one human-readable line to a
#'   character vector that the caller (`checkEtwfeInputs()` or
#'   `checkFetwfeInputs()`) consolidates into a single multi-line
#'   `stop()` message. This collect-all-violations pattern means a user
#'   with K malformed args sees K violations on the first call instead
#'   of one error round-trip per arg.
#' @return A list with three elements:
#'   - `violations`: character vector of human-readable per-arg violation
#'     lines (empty when all inputs are valid).
#'   - `pdata`: the (possibly tibble-coerced) `pdata`. Coercion runs
#'     regardless of downstream violations so the caller does not have
#'     to redo it.
#'   - `indep_count_data_available`: logical scalar; `TRUE` if a valid
#'     `indep_counts` argument was provided, `FALSE` otherwise.
#' @details Within a single argument's check chain (e.g., `time_var`
#'   must be character, length 1, present in `pdata`, integer column),
#'   the first failure short-circuits — downstream checks for that arg
#'   are unsafe to run (e.g., `time_var %in% colnames(pdata)` cannot
#'   sensibly be evaluated if `time_var` is not a length-1 character).
#'   Across DIFFERENT arguments, violations collect independently — if
#'   both `time_var` and `unit_var` are malformed, BOTH violations
#'   appear in the final list. This is "cascade within arg, collect
#'   across args".
#' @keywords internal
#' @noRd
.collect_etwfe_input_violations <- function(
	pdata,
	time_var,
	unit_var,
	treatment,
	response,
	covs = c(),
	indep_counts = NA,
	sig_eps_sq = NA,
	sig_eps_c_sq = NA,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE
) {
	violations <- character(0)

	# `pdata` itself must be a data frame before we can probe its
	# columns. If it is not, every downstream column-existence /
	# column-type check is meaningless; record the violation and skip
	# all pdata-dependent checks.
	pdata_ok <- is.data.frame(pdata)
	if (!pdata_ok) {
		violations <- c(
			violations,
			sprintf(
				"pdata must be a data frame; got %s",
				paste(class(pdata), collapse = "/")
			)
		)
	} else {
		# Tibble -> data frame coercion; harmless when pdata is already
		# a plain data.frame.
		if ("tbl_df" %in% class(pdata)) {
			pdata <- as.data.frame(pdata)
		}
		if (nrow(pdata) < 4) {
			violations <- c(
				violations,
				sprintf(
					"pdata must have at least 4 rows (2 units at 2 times); got %d",
					nrow(pdata)
				)
			)
		}
	}

	# --- time_var --------------------------------------------------------
	if (!is.character(time_var) || length(time_var) != 1) {
		violations <- c(
			violations,
			sprintf(
				"time_var must be a single character string; got %s (length %d)",
				paste(class(time_var), collapse = "/"),
				length(time_var)
			)
		)
	} else if (pdata_ok && !(time_var %in% colnames(pdata))) {
		violations <- c(
			violations,
			sprintf(
				"time_var = '%s' is not a column in pdata; columns are: %s",
				time_var,
				paste(colnames(pdata), collapse = ", ")
			)
		)
	} else if (pdata_ok && !is.integer(pdata[[time_var]])) {
		violations <- c(
			violations,
			sprintf(
				"the time_var column '%s' must be integer; got %s",
				time_var,
				paste(class(pdata[[time_var]]), collapse = "/")
			)
		)
	}

	# --- unit_var --------------------------------------------------------
	if (!is.character(unit_var) || length(unit_var) != 1) {
		violations <- c(
			violations,
			sprintf(
				"unit_var must be a single character string; got %s (length %d)",
				paste(class(unit_var), collapse = "/"),
				length(unit_var)
			)
		)
	} else if (pdata_ok && !(unit_var %in% colnames(pdata))) {
		violations <- c(
			violations,
			sprintf(
				"unit_var = '%s' is not a column in pdata; columns are: %s",
				unit_var,
				paste(colnames(pdata), collapse = ", ")
			)
		)
	} else if (pdata_ok && !is.character(pdata[[unit_var]])) {
		violations <- c(
			violations,
			sprintf(
				"the unit_var column '%s' must be character; got %s",
				unit_var,
				paste(class(pdata[[unit_var]]), collapse = "/")
			)
		)
	}

	# --- treatment -------------------------------------------------------
	if (!is.character(treatment) || length(treatment) != 1) {
		violations <- c(
			violations,
			sprintf(
				"treatment must be a single character string; got %s (length %d)",
				paste(class(treatment), collapse = "/"),
				length(treatment)
			)
		)
	} else if (pdata_ok && !(treatment %in% colnames(pdata))) {
		violations <- c(
			violations,
			sprintf(
				"treatment = '%s' is not a column in pdata; columns are: %s",
				treatment,
				paste(colnames(pdata), collapse = ", ")
			)
		)
	} else if (pdata_ok && !is.integer(pdata[[treatment]])) {
		violations <- c(
			violations,
			sprintf(
				"the treatment column '%s' must be integer; got %s",
				treatment,
				paste(class(pdata[[treatment]]), collapse = "/")
			)
		)
	} else if (pdata_ok && !all(pdata[, treatment] %in% c(0, 1))) {
		bad_vals <- unique(pdata[, treatment][
			!(pdata[, treatment] %in% c(0, 1))
		])
		violations <- c(
			violations,
			sprintf(
				"the treatment column '%s' must contain only 0/1; saw other value(s): %s",
				treatment,
				paste(utils::head(bad_vals, 5), collapse = ", ")
			)
		)
	}

	# --- covs ------------------------------------------------------------
	if (length(covs) > 0) {
		if (!is.character(covs)) {
			violations <- c(
				violations,
				sprintf(
					"covs must be a character vector of column names; got %s",
					paste(class(covs), collapse = "/")
				)
			)
		} else if (pdata_ok && !all(covs %in% colnames(pdata))) {
			missing_covs <- covs[!(covs %in% colnames(pdata))]
			violations <- c(
				violations,
				sprintf(
					"covs contains name(s) not in pdata: %s",
					paste(missing_covs, collapse = ", ")
				)
			)
		} else if (pdata_ok) {
			bad_cov_types <- character(0)
			for (cov in covs) {
				if (
					!(is.numeric(pdata[[cov]]) ||
						is.integer(pdata[[cov]]) ||
						is.factor(pdata[[cov]]))
				) {
					bad_cov_types <- c(
						bad_cov_types,
						sprintf(
							"'%s' (%s)",
							cov,
							paste(class(pdata[[cov]]), collapse = "/")
						)
					)
				}
			}
			if (length(bad_cov_types) > 0) {
				violations <- c(
					violations,
					sprintf(
						"each cov column must be numeric, integer, or factor; offending column(s): %s",
						paste(bad_cov_types, collapse = ", ")
					)
				)
			}
		}
	}

	# --- response --------------------------------------------------------
	if (!is.character(response) || length(response) != 1) {
		violations <- c(
			violations,
			sprintf(
				"response must be a single character string; got %s (length %d)",
				paste(class(response), collapse = "/"),
				length(response)
			)
		)
	} else if (pdata_ok && !(response %in% colnames(pdata))) {
		violations <- c(
			violations,
			sprintf(
				"response = '%s' is not a column in pdata; columns are: %s",
				response,
				paste(colnames(pdata), collapse = ", ")
			)
		)
	} else if (
		pdata_ok &&
			!(is.numeric(pdata[[response]]) || is.integer(pdata[[response]]))
	) {
		violations <- c(
			violations,
			sprintf(
				"the response column '%s' must be numeric or integer; got %s",
				response,
				paste(class(pdata[[response]]), collapse = "/")
			)
		)
	}

	# --- indep_counts ----------------------------------------------------
	# Mirrors original: only validated when not all-NA. If valid, the
	# `indep_count_data_available` flag is TRUE.
	indep_count_data_available <- FALSE
	if (any(!is.na(indep_counts))) {
		if (!is.integer(indep_counts)) {
			violations <- c(
				violations,
				sprintf(
					"indep_counts must be an integer vector; got %s",
					paste(class(indep_counts), collapse = "/")
				)
			)
		} else if (any(indep_counts <= 0)) {
			# Original message preserved verbatim (also surfaced in
			# `etwfeWithSimulatedData()` regression coverage).
			violations <- c(
				violations,
				"At least one cohort in the independent count data has 0 members"
			)
		} else {
			indep_count_data_available <- TRUE
		}
	}

	# --- sig_eps_sq ------------------------------------------------------
	if (any(!is.na(sig_eps_sq))) {
		if (!(is.numeric(sig_eps_sq) || is.integer(sig_eps_sq))) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_sq must be numeric or integer; got %s",
					paste(class(sig_eps_sq), collapse = "/")
				)
			)
		} else if (length(sig_eps_sq) != 1) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_sq must be a single value; got length %d",
					length(sig_eps_sq)
				)
			)
		} else if (sig_eps_sq < 0) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_sq must be non-negative; got %s",
					format(sig_eps_sq)
				)
			)
		}
	}

	# --- sig_eps_c_sq ----------------------------------------------------
	if (any(!is.na(sig_eps_c_sq))) {
		if (!(is.numeric(sig_eps_c_sq) || is.integer(sig_eps_c_sq))) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_c_sq must be numeric or integer; got %s",
					paste(class(sig_eps_c_sq), collapse = "/")
				)
			)
		} else if (length(sig_eps_c_sq) != 1) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_c_sq must be a single value; got length %d",
					length(sig_eps_c_sq)
				)
			)
		} else if (sig_eps_c_sq < 0) {
			violations <- c(
				violations,
				sprintf(
					"sig_eps_c_sq must be non-negative; got %s",
					format(sig_eps_c_sq)
				)
			)
		}
	}

	# --- verbose ---------------------------------------------------------
	if (!is.logical(verbose) || length(verbose) != 1) {
		violations <- c(
			violations,
			sprintf(
				"verbose must be a single logical (TRUE or FALSE); got %s (length %d)",
				paste(class(verbose), collapse = "/"),
				length(verbose)
			)
		)
	}

	# --- alpha -----------------------------------------------------------
	# Match original: only emit the alpha > 0.5 warning when alpha is
	# otherwise fully valid (numeric, length 1, in (0, 1)). The warning
	# is a UX nudge and was unreachable when any alpha check failed.
	alpha_valid <- FALSE
	if (!is.numeric(alpha) || length(alpha) != 1) {
		violations <- c(
			violations,
			sprintf(
				"alpha must be a single numeric; got %s (length %d)",
				paste(class(alpha), collapse = "/"),
				length(alpha)
			)
		)
	} else if (!(alpha > 0 && alpha < 1)) {
		violations <- c(
			violations,
			sprintf(
				"alpha must be in (0, 1); got %s",
				format(alpha)
			)
		)
	} else {
		alpha_valid <- TRUE
	}
	if (alpha_valid && alpha > 0.5) {
		warning(
			"Provided alpha > 0.5; are you sure you didn't mean to enter a smaller alpha? The confidence level will be 1 - alpha."
		)
	}

	# --- add_ridge -------------------------------------------------------
	if (!is.logical(add_ridge) || length(add_ridge) != 1) {
		violations <- c(
			violations,
			sprintf(
				"add_ridge must be a single logical (TRUE or FALSE); got %s (length %d)",
				paste(class(add_ridge), collapse = "/"),
				length(add_ridge)
			)
		)
	}

	list(
		violations = violations,
		pdata = pdata,
		indep_count_data_available = indep_count_data_available
	)
}

# .format_input_violations
#' @title Assemble a multi-line input-violations error message
#' @description Internal helper used by `checkEtwfeInputs()` and
#'   `checkFetwfeInputs()` to render a vector of per-arg violation
#'   lines as a multi-line `stop()` body. Centralising the header /
#'   prefix here means both validators speak with one voice, and the
#'   four entry points (`fetwfe`, `etwfe`, `betwfe`, `twfeCovs`)
#'   continue to produce byte-identical messages — a property the
#'   PR #103 snapshot guardrail asserts.
#' @keywords internal
#' @noRd
.format_input_violations <- function(violations) {
	paste0(
		"Invalid inputs:\n  - ",
		paste(violations, collapse = "\n  - ")
	)
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
#' Equivalently, `1 + (r - 1)*(2*T - r)/2` (algebraic restatement).
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
	if (R == 0) {
		return(f_inds)
	} # No cohorts, no first_inds

	for (r in 1:R) {
		f_inds[r] <- 1 + (r - 1) * (2 * T - r) / 2
	}
	stopifnot(all(f_inds <= n_treats))
	stopifnot(f_inds[1] == 1)
	# Last cohort has T - R treatment effects to estimate. So last first_ind
	# should be at position num_treats - (T - R) + 1 = num_treats - T + R + 1.
	stopifnot(f_inds[R] == n_treats - T + R + 1)

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

#' @title Compute two-sided p-values from estimates and standard errors
#'
#' @description
#' Internal helper used by `summary.*` / `print.*` methods and `catt_df`
#' assemblers across the four estimator classes. Returns the standard
#' `2 * pnorm(-|t|)` p-value for rows where the standard error is positive,
#' and `NA_real_` for rows where the standard error is zero or `NA`. The
#' zero-SE case happens when the estimator zeros a coefficient (selected
#' out under FETWFE / BETWFE bridge regression with `q < 1`), or in
#' degenerate fallback paths where SEs cannot be computed. Routing all
#' p-value computation through this single helper avoids the
#' `0 / 0 = NaN` trap and gives one canonical formula to audit.
#'
#' @param estimates Numeric vector of coefficient estimates.
#' @param ses Numeric vector of standard errors, same length as `estimates`.
#' @return A numeric vector the same length as `estimates`. Each entry is
#'   `2 * pnorm(-abs(estimates[i] / ses[i]))` when `ses[i] > 0` and not
#'   `NA`, and `NA_real_` otherwise.
#' @keywords internal
#' @noRd
.compute_p_values <- function(estimates, ses) {
	ifelse(
		!is.na(ses) & ses > 0,
		2 * pnorm(-abs(estimates / ses)),
		NA_real_
	)
}

#' @title Auto-truncate a panel that has no never-treated units
#'
#' @description
#' Internal helper called by the four public estimator entry points
#' (`fetwfe()`, `etwfe()`, `betwfe()`, `twfeCovs()`) when their
#' `allow_no_never_treated` argument is TRUE (default). Detects whether
#' the input panel `pdata` contains any never-treated units; if so,
#' returns `pdata` unchanged. If not, identifies the latest cohort start
#' time and drops rows at or after that time, issues a warning naming the
#' dropped time periods, and returns the truncated panel. If
#' `allow_no_never_treated = FALSE` and there are no never-treated units,
#' stops with an informative error. If truncation is structurally
#' impossible (would leave fewer than 2 time periods or fewer than 2
#' treated cohorts), stops regardless of the argument.
#'
#' @param pdata A balanced panel `data.frame`.
#' @param time_var Character string; name of the time variable column.
#' @param unit_var Character string; name of the unit identifier column.
#' @param treat_var Character string; name of the treatment indicator column.
#' @param allow_no_never_treated Logical; if `FALSE`, stops with an error
#'   when there are no never-treated units. If `TRUE` (default), auto-
#'   truncates.
#' @details The helper trusts the package's documented absorbing-state
#'   assumption for the treatment column. If a unit's treatment is
#'   non-absorbing (e.g., 0 -> 1 -> 0 -> 1), and the rows that violate
#'   absorbing would be dropped by truncation, this helper will not catch
#'   the violation (the truncated panel sees a clean absorbing history).
#'   Users supplying non-absorbing data are already violating the
#'   package's input contract; the responsibility for catching that sits
#'   at `idCohorts` in the standard pipeline.
#' @return The (possibly truncated) `pdata` data.frame.
#' @keywords internal
#' @noRd
.truncate_if_no_never_treated <- function(
	pdata,
	time_var,
	unit_var,
	treat_var,
	allow_no_never_treated = TRUE
) {
	stopifnot(is.logical(allow_no_never_treated))
	stopifnot(length(allow_no_never_treated) == 1)
	stopifnot(!is.na(allow_no_never_treated))

	# Detect never-treated units.
	treated_any <- tapply(
		pdata[[treat_var]],
		pdata[[unit_var]],
		function(x) any(x == 1)
	)
	if (any(!treated_any)) {
		return(pdata)
	}

	# All units are eventually treated.
	if (!isTRUE(allow_no_never_treated)) {
		stop(
			"No never-treated units detected in data to fit model; ",
			"estimating treatment effects is not possible. Set ",
			"`allow_no_never_treated = TRUE` to enable automatic panel ",
			"truncation, which drops time periods at and after the latest ",
			"cohort's start time so the latest-cohort units become the ",
			"never-treated comparison group."
		)
	}

	# Identify each unit's first-treatment time.
	unit_first_treat <- tapply(
		seq_len(nrow(pdata)),
		pdata[[unit_var]],
		function(idx) {
			treated_rows <- idx[pdata[[treat_var]][idx] == 1]
			min(pdata[[time_var]][treated_rows])
		}
	)
	r_max <- max(unit_first_treat)
	times_sorted <- sort(unique(pdata[[time_var]]))

	# If `r_max == times[1]`, every unit was treated from the very first
	# time period (otherwise `r_max` would be larger). There is no
	# pre-treatment period to retain; the existing pipeline's first-period
	# filter in `idCohorts()` already reports this case with a clearer
	# message ("All units were treated in the first time period..."), so
	# defer to it rather than rewrite the diagnostic.
	if (r_max <= times_sorted[1]) {
		return(pdata)
	}

	retained_times <- times_sorted[times_sorted < r_max]
	dropped_times <- times_sorted[times_sorted >= r_max]

	# Impossible-truncation guards.
	if (length(retained_times) < 2) {
		stop(
			"Cannot estimate treatment effects: panel has no never-treated ",
			"units, and truncating to the largest sub-panel where any ",
			"units are still untreated would leave only ",
			length(retained_times),
			" time period(s) (need at least 2)."
		)
	}

	# After truncation, retained cohorts are those with r_i < r_max.
	# FETWFE / ETWFE / BETWFE / twfeCovs each require R >= 2 (enforced
	# inside `prep_for_etwfe_core()` in R/core_funcs.R with a user-facing
	# error); the guard here matches that constraint so the user gets a
	# clear truncation-specific message rather than the deeper R >= 2
	# error after auto-truncation.
	retained_cohorts <- setdiff(unique(unit_first_treat), r_max)
	if (length(retained_cohorts) < 2) {
		stop(
			"Cannot estimate treatment effects: panel has no never-treated ",
			"units, and the truncated panel would have only ",
			length(retained_cohorts),
			" treated cohort(s). FETWFE/ETWFE/BETWFE/twfeCovs each require ",
			"at least 2 treated cohorts."
		)
	}

	# Truncate.
	truncated <- pdata[pdata[[time_var]] %in% retained_times, , drop = FALSE]
	warning(
		"No never-treated units in input data; auto-truncated panel by ",
		"dropping time periods at or after ",
		r_max,
		" (dropped: ",
		paste(dropped_times, collapse = ", "),
		"). The units that started treatment at ",
		r_max,
		" serve as the never-treated comparison group in the truncated ",
		"panel. Set `allow_no_never_treated = FALSE` to disable this ",
		"behavior."
	)
	truncated
}

#' @title Cohort-r treatment-coefficient index range
#' @description
#' For cohort `r` (one of `1:R`), returns the integer range
#' `first_ind_r:last_ind_r` indexing into a length-`num_treats`
#' treatment-coefficient vector. The closed-form is
#' `first_ind_r = first_inds[r]` and
#' `last_ind_r = if (r < R) first_inds[r + 1] - 1 else num_treats`.
#' Used inside `for (r in 1:R)` cohort loops in `getCohortATTsFinal`,
#' `getSecondVarTermDataApp`, `getSecondVarTermOLS`,
#' `prep_for_etwfe_regression`, and `getActualCohortTes`. Consolidated by
#' GitHub #83.
#' @param r Integer; the cohort index, `1 <= r <= R`.
#' @param R Integer; total number of treated cohorts.
#' @param first_inds Integer vector; the per-cohort starting indices
#'   (length `R`), as returned by `getFirstInds()`.
#' @param num_treats Integer; the total number of treatment coefficients
#'   (one per (cohort, time) pair), as returned by `getNumTreats()`.
#' @return Integer vector; the range `first_ind_r:last_ind_r`.
#' @seealso `getFirstInds()`, `getNumTreats()`.
#' @keywords internal
#' @noRd
.cohort_block_inds <- function(r, R, first_inds, num_treats) {
	first_ind_r <- first_inds[r]
	last_ind_r <- if (r < R) first_inds[r + 1] - 1 else num_treats
	first_ind_r:last_ind_r
}

#' @title Multinomial covariance matrix from cohort-membership probabilities
#' @description
#' Returns the `length(probs) x length(probs)` covariance matrix of the
#' multinomial random variable with cell probabilities `probs` and a
#' single trial (or, equivalently, the per-observation covariance of the
#' indicator vector for cohort membership). The closed-form is
#' `S = -outer(probs, probs); diag(S) <- probs * (1 - probs)`.
#' Used inside the variance-term-2 machinery at `getSecondVarTermDataApp`
#' (R/variance_machinery.R), `getSecondVarTermOLS` (R/variance_machinery.R),
#' and `.event_study_var2_fetwfe` (R/event_study.R). Consolidated by
#' GitHub #83.
#' @param probs Numeric vector of length `R`; cohort-membership
#'   probabilities P(W = r) for r in 1:R. Should sum to less than 1
#'   (the residual mass is the never-treated probability).
#' @return Numeric matrix of dimensions `length(probs) x length(probs)`;
#'   the symmetric covariance matrix.
#' @keywords internal
#' @noRd
.multinomial_cov <- function(probs) {
	S <- -outer(probs, probs)
	diag(S) <- probs * (1 - probs)
	S
}

#' @title Run the shared input-prep pipeline for the four estimator entry points
#' @description
#' Internal helper that consolidates the verbatim Steps 3-5 (input
#' validation + auto-truncation + design-matrix prep) shared across
#' `fetwfe()`, `etwfe()`, `betwfe()`, and `twfeCovs()`. Branches on
#' `estimator_type` to pick the correct validator
#' (`checkFetwfeInputs()` for fetwfe/betwfe, `checkEtwfeInputs()` for
#' etwfe/twfeCovs), then calls `.truncate_if_no_never_treated()` and
#' `prep_for_etwfe_core()`. The four caller-side rewrites consume the
#' returned list with verbose per-field unpacking. Consolidated by
#' GitHub #79.
#' @param pdata,time_var,unit_var,treatment,response,covs,indep_counts,sig_eps_sq,sig_eps_c_sq,verbose,alpha,add_ridge
#'   The shared arguments validated by both `checkEtwfeInputs()` and
#'   `checkFetwfeInputs()`. See either validator's source for the
#'   per-argument contracts.
#' @param lambda.max,lambda.min,q The fetwfe/betwfe-specific bridge
#'   penalty arguments. Forwarded to `checkFetwfeInputs()` when
#'   `estimator_type == "fetwfe"`. When `estimator_type == "etwfe"`,
#'   these must remain `NA` (their defaults); a `stopifnot()` guards
#'   caller-side misuse.
#' @param allow_no_never_treated Logical; forwarded to
#'   `.truncate_if_no_never_treated()`. Defaults to `TRUE`.
#' @param estimator_type Character scalar; one of `"fetwfe"` or
#'   `"etwfe"`. Selects which validator to call. `"fetwfe"` covers
#'   both `fetwfe()` and `betwfe()` (they share `checkFetwfeInputs()`);
#'   `"etwfe"` covers both `etwfe()` and `twfeCovs()` (they share
#'   `checkEtwfeInputs()`).
#' @return A list with exactly 14 elements:
#'   `pdata`, `covs`, `X_ints`, `y`, `y_mean`, `N`, `T`, `d`, `p`,
#'   `in_sample_counts`, `num_treats`, `first_inds`, `R`,
#'   `indep_count_data_available`. Each caller unpacks these into local
#'   variables before invoking its `_core()` function.
#' @keywords internal
#' @noRd
.run_estimator_input_prep <- function(
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
	q = NA,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	estimator_type
) {
	stopifnot(
		is.character(estimator_type),
		length(estimator_type) == 1,
		estimator_type %in% c("fetwfe", "etwfe")
	)

	# Defensive caller-error guard: etwfe/twfeCovs callers must not
	# forward bridge-penalty arguments.
	if (estimator_type == "etwfe") {
		stopifnot(is.na(q) && is.na(lambda.max) && is.na(lambda.min))
	}

	# Step 3: input validation.
	if (estimator_type == "fetwfe") {
		ret <- checkFetwfeInputs(
			pdata = pdata,
			time_var = time_var,
			unit_var = unit_var,
			treatment = treatment,
			response = response,
			covs = covs,
			indep_counts = indep_counts,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			lambda.max = lambda.max,
			lambda.min = lambda.min,
			q = q,
			verbose = verbose,
			alpha = alpha,
			add_ridge = add_ridge
		)
	} else {
		ret <- checkEtwfeInputs(
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
	}

	pdata <- ret$pdata
	indep_count_data_available <- ret$indep_count_data_available

	# Step 4: auto-truncation when no never-treated units.
	pdata <- .truncate_if_no_never_treated(
		pdata,
		time_var = time_var,
		unit_var = unit_var,
		treat_var = treatment,
		allow_no_never_treated = allow_no_never_treated
	)

	# Step 5: design-matrix prep.
	res1 <- prep_for_etwfe_core(
		pdata = pdata,
		response = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs,
		verbose = verbose,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts
	)

	list(
		pdata = res1$pdata,
		covs = res1$covs,
		X_ints = res1$X_ints,
		y = res1$y,
		y_mean = res1$y_mean,
		N = res1$N,
		T = res1$T,
		d = res1$d,
		p = res1$p,
		in_sample_counts = res1$in_sample_counts,
		num_treats = res1$num_treats,
		first_inds = res1$first_inds,
		R = res1$R,
		indep_count_data_available = indep_count_data_available
	)
}

#' @title Select the ATT/SE/cohort-probs branch from a `_core()` result
#' @description
#' Internal helper that consolidates the verbatim Step 8 ("pick from
#' in-sample vs indep") shared across the four estimator entry points.
#' Reads `res$indep_*` vs `res$in_sample_*` based on
#' `indep_count_data_available`. The SE-non-NA `stopifnot()` guards
#' that the four pre-consolidation entry points carried differ in
#' shape between the bridge and OLS estimators, so the helper takes a
#' `q` argument to discriminate. For fetwfe/betwfe callers (non-NA
#' `q`), the SE-non-NA guard runs only when `q < 1 && res$calc_ses` in
#' both branches. For etwfe/twfeCovs callers (`q == NA`), the SE-non-NA
#' guard on `res$indep_att_se` runs unconditionally in the indep
#' branch (matching pre-#79 etwfe/twfeCovs behavior: rank-deficient
#' fits with `indep_counts` still trip this guard); the in-sample
#' branch carries no SE-non-NA guard for these callers. Consolidated
#' by GitHub #79.
#' @param res A list returned by `fetwfe_core()`, `etwfe_core()`,
#'   `betwfe_core()`, or `twfeCovs_core()`. Required fields are
#'   `att_hat`, `att_se`, `cohort_probs`, `cohort_probs_overall`,
#'   `indep_att_hat`, `indep_att_se`, `indep_cohort_probs`,
#'   `indep_cohort_probs_overall`, `in_sample_att_hat`,
#'   `in_sample_att_se`, and `calc_ses`.
#' @param indep_count_data_available Logical scalar; if `TRUE`, the
#'   indep_* branch is selected.
#' @param q Numeric scalar; the bridge penalty exponent. `NA` (default)
#'   for etwfe/twfeCovs callers; non-NA for fetwfe/betwfe callers.
#'   Selects which SE-non-NA guards run (see Description).
#' @return A list with four elements: `att_hat`, `att_se`,
#'   `cohort_probs`, `cohort_probs_overall`.
#' @keywords internal
#' @noRd
.select_att_branch <- function(
	res,
	indep_count_data_available,
	q = NA
) {
	uses_bridge <- !is.na(q)

	if (indep_count_data_available) {
		stopifnot(!is.na(res$indep_att_hat))

		if (uses_bridge) {
			if (q < 1 && res$calc_ses) {
				stopifnot(!is.na(res$indep_att_se))
			}
		} else {
			# Issue #84 item 8: gate the etwfe/twfeCovs unconditional
			# guard on `calc_ses`. The pre-#79 behavior stopped()-d on
			# `is.na(res$indep_att_se)` even when SEs had been
			# intentionally not computed (e.g., a future code path
			# where calc_ses = FALSE flows through this branch with
			# valid indep counts). Currently `calc_ses = FALSE` is
			# only reachable via the `q >= 1` bridge regime, so no
			# live code path can trip the un-gated guard today — but
			# the explicit gate hardens against a future refactor that
			# routes a calc_ses = FALSE etwfe/twfeCovs fit through
			# this branch.
			if (res$calc_ses) {
				stopifnot(!is.na(res$indep_att_se))
			}
		}

		stopifnot(all(!is.na(res$indep_cohort_probs)))
		att_hat <- res$indep_att_hat
		att_se <- res$indep_att_se
		cohort_probs <- res$indep_cohort_probs
		cohort_probs_overall <- res$indep_cohort_probs_overall
	} else {
		if (uses_bridge) {
			stopifnot(!is.na(res$in_sample_att_hat))

			if (q < 1 && res$calc_ses) {
				stopifnot(!is.na(res$in_sample_att_se))
			}
		}

		att_hat <- res$in_sample_att_hat
		att_se <- res$in_sample_att_se
		cohort_probs <- res$cohort_probs
		cohort_probs_overall <- res$cohort_probs_overall
	}

	list(
		att_hat = att_hat,
		att_se = att_se,
		cohort_probs = cohort_probs,
		cohort_probs_overall = cohort_probs_overall
	)
}

#' @title 4-way grpreg::gBridge dispatch with lambda-path diagnostics
#' @description Dispatches `grpreg::gBridge()` based on which of `lambda.max`
#'   and `lambda.min` were user-supplied (NA -> leave as default; non-NA ->
#'   pass through). Returns the `fit` object plus the four lambda-path
#'   diagnostic locals (`lambda.max`, `lambda.min`, `lambda.max_model_size`,
#'   `lambda.min_model_size`) used by `fetwfe_core()` and `betwfe_core()`
#'   downstream.
#'
#'   Extracted from a byte-identical 46-line block previously duplicated
#'   across `R/fetwfe_core.R` and `R/betwfe_core.R`. The duplication was
#'   a documented drift risk (issue #119); future fixes to the lambda
#'   path now have a single landing site.
#' @param X_final_scaled Numeric matrix; the design matrix (typically
#'   GLS-then-fusion transformed and scaled).
#' @param y_final Numeric vector; the GLS-transformed response.
#' @param q Numeric scalar; the bridge-penalty exponent passed as `gamma` to
#'   `gBridge()`.
#' @param lambda.max Numeric scalar or `NA`; if non-NA, pass through to
#'   `gBridge()` as the maximum lambda. If NA, `gBridge()` picks one.
#' @param lambda.min Numeric scalar or `NA`; if non-NA, pass through to
#'   `gBridge()`. If NA, `gBridge()` picks one.
#' @param nlambda Integer; the lambda-grid length passed to `gBridge()`.
#' @param verbose Logical; if `TRUE`, emit progress messages around the
#'   `gBridge()` fit.
#' @return A list containing:
#'   \item{fit}{The `gBridge` fit object.}
#'   \item{lambda.max}{`max(fit$lambda)` -- the realised maximum of the
#'     lambda grid.}
#'   \item{lambda.min}{`min(fit$lambda)` -- the realised minimum.}
#'   \item{lambda.max_model_size}{Number of nonzero `fit$beta` columns at
#'     `lambda.max` (largest lambda, smallest model).}
#'   \item{lambda.min_model_size}{Number of nonzero `fit$beta` columns at
#'     `lambda.min` (smallest lambda, largest model).}
#' @keywords internal
#' @noRd
.fit_bridge_with_lambda_path <- function(
	X_final_scaled,
	y_final,
	q,
	lambda.max,
	lambda.min,
	nlambda,
	verbose
) {
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

	list(
		fit = fit,
		lambda.max = max(fit$lambda),
		lambda.min = min(fit$lambda),
		lambda.max_model_size = sum(fit$beta[, ncol(fit$beta)] != 0),
		lambda.min_model_size = sum(fit$beta[, 1] != 0)
	)
}

#' @title Compute treatment-effect and treatment-interaction indices
#' @description Returns the integer index vectors `treat_inds` (base
#'   treatment-effect columns) and `treat_int_inds` (covariate x
#'   treatment interaction columns) for the ETWFE / FETWFE / BETWFE
#'   design matrix.
#'
#'   The math is `treat_inds = getTreatInds(R, T, d, num_treats)` plus
#'   the contract `max(treat_inds) == R + T - 1 + d + R*d + (T-1)*d +
#'   num_treats` when `d > 0` (and `max(treat_inds) == R + T - 1 +
#'   num_treats` when `d == 0`). The treatment-interaction block lives
#'   immediately above `treat_inds` and runs `(max(treat_inds) + 1):p`
#'   when `d > 0`; empty otherwise.
#'
#'   Extracted from a byte-identical 17-line block previously duplicated
#'   across `R/etwfe_core.R`, `R/fetwfe_core.R`, and `R/betwfe_core.R`
#'   (issue #119).
#' @param R Integer; number of treated cohorts.
#' @param T Integer; number of time periods.
#' @param d Integer; number of (time-invariant) covariates.
#' @param num_treats Integer; total number of base treatment-effect
#'   parameters in the design.
#' @param p Integer; total number of columns in the design matrix.
#' @return A list with:
#'   \item{treat_inds}{Integer vector of base treatment-effect column
#'     indices, of length `num_treats`.}
#'   \item{treat_int_inds}{Integer vector of treatment x covariate
#'     interaction column indices, of length `num_treats * d` when
#'     `d > 0` and empty when `d == 0`.}
#' @keywords internal
#' @noRd
.compute_treat_inds <- function(R, T, d, num_treats, p) {
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
	list(treat_inds = treat_inds, treat_int_inds = treat_int_inds)
}
