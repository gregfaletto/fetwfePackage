# Input preparation and panel/design setup for the ETWFE / bridge fitting
# core. Split out of R/core_funcs.R for issue #186 (no behavior change). Holds
# the input validation, design-matrix prep, and cohort-probability helpers that
# run before the GLS step (see R/gls_machinery.R).

#' @title Validate that `in_sample_counts` is fully named with unique names
#' @description `prep_for_etwfe_core()` and `check_etwfe_core_inputs()` both
#'   require `in_sample_counts` to have one uniquely-named entry per cohort.
#'   Factored out of the two byte-identical validation blocks (#325) so the
#'   error message lives in a single place.
#' @param in_sample_counts The (named) count vector to validate.
#' @return `invisible(NULL)`; errors via `stop()` if the names are missing or
#'   non-unique.
#' @keywords internal
#' @noRd
.validate_in_sample_counts_named_unique <- function(in_sample_counts) {
	nm <- names(in_sample_counts)
	if (
		length(nm) != length(in_sample_counts) ||
			length(nm) != length(unique(nm))
	) {
		stop(
			"in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)",
			call. = FALSE
		)
	}
	invisible(NULL)
}

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
#' @param N,T,G,d,p,num_treats Integers giving key problem dimensions:
#'   number of units, time periods, treated cohorts, covariates,
#'   total parameters, and base treatment-effect parameters, respectively.
#' @param add_ridge Logical.  Whether to append rows that implement a small
#'   L2 penalty on the \emph{untransformed} coefficients.  Default \code{FALSE}.
#' @param first_inds Integer vector (length \(G\)).  Column indices of the
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
#'   are constructed).  Required (no default): the prior default of
#'   \code{TRUE} produced GitHub #74 (BETWFE silently inherited the FETWFE
#'   augmentation basis).  All four caller sites
#'   (\code{fetwfe_core}, \code{etwfe_core}, \code{betwfe_core},
#'   \code{twfeCovs_core}) must pass this argument explicitly.
#' @param is_twfe_covs Logical.  If \code{TRUE}, columns will be removed and
#'   consolidated as required for `twfeCovs()`.  Default \code{FALSE}.
#' @param fusion_structure Character; one of `"cohort"` (default) or
#'   `"event_study"`. Forwarded to `.append_ridge_rows()` (and thence the
#'   full inverse-fusion matrix) when `add_ridge = TRUE && is_fetwfe`.
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   user-supplied already-inverted treatment-effect fusion block
#'   (`solve(fusion_matrix)`, #236), forwarded to `.append_ridge_rows()`.
#'   Default `NULL`.
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
#'     Computes \eqn{\hat\pi_g = n_g / \sum_{s=1}^G n_s} from the in-sample
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
#'     \item{\code{cohort_probs}}{Vector of \(\hat\pi_g \mid \text{treated}\)
#'       from the estimation sample.}
#'     \item{\code{cohort_probs_overall}}{Vector of unconditional
#'       probabilities \(P(W=g)\).}
#'     \item{\code{indep_cohort_probs}, \code{indep_cohort_probs_overall}}{Same
#'       probabilities if an independent split was provided; otherwise \code{NA}.}
#'     \item{\code{sig_eps_sq}, \code{sig_eps_c_sq}}{Possibly estimated
#'       variance components carried forward.}
#'     \item{\code{lambda_ridge}}{Numeric scalar.  Value of the ridge penalty
#'       used (or \code{NA} if none).}
#'     \item{\code{twfe_covs_treat_inds}}{Integer vector of treatment-column
#'       indices in the collapsed twfeCovs design (from
#'       \code{.collapse_design_for_twfe_covs()}); \code{NULL} on every other
#'       path (#337).}
#'   }
#' @keywords internal
#' @noRd
prep_for_etwfe_regression <- function(
	verbose,
	sig_eps_sq,
	sig_eps_c_sq,
	y,
	X_ints,
	X_mod,
	N,
	T,
	G,
	d,
	p,
	num_treats,
	add_ridge,
	first_inds,
	in_sample_counts,
	indep_count_data_available,
	indep_counts = NA,
	is_fetwfe,
	is_twfe_covs = FALSE,
	fusion_structure = "cohort",
	d_inv_treat = NULL,
	gls = TRUE
) {
	if (gls) {
		# Default: estimate the variance components (REML when not supplied) and
		# GLS-whiten the design. This is where `estOmegaSqrtInv()` hard-stops at
		# `p >= N(T - 1)` (its REML guard) -- the reason `gls = FALSE` exists for the
		# high-dimensional path.
		gls_res <- .estimate_variance_and_gls(
			y = y,
			X_ints = X_ints,
			X_mod = X_mod,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			N = N,
			T = T,
			p = p,
			verbose = verbose
		)
		y_final <- gls_res$y_gls
		X_final <- gls_res$X_gls
		sig_eps_sq <- gls_res$sig_eps_sq
		sig_eps_c_sq <- gls_res$sig_eps_c_sq
	} else {
		# gls = FALSE (#307): skip whitening AND variance-component estimation
		# entirely. The high-dimensional (`p >= NT`) path, where REML cannot estimate
		# Omega -- and whitening buys efficiency, not validity: the `debiasedATT()`
		# cluster-robust sandwich SE needs no Omega (paper Decision D1 /
		# gregfaletto/fetwfe#90). Use the un-whitened (fusion-transformed) design
		# `X_mod`; leave the variance components NA (`calc_ses` is forced FALSE in
		# `fetwfe_core()`).
		y_final <- as.numeric(y)
		X_final <- X_mod
		sig_eps_sq <- NA_real_
		sig_eps_c_sq <- NA_real_
	}

	# Treatment-column indices in the (possibly collapsed) design. Only the
	# twfeCovs collapse re-lays-out the columns; every other path derives its own
	# treatment indices, so this stays NULL there. (#337)
	twfe_covs_treat_inds <- NULL
	if (is_twfe_covs) {
		coll <- .collapse_design_for_twfe_covs(
			X_gls = X_final,
			N = N,
			T = T,
			G = G,
			d = d,
			num_treats = num_treats,
			first_inds = first_inds
		)
		X_final <- coll$X_collapsed
		p <- coll$p_short
		twfe_covs_treat_inds <- coll$treat_inds
	}

	X_final_scaled <- my_scale(X_final)
	stopifnot(ncol(X_final_scaled) == p)
	scale_center <- attr(X_final_scaled, "scaled:center")
	scale_scale <- attr(X_final_scaled, "scaled:scale")

	ridge <- .append_ridge_rows(
		X_scaled = X_final_scaled,
		y_gls = y_final,
		p = p,
		add_ridge = add_ridge,
		is_fetwfe = is_fetwfe,
		first_inds = first_inds,
		T = T,
		G = G,
		d = d,
		num_treats = num_treats,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		fusion_structure = fusion_structure,
		d_inv_treat = d_inv_treat
	)
	X_final_scaled <- ridge$X_scaled
	y_final <- ridge$y_final
	lambda_ridge <- ridge$lambda_ridge

	probs <- .compute_cohort_probs(
		in_sample_counts = in_sample_counts,
		indep_counts = indep_counts,
		N = N,
		G = G,
		indep_count_data_available = indep_count_data_available
	)

	list(
		X_final_scaled = X_final_scaled,
		X_final = X_final,
		y_final = y_final,
		scale_center = scale_center,
		scale_scale = scale_scale,
		cohort_probs = probs$cohort_probs,
		cohort_probs_overall = probs$cohort_probs_overall,
		indep_cohort_probs = probs$indep_cohort_probs,
		indep_cohort_probs_overall = probs$indep_cohort_probs_overall,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda_ridge = lambda_ridge,
		twfe_covs_treat_inds = twfe_covs_treat_inds
	)
}


#' Prepare Data & Design Matrix for ETWFE/TWFE Workflows
#'
#' @description
#' A helper that converts raw **panel data** into the core objects required by
#' `etwfe_core()`, `twfeCovs_core()`, and related fitting routines.  The function
#' (i) keeps only the relevant variables, (ii) one-hot-encodes or otherwise
#' expands any *factor* covariates, (iii) builds the stacked design matrix
#' `X_ints` and centred response `y` via \code{prepXints()}, and (iv) performs a
#' battery of **sanity checks** on cohort counts before estimation proceeds.
#'
#' @param pdata A **data.frame** in long (unit x time) format containing at
#'   least the columns named in `response`, `time_var`, `unit_var`,
#'   `treatment`, and `covs`.
#' @param response Character; name of the response variable column.
#' @param time_var Character; name of the time variable column.
#' @param unit_var Character; name of the unit variable column.
#' @param treatment Character; name of the treatment indicator column.
#' @param covs Character vector of additional covariate column names.  Factor
#'   covariates are expanded to dummies by `processFactors()`.  May be empty.
#' @param verbose Logical; print progress/timing information?
#' @param indep_count_data_available Logical; `TRUE` when *independent* cohort
#'   counts are supplied in `indep_counts` and should be validated.
#' @param indep_counts Optional integer vector with length
#'   `1 + G` (never-treated plus `G` treated cohorts) giving cohort sizes in an
#'   independent sample.  Only used when
#'   `indep_count_data_available = TRUE`.
#'
#' @return A named **list** ready to be passed into the "_core" estimators:
#' \describe{
#'   \item{pdata}{The (possibly modified) data frame after factor processing.}
#'   \item{covs}{Updated character vector of covariate names after dummy-expansion.}
#'   \item{X_ints}{Design matrix with unit FE, time FE, covariates,
#'     treatment dummies and their interactions (dimensions \(N T \times p\)).}
#'   \item{y}{Centred response vector of length \(N T\).}
#'   \item{y_mean}{Mean of the original (pre-centering) response. Preserved so
#'     downstream methods (`augment()`, `predict()`) can return fitted values
#'     on the original response scale.}
#'   \item{N, T}{Integers - number of unique units and time periods.}
#'   \item{d}{Integer - number of *raw* covariates after processing.}
#'   \item{p}{Integer - total number of columns in `X_ints`.}
#'   \item{in_sample_counts}{Named integer vector of length `1 + G` with cohort
#'     sizes in the estimation sample (first entry = never-treated).}
#'   \item{num_treats}{Integer - total number of base treatment-effect
#'     parameters in `X_ints`.}
#'   \item{first_inds}{Integer vector (length `G`) giving the first column
#'     index of each cohort's treatment-effect block inside `X_ints`.}
#'   \item{G}{Integer - number of treated cohorts detected.}
#' }
#'
#' @details
#' The routine executes the following steps:
#' \enumerate{
#'   \item **Column subset** - drops all variables not explicitly required.
#'   \item **Factor processing** - calls `processFactors()` so that every
#'     categorical covariate becomes a set of \(0/1\) dummies.
#'   \item **Design-matrix construction** - hands the cleaned data to
#'     `prepXints()`, which adds unit and time dummies, interactions, etc.
#'   \item **Cohort diagnostics** - verifies that
#'     \eqn{G \ge 1}, at least one never-treated unit exists, cohort names are
#'     unique, and (if provided) `indep_counts` are dimensionally consistent.
#' }
#'
#' Any violation of the checks triggers an informative `stop()` message so that
#' higher-level functions fail fast.
#'
#' @keywords internal
#' @noRd
prep_for_etwfe_core <- function(
	pdata,
	response,
	time_var,
	unit_var,
	treatment,
	covs,
	verbose,
	indep_count_data_available,
	indep_counts
) {
	# Subset pdata to include only the key columns
	pdata <- pdata[, c(response, time_var, unit_var, treatment, covs)]

	# Process any factor covariates:
	if (length(covs) > 0) {
		pf_res <- processFactors(pdata, covs)
		pdata <- pf_res$pdata
		covs <- pf_res$covs
	}

	# Proceed to generate design matrix and other objects
	res <- prepXints(
		data = pdata,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs,
		response = response,
		verbose = verbose
	)

	X_ints <- res$X_ints
	y <- res$y
	y_mean <- res$y_mean
	N <- res$N
	T <- res$T
	d <- res$d
	p <- res$p
	in_sample_counts <- res$in_sample_counts
	num_treats <- res$num_treats
	first_inds <- res$first_inds
	first_year <- res$first_year

	rm(res)

	G <- length(in_sample_counts) - 1
	stopifnot(G <= T - 1)
	if (G < 1) {
		stop(
			"No treated cohorts detected in data. fetwfe, etwfe, betwfe, and twfeCovs require at least one treated cohort (one cohort of units that adopt treatment, plus never-treated units)."
		)
	}
	stopifnot(N >= G + 1)
	stopifnot(sum(in_sample_counts) == N)
	stopifnot(all(in_sample_counts >= 0))
	stopifnot(is.integer(in_sample_counts))
	if (in_sample_counts[1] == 0) {
		stop(
			"No never-treated units detected in data to fit model; estimating treatment effects is not possible"
		)
	}
	.validate_in_sample_counts_named_unique(in_sample_counts)
	if (indep_count_data_available) {
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
	}

	return(list(
		pdata = pdata,
		covs = covs,
		X_ints = X_ints,
		y = y,
		y_mean = y_mean,
		N = N,
		T = T,
		d = d,
		p = p,
		in_sample_counts = in_sample_counts,
		num_treats = num_treats,
		first_inds = first_inds,
		first_year = first_year,
		G = G
	))
}

#' @title Collapse the twfeCovs design matrix
#' @description Drop treatment-effect interaction columns and the
#'   covariate-by-time / covariate-by-cohort interaction columns, then
#'   collapse the per-cohort treatment dummies into one column per
#'   cohort. Called from the orchestrator only when
#'   `is_twfe_covs = TRUE`.
#' @param X_gls Numeric matrix; GLS-transformed design.
#' @param N,T,G,d Integers; units, time periods, treated cohorts,
#'   covariates.
#' @param num_treats Integer; treatment-column count in the pre-collapse
#'   design.
#' @param first_inds Integer vector; first-index-within-cohort offsets.
#' @return List with `X_collapsed` (post-collapse design), `p_short`
#'   (column count `= G + T - 1 + d + G`), and `treat_inds` (indices of the
#'   trailing G treatment columns in `X_collapsed`). `treat_inds` is the single
#'   source of truth for the collapsed-design treatment-column layout, consumed
#'   by `.ols_estimator_core()` rather than re-derived there (#337).
#' @keywords internal
#' @noRd
.collapse_design_for_twfe_covs <- function(
	X_gls,
	N,
	T,
	G,
	d,
	num_treats,
	first_inds
) {
	stopifnot(nrow(X_gls) == N * T)

	X <- X_gls[, 1:(G + T - 1 + d * (1 + G + T - 1) + num_treats)]

	# Treatment columns in the (pre-collapse) GLS design. Use the validated helper
	# getTreatInds() rather than an inline formula, so the pre-collapse treatment
	# block matches the etwfe path exactly and cannot drift. This fixes a previous
	# off-by-one: the old inline `first_treat_ind <- G + T - 1 + d*(1 + G + T - 1)`
	# was the *count* of pre-treatment columns, so `first_treat_ind:(...)` started
	# one column too low -- it pulled the last covariate-interaction column into the
	# treatment block and dropped the last treatment column (#339).
	treat_inds <- getTreatInds(G = G, T = T, d = d, num_treats = num_treats)
	stopifnot(length(treat_inds) == num_treats)

	X <- X[, c(1:(G + T - 1 + d), treat_inds)]
	stopifnot(nrow(X) == N * T)

	treat_inds_mat <- matrix(as.numeric(NA), nrow = N * T, ncol = G)
	for (g in 1:G) {
		inds_g <- .cohort_block_inds(g, G, first_inds, num_treats)
		cols_g <- G + T - 1 + d + inds_g
		treat_inds_mat[, g] <- rowSums(X[, cols_g, drop = FALSE])
	}
	stopifnot(all(!is.na(treat_inds_mat)))

	# The collapsed per-cohort treatment columns are cbind'd AFTER the
	# non-treatment block (cohort FE + time FE + covariate main effects), so the
	# treatment block is the trailing G columns. Derive treat_inds from that
	# actual placement and return it as the single source of truth, so the
	# "treatment = trailing G columns" contract is not independently re-encoded by
	# callers such as .ols_estimator_core() (#337).
	n_non_treat <- G + T - 1 + d
	X_collapsed <- cbind(X[, 1:n_non_treat], treat_inds_mat)

	p_short <- ncol(X_collapsed)
	treat_inds_collapsed <- (n_non_treat + 1):p_short

	# p_short == n_non_treat + G already implies length(treat_inds_collapsed) == G.
	stopifnot(p_short == n_non_treat + G)

	list(
		X_collapsed = X_collapsed,
		p_short = p_short,
		treat_inds = treat_inds_collapsed
	)
}

#' @title Compute in-sample (and optionally independent) cohort probabilities
#' @description Convert raw cohort counts into two normalizations:
#'   within-treated (`cohort_probs`, sums to 1) and overall
#'   (`cohort_probs_overall`, equals `count_g / N`). If
#'   `indep_count_data_available = TRUE`, do the same for `indep_counts`;
#'   otherwise the indep-side outputs are NA.
#' @param in_sample_counts Integer vector of length `G + 1`; counts of
#'   never-treated + treated cohorts in the working panel.
#' @param indep_counts Integer vector of length `G + 1` or NA;
#'   independent-sample counts.
#' @param N Integer; total units.
#' @param G Integer; treated cohorts.
#' @param indep_count_data_available Logical.
#' @return List with `cohort_probs`, `cohort_probs_overall`,
#'   `indep_cohort_probs`, `indep_cohort_probs_overall`.
#' @keywords internal
#' @noRd
.compute_cohort_probs <- function(
	in_sample_counts,
	indep_counts,
	N,
	G,
	indep_count_data_available
) {
	cohort_probs <- in_sample_counts[2:(G + 1)] /
		sum(in_sample_counts[2:(G + 1)])

	stopifnot(all(!is.na(cohort_probs)))
	stopifnot(all(cohort_probs >= 0))
	stopifnot(all(cohort_probs <= 1))
	stopifnot(length(cohort_probs) == G)
	stopifnot(abs(sum(cohort_probs) - 1) < 1e-6)

	cohort_probs_overall <- in_sample_counts[2:(G + 1)] / N

	stopifnot(
		abs(1 - sum(cohort_probs_overall) - in_sample_counts[1] / N) < 1e-6
	)

	if (indep_count_data_available) {
		indep_cohort_probs <- indep_counts[2:(G + 1)] /
			sum(indep_counts[2:(G + 1)])

		stopifnot(all(!is.na(indep_cohort_probs)))
		stopifnot(all(indep_cohort_probs >= 0))
		stopifnot(all(indep_cohort_probs <= 1))
		stopifnot(length(indep_cohort_probs) == G)
		stopifnot(abs(sum(indep_cohort_probs) - 1) < 1e-6)

		indep_cohort_probs_overall <- indep_counts[2:(G + 1)] / N

		stopifnot(
			abs(
				1 -
					sum(indep_cohort_probs_overall) -
					indep_counts[1] / N
			) <
				1e-6
		)
	} else {
		indep_cohort_probs <- NA
		indep_cohort_probs_overall <- NA
	}

	list(
		cohort_probs = cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs = indep_cohort_probs,
		indep_cohort_probs_overall = indep_cohort_probs_overall
	)
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
#' @param in_sample_counts Integer named vector. Length `G+1`.
#'   The first element must be the number of never-treated units; the
#'   remaining `G` elements give the number of units in each treated
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
#'   \item Confirming that \eqn{1 \le G \le T-1}.
#'   \item Basic sanity checks for all numeric scalars (\code{alpha} inside \eqn{(0,1)}, non-negative variances, etc.).
#'   \item All structural requirements on \code{indep_counts} when supplied.
#' }
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{G}}{Integer. The number of treated cohorts `=length(in_sample_counts)-1`).}
#'     \item{\code{c_names}}{Character vector of cohort names (length `G`).}
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
	G <- length(in_sample_counts) - 1

	c_names <- names(in_sample_counts)[2:(G + 1)]

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
	.validate_in_sample_counts_named_unique(in_sample_counts)

	# Mirror prep_for_etwfe_core()'s friendly message (#208) so that whichever
	# validator fires first, the user sees the same actionable text rather than
	# an opaque `stopifnot` failure.
	if (G < 1) {
		stop(
			"No treated cohorts detected in data. fetwfe, etwfe, betwfe, and twfeCovs require at least one treated cohort (one cohort of units that adopt treatment, plus never-treated units)."
		)
	}
	stopifnot(G <= T - 1)

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
		G = G,
		c_names = c_names,
		indep_count_data_available = indep_count_data_available
	))
}
