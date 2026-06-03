# Data simulation pipeline: `simulateData()` and its core
# (`simulateDataCore()`). Moved from R/gen_funcs.R in 1.9.25.

#' Generate Random Panel Data for FETWFE Simulations
#'
#' @description
#' Generates a random panel data set for simulation studies of the fused extended two-way fixed
#' effects (FETWFE) estimator by taking an object of class  \code{"FETWFE_coefs"} (produced by
#' \code{genCoefs()}) and using it to simulate data. The function creates a balanced panel
#' with \eqn{N} units over \eqn{T} time periods, assigns treatment status across \eqn{G}
#' treated cohorts (with equal marginal probabilities for treatment and non-treatment), and
#' constructs a design matrix along with the corresponding outcome. The covariates are
#' generated according to the specified \code{distribution}: by default, covariates are drawn
#' from a normal distribution; if \code{distribution = "uniform"}, they are drawn uniformly
#' from \eqn{[-\sqrt{3}, \sqrt{3}]}. When \eqn{d = 0} (i.e. no covariates), no
#' covariate-related columns or interactions are generated. See the simulation studies section of
#' Faletto (2025) for details.
#'
#' @param coefs_obj An object of class \code{"FETWFE_coefs"} containing the coefficient vector
#' and simulation parameters.
#' @param N Integer. Number of units in the panel.
#' @param sig_eps_sq Numeric. Variance of the idiosyncratic (observation-level) noise.
#' @param sig_eps_c_sq Numeric. Variance of the unit-level random effects.
#'   Must be non-negative; `0` is allowed (yields a panel with no unit-level
#'   random effects).
#' @param distribution Character. Distribution to generate covariates.
#'   Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
#'   from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#' @param guarantee_rank_condition (Optional). Logical. If TRUE, the returned
#' data set is guaranteed to have at least `d + 1` units per cohort, which is
#' necessary for the final design matrix to have full column rank. Default is
#' FALSE, in which case no such condition is enforced.
#'
#' @return An object of class \code{"FETWFE_simulated"}, which is a list containing:
#' \describe{
#'   \item{pdata}{A dataframe containing generated data that can be passed to \code{fetwfe()}.}
#'   \item{X}{The design matrix \eqn{X}, with \eqn{p} columns with interactions.}
#'   \item{y}{A numeric vector of length \eqn{N \times T} containing the generated responses.}
#'   \item{covs}{A character vector containing the names of the generated features (if \eqn{d > 0}),
#'          or simply an empty vector (if \eqn{d = 0})}
#'   \item{time_var}{The name of the time variable in pdata}
#'   \item{unit_var}{The name of the unit variable in pdata}
#'   \item{treatment}{The name of the treatment variable in pdata}
#'   \item{response}{The name of the response variable in pdata}
#'   \item{coefs}{The coefficient vector \eqn{\beta} used for data generation.}
#'   \item{first_inds}{A vector of indices indicating the first treatment effect for each treated
#'         cohort.}
#'   \item{N_UNTREATED}{The number of never-treated units.}
#'   \item{assignments}{A vector of counts (of length \eqn{G+1}) indicating how many units fall into
#'         the never-treated group and each of the \eqn{G} treated cohorts.}
#'   \item{indep_counts}{Independent cohort assignments (for auxiliary purposes).}
#'   \item{p}{The number of columns in the design matrix \eqn{X}.}
#'   \item{N}{Number of units.}
#'   \item{T}{Number of time periods.}
#'   \item{G}{Number of treated cohorts.}
#'   \item{R}{Deprecated alias for \code{G}, retained for backward
#'         compatibility; populated with the same value. Use \code{G}. Will be
#'         removed in a future release.}
#'   \item{d}{Number of covariates.}
#'   \item{sig_eps_sq}{The idiosyncratic noise variance.}
#'   \item{sig_eps_c_sq}{The unit-level noise variance.}
#' }
#'
#' @details
#' This function extracts simulation parameters from the \code{FETWFE_coefs} object and passes them,
#' along with additional simulation parameters, to the internal function \code{simulateDataCore()}.
#' It validates that all necessary components are returned and assigns the S3 class
#' \code{"FETWFE_simulated"} to the output.
#'
#' The simulated panel is fully determined by \code{coefs_obj$seed} (the seed
#' passed to \code{genCoefs()}): \code{simulateData()} re-seeds from it on every
#' call, so re-calling \code{simulateData()} on the same \code{coefs_obj}
#' returns the identical panel. To vary the panel across Monte Carlo
#' replications, vary the \code{seed} passed to \code{genCoefs()} (re-calling
#' \code{simulateData()} alone does not).
#'
#' The argument \code{distribution} controls the generation of covariates. For
#' \code{"gaussian"}, covariates are drawn from \code{rnorm}; for \code{"uniform"},
#' they are drawn from \code{runif} on the interval \eqn{[-\sqrt{3}, \sqrt{3}]} (which ensures that
#' the covariates have unit variance regardless of which distribution is chosen).
#'
#' When \eqn{d = 0} (i.e. no covariates), the function omits any covariate-related columns
#' and their interactions.
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
#'   coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#' }
#'
#' @export
simulateData <- function(
	coefs_obj,
	N,
	sig_eps_sq,
	sig_eps_c_sq,
	distribution = "gaussian",
	guarantee_rank_condition = FALSE
) {
	if (!inherits(coefs_obj, "FETWFE_coefs")) {
		stop("coefs_obj must be an object of class 'FETWFE_coefs'")
	}
	if (!is.numeric(N) || length(N) != 1 || N <= 0) {
		stop("N must be a positive numeric value")
	}
	if (!is.numeric(sig_eps_sq) || length(sig_eps_sq) != 1 || sig_eps_sq <= 0) {
		stop("sig_eps_sq must be a positive numeric value")
	}
	if (
		!is.numeric(sig_eps_c_sq) ||
			length(sig_eps_c_sq) != 1 ||
			sig_eps_c_sq < 0
	) {
		stop("sig_eps_c_sq must be a non-negative numeric value")
	}

	# Extract parameters from the coefs object
	G <- coefs_obj$G
	T <- coefs_obj$T
	d <- coefs_obj$d
	beta <- coefs_obj$beta
	seed <- coefs_obj$seed
	# Backward-compat: saved-from-v1.13.x FETWFE_coefs objects don't have
	# these slots; coerce a missing $assignment_type to "marginal" so the
	# upgrade path doesn't crash on match.arg(NULL, c(...)) downstream.
	assignment_type <- if (is.null(coefs_obj$assignment_type)) {
		"marginal"
	} else {
		coefs_obj$assignment_type
	}
	assignment_coefs <- coefs_obj$assignment_coefs

	sim_data <- simulateDataCore(
		N = N,
		T = T,
		G = G,
		d = d,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		beta = beta,
		seed = seed,
		gen_ints = TRUE,
		distribution = distribution,
		guarantee_rank_condition = guarantee_rank_condition,
		assignment_type = assignment_type,
		assignment_coefs = assignment_coefs
	)

	required_fields <- c(
		"pdata",
		"X",
		"y",
		"covs",
		"time_var",
		"unit_var",
		"treatment",
		"response",
		"coefs",
		"first_inds",
		"N_UNTREATED",
		"assignments",
		"indep_counts",
		"p",
		"N",
		"T",
		"G",
		"R",
		"d",
		"sig_eps_sq",
		"sig_eps_c_sq"
	)
	missing_fields <- setdiff(required_fields, names(sim_data))
	if (length(missing_fields) > 0) {
		stop(paste(
			"simulateDataCore did not return expected components:",
			paste(missing_fields, collapse = ", ")
		))
	}

	if (!inherits(sim_data, "FETWFE_simulated")) {
		stop("sim_data must be an object of class 'FETWFE_simulated'")
	}

	return(sim_data)
}


#' Generate Random Panel Data for FETWFE Simulations (core)
#'
#' @description
#' Generates a random panel data set for simulation studies of the fused extended two-way fixed
#' effects (FETWFE) estimator. The function creates a balanced panel with \eqn{N} units over \eqn{T}
#' time periods, assigns treatment status across \eqn{G} treated cohorts (with equal marginal
#' probabilities for treatment and non-treatment), and constructs a design matrix along with the
#' corresponding outcome. When \code{gen_ints = TRUE} the full design matrix is returned (including
#' interactions between covariates and fixed effects and treatment indicators). When
#' \code{gen_ints = FALSE} the design matrix is generated in a simpler format (with no interactions)
#' as expected by \code{fetwfe()}. Moreover, the covariates are generated according to the
#' specified \code{distribution}: by default, covariates are drawn from a normal distribution;
#' if \code{distribution = "uniform"}, they are drawn uniformly from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' When \eqn{d = 0} (i.e. no covariates), no covariate-related columns or interactions are
#' generated.
#'
#' See the simulation studies section of Faletto (2025) for details.
#'
#' @param N Integer. Number of units in the panel.
#' @param T Integer. Number of time periods.
#' @param G Integer. Number of treated cohorts (with treatment starting in periods 2 to T).
#' @param d Integer. Number of time-invariant covariates.
#' @param sig_eps_sq Numeric. Variance of the idiosyncratic (observation-level) noise.
#' @param sig_eps_c_sq Numeric. Variance of the unit-level random effects.
#'   Must be non-negative; `0` is allowed (yields a panel with no unit-level
#'   random effects).
#' @param beta Numeric vector. Coefficient vector for data generation. Its required length depends
#' on the value of \code{gen_ints}:
#'   \itemize{
#'     \item If \code{gen_ints = TRUE} and \code{d > 0}, the expected length is
#'       \eqn{p = G + (T-1) + d + dG + d(T-1) + num\_treats + num\_treats \times d}, where
#'       \eqn{num\_treats = T \times G - \frac{G(G+1)}{2}}.
#'     \item If \code{gen_ints = TRUE} and \code{d = 0}, the expected length is
#'       \eqn{p = G + (T-1) + num\_treats}.
#'     \item If \code{gen_ints = FALSE}, the expected length is
#'       \eqn{p = G + (T-1) + d + num\_treats}.
#'   }
#' @param seed (Optional) Integer. Seed for reproducibility.
#' @param gen_ints Logical. If \code{TRUE}, generate the full design matrix with interactions;
#'   if \code{FALSE} (the default), generate a design matrix without any interaction terms.
#' @param distribution Character. Distribution to generate covariates.
#'   Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
#'   from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#' @param guarantee_rank_condition (Optional). Logical. If TRUE, the returned
#' data set is guaranteed to have at least `d + 1` units per cohort, which is
#' necessary for the final design matrix to have full column rank. Default is
#' FALSE, in which case no such condition is enforced.
#' @param assignment_type Character. One of \code{"marginal"} (default),
#'   \code{"multinomial"}, or \code{"ordered"}. Selects the cohort-assignment
#'   DGP. \code{"marginal"} preserves the pre-1.14.0 behavior. The
#'   non-marginal types require a non-NULL \code{assignment_coefs} argument
#'   (typically pulled from a \code{FETWFE_coefs} object built with the
#'   matching \code{assignment_type}).
#' @param assignment_coefs Optional list returned by
#'   \code{.gen_assignment_coefs()} (an internal helper). Required when
#'   \code{assignment_type != "marginal"}.
#' @param R Deprecated. The former name for \code{G}; still accepted with a
#'   deprecation warning, and will be removed in a future release. Use
#'   \code{G}.
#'
#' @return An object of class \code{"FETWFE_simulated"}, which is a list containing:
#' \describe{
#'   \item{pdata}{A dataframe containing generated data that can be passed to \code{fetwfe()}.}
#'   \item{X}{The design matrix. When \code{gen_ints = TRUE}, \eqn{X} has \eqn{p} columns with
#'     interactions; when \code{gen_ints = FALSE}, \eqn{X} has no interactions.}
#'   \item{y}{A numeric vector of length \eqn{N \times T} containing the generated responses.}
#'   \item{covs}{A character vector containing the names of the generated features (if \eqn{d > 0}),
#'          or simply an empty vector (if \eqn{d = 0})}
#'   \item{time_var}{The name of the time variable in pdata}
#'   \item{unit_var}{The name of the unit variable in pdata}
#'   \item{treatment}{The name of the treatment variable in pdata}
#'   \item{response}{The name of the response variable in pdata}
#'   \item{coefs}{The coefficient vector \eqn{\beta} used for data generation.}
#'   \item{first_inds}{A vector of indices indicating the first treatment effect for each treated
#'   cohort.}
#'   \item{N_UNTREATED}{The number of never-treated units.}
#'   \item{assignments}{A vector of counts (of length \eqn{G+1}) indicating how many units fall into
#'         the never-treated group and each of the \eqn{G} treated cohorts.}
#'   \item{indep_counts}{Independent cohort assignments (for auxiliary purposes).}
#'   \item{p}{The number of columns in the design matrix \eqn{X}.}
#'   \item{N}{Number of units.}
#'   \item{T}{Number of time periods.}
#'   \item{G}{Number of treated cohorts.}
#'   \item{R}{Deprecated alias for \code{G}, retained for backward
#'         compatibility; populated with the same value. Use \code{G}. Will be
#'         removed in a future release.}
#'   \item{d}{Number of covariates.}
#'   \item{sig_eps_sq}{The idiosyncratic noise variance.}
#'   \item{sig_eps_c_sq}{The unit-level noise variance.}
#' }
#'
#' @details
#' When \code{gen_ints = TRUE}, the function constructs the design matrix by first generating
#' base fixed effects and a long-format covariate matrix (via \code{generateBaseEffects()}), then
#' appending interactions between the covariates and cohort/time fixed effects (via
#' \code{generateFEInts()}) and finally treatment indicator columns and treatment-covariate
#' interactions (via \code{genTreatVarsSim()} and \code{genTreatInts()}). When
#' \code{gen_ints = FALSE}, the design matrix consists only of the base fixed effects, covariates,
#' and treatment indicators.
#'
#' The argument \code{distribution} controls the generation of covariates. For
#' \code{"gaussian"}, covariates are drawn from \code{rnorm}; for \code{"uniform"},
#' they are drawn from \code{runif} on the interval \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' When \eqn{d = 0} (i.e. no covariates), the function omits any covariate-related columns
#' and their interactions.
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' \dontrun{
#'   # Set simulation parameters
#'   N <- 100           # Number of units in the panel
#'   T <- 5             # Number of time periods
#'   G <- 3             # Number of treated cohorts
#'   d <- 2             # Number of time-invariant covariates
#'   sig_eps_sq <- 1    # Variance of observation-level noise
#'   sig_eps_c_sq <- 0.5  # Variance of unit-level random effects
#'
#'   # Generate coefficient vector using genCoefsCore()
#'   # (Here, density controls sparsity and eff_size scales nonzero entries)
#'   coefs_core <- genCoefsCore(G = G, T = T, d = d, density = 0.2, eff_size = 2, seed = 123)
#'
#'   # Now simulate the data. Setting gen_ints = TRUE generates the full design
#'   matrix with interactions.
#'   sim_data <- simulateDataCore(
#'     N = N,
#'     T = T,
#'     G = G,
#'     d = d,
#'     sig_eps_sq = sig_eps_sq,
#'     sig_eps_c_sq = sig_eps_c_sq,
#'     beta = coefs_core$beta,
#'     seed = 456,
#'     gen_ints = TRUE,
#'     distribution = "gaussian"
#'   )
#'
#'   # Examine the returned list:
#'   str(sim_data)
#' }
#'
#' @export
simulateDataCore <- function(
	N,
	T,
	G = NULL,
	d,
	sig_eps_sq,
	sig_eps_c_sq,
	beta,
	seed = NULL,
	gen_ints = FALSE,
	distribution = "gaussian",
	guarantee_rank_condition = FALSE,
	assignment_type = "marginal",
	assignment_coefs = NULL,
	R = NULL
) {
	# Resolve the canonical cohort count (G), mapping the deprecated `R`
	# alias and warning if it is supplied (#41). The body below uses `G`
	# directly.
	G <- .resolve_cohort_count_arg(G, R, "simulateDataCore")

	if (!is.null(seed)) {
		set.seed(seed)
	}

	# Normalize assignment_type defensively (backward-compat for callers
	# that pass NULL directly).
	if (is.null(assignment_type)) {
		assignment_type <- "marginal"
	}
	assignment_type <- match.arg(
		assignment_type,
		c("marginal", "multinomial", "ordered")
	)
	if (assignment_type != "marginal" && is.null(assignment_coefs)) {
		stop(sprintf(
			"assignment_type = '%s' requires non-NULL assignment_coefs",
			assignment_type
		))
	}

	res <- testGenRandomDataInputs(
		beta = beta,
		G = G,
		T = T,
		d = d,
		N = N,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq
	)

	num_treats <- res$num_treats
	p_expected <- res$p_expected

	rm(res)

	# Generate base effects and covariates. Marginal path keeps the
	# pre-1.14.0 byte-identical behavior; covariate path samples per-unit
	# cohort labels from the propensity-score model in `assignment_coefs`
	# and sorts by W to preserve the row-block invariant.
	if (assignment_type == "marginal") {
		res_base <- generateBaseEffects(
			N = N,
			d = d,
			T = T,
			G = G,
			distribution = distribution,
			guarantee_rank_condition = guarantee_rank_condition
		)
	} else {
		res_base <- .generateBaseEffectsCovariate(
			N = N,
			d = d,
			T = T,
			G = G,
			distribution = distribution,
			guarantee_rank_condition = guarantee_rank_condition,
			assignment_coefs = assignment_coefs
		)
	}
	cohort_fe <- res_base$cohort_fe
	time_fe <- res_base$time_fe
	X_long <- res_base$X_long
	assignments <- res_base$assignments
	cohort_inds <- res_base$cohort_inds

	stopifnot(ncol(cohort_fe) == G)
	stopifnot(is.matrix(X_long))

	# Base matrix: cohort FE, time FE, and covariates (if any)
	X_base <- if (d > 0) {
		cbind(cohort_fe, time_fe, X_long)
	} else {
		cbind(cohort_fe, time_fe)
	}

	# Auxiliary cohort counts (indep_counts). Marginal path uses the
	# original uniform multinomial draw. Non-marginal path draws an
	# independent X' + W' under the SAME DGP so indep_counts is drawn
	# from the same distribution as the main panel (BL2 fix per
	# round-1 review). Same rank-condition retry guard applies.
	if (assignment_type == "marginal") {
		indep_assignments <- genAssignments(
			N = N,
			G = G,
			guarantee_rank_condition = guarantee_rank_condition,
			d = d
		)
	} else {
		indep_assignments <- .draw_indep_assignments_covariate(
			N = N,
			d = d,
			T = T,
			G = G,
			distribution = distribution,
			guarantee_rank_condition = guarantee_rank_condition,
			assignment_coefs = assignment_coefs
		)
	}

	if (d > 0) {
		res_ints <- generateFEInts(X_long, cohort_fe, time_fe, N, T, G, d)
		X_ints1 <- cbind(X_base, res_ints$X_long_cohort, res_ints$X_long_time)
	} else {
		X_ints1 <- X_base
	}

	first_inds_test <- getFirstInds(G = G, T = T)
	res_treat <- genTreatVarsSim(
		num_treats,
		N,
		T,
		G,
		assignments,
		cohort_inds,
		N_UNTREATED = assignments[1],
		first_inds_test = first_inds_test
	)
	treat_mat_long <- res_treat$treat_mat_long
	first_inds <- res_treat$first_inds

	X_ints2 <- cbind(X_ints1, treat_mat_long)

	if (d > 0) {
		stopifnot(ncol(cohort_fe) == G)
		X_long_treat <- genTreatInts(
			treat_mat_long = treat_mat_long,
			X_long = X_long,
			n_treats = num_treats,
			cohort_fe = cohort_fe,
			N = N,
			T = T,
			G = G,
			d = d,
			N_UNTREATED = assignments[1]
		)

		X_final <- cbind(X_ints2, X_long_treat)
	} else {
		X_final <- X_ints2
	}

	if (ncol(X_final) != p_expected) {
		stop(
			"Constructed design matrix with interactions has incorrect number of columns."
		)
	}

	unit_res <- rnorm(N, mean = 0, sd = sqrt(sig_eps_c_sq))
	y <- X_final %*%
		beta +
		rep(unit_res, each = T) +
		rnorm(N * T, mean = 0, sd = sqrt(sig_eps_sq))
	y <- y - mean(y)

	if (gen_ints) {
		X_ret <- X_final
	} else {
		# Return X with no interactions
		# Expected number of columns: p = G + (T - 1) + d + num_treats
		p_expected <- G + (T - 1) + d + num_treats

		X_ret <- cbind(X_base, treat_mat_long)
		if (ncol(X_ret) != p_expected) {
			stop(
				"Constructed design matrix without interactions has incorrect number of columns."
			)
		}
	}

	# Prepare dataframe for `fetwfe()`

	# We know that when gen_ints = FALSE, the design matrix X is:
	# X = [cohort_fe, time_fe, X_long, treat_mat_long]
	# The base part (cohort_fe, time_fe, X_long) has (G + (T-1) + d) columns.
	base_cols <- G + (T - 1) + d

	# The treatment dummy block is in columns (base_cols + 1) : (base_cols + num_treats)
	treat_dummy <- cbind(X_base, treat_mat_long)[,
		(base_cols + 1):(base_cols + num_treats),
		drop = FALSE
	]
	# For each row, the observed treatment indicator is 1 if any entry in treat_dummy is 1.
	treatment_vec <- as.integer(apply(treat_dummy, 1, function(x) {
		any(x == 1)
	}))

	# Extract the covariate columns from the base part: they are the last d columns of the base part.
	if (d > 0) {
		cov_cols <- seq(from = (G + (T - 1) + 1), to = (G + (T - 1) + d))
		covariates <- cbind(X_base, treat_mat_long)[, cov_cols, drop = FALSE]
	}

	# Construct a data frame with the panel structure.
	df_panel <- data.frame(
		time = rep(1:T, times = N),
		unit = rep(sprintf("unit%02d", 1:N), each = T),
		treatment = treatment_vec,
		y = as.numeric(y)
	)

	cov_names <- c()

	if (d > 0) {
		# Add covariate columns with names "cov1", "cov2", ...
		for (j in seq_len(d)) {
			cov_name_j <- paste0("cov", j)
			df_panel[[cov_name_j]] <- covariates[, j]
			cov_names <- c(cov_names, cov_name_j)
		}
	}

	stopifnot(length(cov_names) == d)

	# Ensure that time is integer and unit is character.
	df_panel$time <- as.integer(df_panel$time)
	df_panel$unit <- as.character(df_panel$unit)
	df_panel$treatment <- as.integer(df_panel$treatment)

	# confirm that outputs satisfy input requirements of fetwfe()
	# (Plug in default values for arguments not generated here)
	checkFetwfeInputs(
		pdata = df_panel,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = cov_names,
		indep_counts = indep_assignments,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		lambda.max = NA,
		lambda.min = NA,
		q = 0.5,
		verbose = FALSE,
		alpha = 0.05,
		add_ridge = FALSE
	)

	ret <- list(
		pdata = df_panel,
		X = X_ret,
		y = y,
		covs = cov_names,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		coefs = beta,
		first_inds = first_inds,
		N_UNTREATED = assignments[1],
		assignments = assignments,
		indep_counts = indep_assignments,
		p = p_expected,
		N = N,
		T = T,
		G = G,
		R = G,
		d = d,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq
	)

	class(ret) <- "FETWFE_simulated"

	return(ret)
}
