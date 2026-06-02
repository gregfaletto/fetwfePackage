#' @title Two-way fixed effects with covariates and separate treatment effects
#' for each cohort
#'
#' @description **WARNING: This function should NOT be used for estimation. It
#' is a biased estimator of treatment effects.** Implementation of two-way fixed
#' effects with covariates and separate treatment effects for each cohort.
#' Estimates overall ATT as well as CATT (cohort average treatment effects on
#' the treated units). It is implemented only for the sake of the simulation
#' studies in Faletto (2025). This estimator is only unbiased under the
#' assumptions that treatment effects are homogeneous across covariates and are
#' identical within cohorts across all times since treatment.
#'
#' @param pdata Dataframe; the panel data set. Each row should represent an
#' observation of a unit at a time. Should contain columns as described below.
#' @param time_var Character; the name of a single column containing a variable
#' for the time period. This column is expected to contain integer values (for
#' example, years). Recommended encodings for dates include format YYYY, YYYYMM,
#' or YYYYMMDD, whichever is appropriate for your data.
#' @param unit_var Character; the name of a single column containing a variable
#' for each unit. This column is expected to contain character values (i.e. the
#' "name" of each unit).
#' @param treatment Character; the name of a single column containing a variable
#' for the treatment dummy indicator. This column is expected to contain integer
#' values, and in particular, should equal 0 if the unit was untreated at that
#' time and 1 otherwise. Treatment should be an absorbing state; that is, if
#' unit `i` is treated at time `t`, then it must also be treated at all times
#' `t` + 1, ..., `T`. Any units treated in the first time period will be removed
#' automatically. Please make sure yourself that at least some units remain
#' untreated at the final time period ("never-treated units").
#' @param response Character; the name of a single column containing the
#' response for each unit at each time. The response must be an integer or
#' numeric value.
#' @param covs (Optional.) Either a character vector containing the names of
#' the columns for covariates (e.g., `covs = c("x1", "x2")`), or a one-sided
#' formula (e.g., `covs = ~ x1 + x2`) -- the formula form mirrors the
#' convention used by `did::att_gt(xformla = ...)`. Only additive bare
#' variable names are supported in the formula form; for derived variables,
#' compute them in the data frame first and pass via the character-vector
#' form. All of these columns are expected to contain integer, numeric, or
#' factor values, and any categorical values will be automatically encoded
#' as binary indicators. If no covariates are provided, the treatment effect
#' estimation will proceed, but it will only be valid under unconditional
#' versions of the parallel trends and no anticipation assumptions. Default
#' is c().
#' @param indep_counts (Optional.) Integer; a vector. If you have a sufficiently
#' large number of units, you can optionally randomly split your data set in
#' half (with `N` units in each data set). The data for half of the units should
#' go in the `pdata` argument provided above. For the other `N` units, simply
#' provide the counts for how many units appear in the untreated cohort plus
#' each of the other `G` cohorts in this argument `indep_counts`. The benefit
#' of doing this is that the standard error for the average treatment effect
#' will be (asymptotically) exact instead of conservative. The length of
#' `indep_counts` must equal 1 plus the number of treated cohorts in `pdata`.
#' All entries of `indep_counts` must be strictly positive (if you are concerned
#' that this might not work out, maybe your data set is on the small side and
#' it's best to just leave your full data set in `pdata`). The sum of all the
#' counts in `indep_counts` must match the total number of units in `pdata`.
#' Default is NA (in which case the conservative standard error formula will
#' be used).
#' @param sig_eps_sq (Optional.) Numeric; the variance of the row-level IID
#' noise assumed to apply to each observation. See Section 2 of Faletto (2025)
#' for details. It is best to provide this variance if it is known (for example,
#' if you are using simulated data). If this variance is unknown, this argument
#' can be omitted, and the variance will be estimated by
#' REML on the linear mixed-effects model `y ~ X + (1 | unit)` via
#' `lme4::lmer` (Bates et al. 2015; Patterson & Thompson 1971). Default is NA.
#' @param sig_eps_c_sq (Optional.) Numeric; the variance of the unit-level IID
#' noise (random effects) assumed to apply to each observation. See Section 2 of
#' Faletto (2025) for details. It is best to provide this variance if it is
#' known (for example, if you are using simulated data). If this variance is
#' unknown, this argument can be omitted, and the variance will be estimated
#' by REML via `lme4::lmer` on the
#' linear mixed-effects model `y ~ X + (1 | unit)` (Bates et al. 2015;
#' Patterson & Thompson 1971). Default is NA.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @param allow_no_never_treated (Optional.) Logical; if `TRUE` (default) and
#' the input panel contains no never-treated units, the panel is auto-truncated
#' by dropping time periods at and after the latest cohort's start time --- the
#' units in that latest cohort then serve as the never-treated comparison group
#' in the retained sub-panel --- with a warning naming the dropped periods. If
#' `FALSE`, the estimator stops with an error in this case (the package's
#' behavior prior to version 1.5.6). The argument has no effect when the input
#' already contains never-treated units. Default is `TRUE`.
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF) (asymptotically exact for the package's default cohort
#' sample-proportions estimator); `"conservative"` returns the
#' Cauchy-Schwarz upper bound from Theorem (c) (use only when the
#' propensity-score estimator violates (Psi-IF)); `"cluster"` is an
#' *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support (see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' Default is `"default"`. v1.12.0 introduced the tight Gaussian default;
#' versions <= 1.11.7 used the conservative Cauchy-Schwarz formula as
#' the default.
#' @param ci_type Character; one of `"simultaneous"` (default) or
#'   `"pointwise"`. Controls the confidence-interval bounds reported for the
#'   cohort-specific ATTs (in `catt_df`). `"simultaneous"` reports parametric
#'   simultaneous (family-wise, uniform) bands computed via [simultaneousCIs()]:
#'   the band covers all cohort effects jointly with probability `1 - alpha`,
#'   matching the default presentation of `did::aggte(cband = TRUE)`.
#'   `"pointwise"` reports per-effect Wald intervals (each covers its own
#'   effect with probability `1 - alpha`, no joint guarantee --- the behavior
#'   of versions <= 1.15.1). Only the interval bounds change; the standard
#'   errors (`se`) and per-cohort p-values (`p_value`) are identical under both
#'   settings. `twfeCovs` estimates a single pooled effect per cohort, so only
#'   the cohort family is affected (it has no event-study surface). When
#'   standard errors are unavailable (e.g., a rank-deficient design) the bounds
#'   are `NA` under both settings. Default is `"simultaneous"`.
#' @return An object of class \code{twfeCovs} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{A standard error for the ATT. If the Gram matrix is not
#' invertible, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA`. Standard post-OLS interpretation; `twfeCovs` does
#' not perform selection.}
#' \item{catt_hats}{A named vector containing the
#' estimated average treatment effects for each cohort.} \item{catt_ses}{A named
#' vector containing the (asymptotically exact) standard errors for
#' the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names
#' (`cohort`), average treatment effects (`estimate`), standard errors (`se`),
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), and
#' per-cohort p-values (`p_value`). No `selected` column; `twfeCovs` does
#' not perform selection. The `catt_df` S3 class makes `[[` / `$` / `[`
#' access on the pre-1.11.0 Title-Case column names (`Cohort`, `Estimated
#' TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`) `stop()` with a
#' migration message pointing to the new name. See `NEWS.md` for the
#' rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.} \item{X_ints}{The design matrix created containing all
#' interactions, time and cohort dummies, etc.} \item{y}{The vector of
#' responses, containing `nrow(X_ints)` entries.} \item{X_final}{The design
#' matrix after applying the change in coordinates to fit the model and also
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{y_final}{The final response after
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{N}{The final number of units that
#' were in the  data set used for estimation (after any units may have been
#' removed because they were treated in the first time period).} \item{T}{The
#' number of time periods in the final data set.} \item{G}{The final number of
#' treated cohorts that appear in the final data set.} \item{R}{Deprecated alias
#' for \code{G}, retained for backward compatibility; populated with the same
#' value. Use \code{G}. Will be removed in a future release.} \item{d}{The final
#' number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#' \item{y_mean}{Numeric scalar; mean of the original (pre-centering) response.
#' Stored so downstream methods (`augment()`, `predict()`) can return fitted
#' values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the response column name in
#' the original `pdata`. Reserved for future `augment()` / `predict()` methods.}
#' \item{time_var, unit_var, treatment}{Character scalars; the corresponding
#' arguments the user passed.}
#' \item{covs}{Character vector; the original `covs` argument (pre-factor-
#' expansion).}
#' \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"`, `"conservative"`, or `"cluster"`).}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
#' \item{internal}{A list containing internal outputs that are typically
#'   not needed for interpretation, packaged here for parity with
#'   `fetwfe()` so downstream consumers can use a single canonical
#'   access path across all four estimator classes (#144). The first
#'   five sub-slots (`X_ints`, `y`, `X_final`, `y_final`, `calc_ses`)
#'   are also duplicated at top level for backward compat;
#'   `variance_components` and `first_year` live only under
#'   `$internal`:
#'   \describe{
#'     \item{X_ints}{The design matrix containing all interactions,
#'       time and cohort dummies, etc. Same value as top-level `X_ints`.}
#'     \item{y}{The vector of responses. Same as top-level `y`.}
#'     \item{X_final}{The design matrix after the change-of-coordinates
#'       step. Same as top-level `X_final`.}
#'     \item{y_final}{The transformed response vector. Same as top-level
#'       `y_final`.}
#'     \item{calc_ses}{Logical indicating whether standard errors were
#'       calculated. Same as top-level `calc_ses`.}
#'     \item{variance_components}{A list exposing the two variance pieces
#'       (`att_var_1`, `att_var_2`) plus paper-notation counterparts
#'       (`V_1`, `V_2`) and unit-scaled variance estimators
#'       (`tilde_v_N`, `hat_v_N`, `tilde_v_N_C`, `tilde_v_N_C_pi_hat`,
#'       `tilde_v_N_C_pi_hat_cons`, `tilde_v_N_cons`). The Wald CI is
#'       `[hat_T_N +- qnorm(1-alpha/2) * sqrt(tilde_v_N / N)]` (paper Eq.
#'       `conf.int.form`). New in v1.12.0 (issue #141 + #146).}
#'     \item{first_year}{Integer or numeric scalar; the first (earliest)
#'       `time_var` value in the panel after `idCohorts()` processing.
#'       Consumed by `eventStudy()` to map `cohort_probs`' cohort labels
#'       (treatment-start years) to 1-based panel-time-index offsets when
#'       the labels are integer-coercible. New in v1.13.3 (issue #174).}
#'   }
#' }
#' @author Gregory Faletto
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting
#' Linear Mixed-Effects Models Using lme4. \emph{Journal of Statistical
#' Software}, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.
#'
#' Patterson, H. D., & Thompson, R. (1971). Recovery of inter-block
#' information when block sizes are unequal. \emph{Biometrika}, 58(3),
#' 545-554.
#'
#' Pinheiro, J. C., & Bates, D. M. (2000). \emph{Mixed-Effects Models in
#' S and S-PLUS}. Springer.
#' @examples
#' \dontrun{
#' library(bacondecomp)
#'
#' data(castle)
#'
#' # Response: the log homicide rate. Treatment: `cdl` records the share of
#' # the year the castle-doctrine law was in effect, so `cdl > 0` gives the
#' # absorbing 0/1 treatment indicator.
#' castle$l_homicide <- log(castle$homicide)
#' castle$treated <- as.integer(castle$cdl > 0)
#'
#' # No `covs` here: twfeCovs is pure OLS (no bridge penalty), and castle's
#' # smallest adoption cohorts contain a single state, so the design is
#' # rank-deficient once any covariate is added.
#' res <- twfeCovs(
#'     pdata = castle,
#'     time_var = "year",
#'     unit_var = "state",
#'     treatment = "treated",
#'     response = "l_homicide",
#'     verbose = TRUE)
#'
#' # Print results
#' print(res, max_cohorts = Inf)
#' }
#' @export
twfeCovs <- function(
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
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default",
	ci_type = c("simultaneous", "pointwise")
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	ci_type <- match.arg(ci_type)

	# Normalize `covs` to a character vector if a one-sided formula was
	# supplied (#28).
	covs <- .process_covs_input(covs)

	covs_orig <- covs

	# Steps 3-5: input validation + auto-truncation + design-matrix prep.
	prep <- .run_estimator_input_prep(
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
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		estimator_type = "etwfe"
	)

	pdata <- prep$pdata
	covs <- prep$covs
	X_ints <- prep$X_ints
	y <- prep$y
	y_mean <- prep$y_mean
	N <- prep$N
	T <- prep$T
	d <- prep$d
	p <- prep$p
	in_sample_counts <- prep$in_sample_counts
	num_treats <- prep$num_treats
	first_inds <- prep$first_inds
	first_year <- prep$first_year
	R <- prep$R
	indep_count_data_available <- prep$indep_count_data_available

	rm(prep)

	.check_cohort_rank_for_ols(
		in_sample_counts = in_sample_counts,
		R = R,
		d = d,
		add_ridge = add_ridge
	)

	res <- twfeCovs_core(
		X_ints = X_ints,
		y = y,
		in_sample_counts = in_sample_counts,
		N = N,
		T = T,
		d = d,
		p = p,
		num_treats = num_treats,
		first_inds = first_inds,
		indep_counts = indep_counts,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		se_type = se_type
	)

	att_branch <- .select_att_branch(
		res,
		indep_count_data_available = indep_count_data_available
	)
	att_hat <- att_branch$att_hat
	att_se <- att_branch$att_se
	cohort_probs <- att_branch$cohort_probs
	cohort_probs_overall <- att_branch$cohort_probs_overall

	att_p_value <- .compute_p_values(att_hat, att_se)

	variance_components <- .build_variance_components(
		att_var_1 = att_branch$att_var_1,
		att_var_2 = att_branch$att_var_2,
		N = res$N,
		T = res$T,
		se_type = se_type,
		indep_counts_used = indep_count_data_available
	)

	out <- list(
		att_hat = att_hat,
		att_se = att_se,
		att_p_value = att_p_value,
		catt_hats = res$catt_hats,
		catt_ses = res$catt_ses,
		cohort_probs = cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		catt_df = res$catt_df,
		beta_hat = res$beta_hat,
		treat_inds = res$treat_inds,
		treat_int_inds = res$treat_int_inds,
		sig_eps_sq = res$sig_eps_sq,
		sig_eps_c_sq = res$sig_eps_c_sq,
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		N = res$N,
		T = res$T,
		G = res$R,
		R = res$R,
		d = res$d,
		p = res$p,
		calc_ses = res$calc_ses,
		se_type = se_type,
		indep_counts_used = indep_count_data_available,
		y_mean = y_mean,
		response_col_name = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs_orig,
		ci_type = ci_type
	)
	# Add internal outputs in a separate list for parity with `fetwfe()` (#144).
	# The first five sub-slots (`X_ints`, `y`, `X_final`, `y_final`,
	# `calc_ses`) are also duplicated at top level for backward compat;
	# `variance_components` and `first_year` live only under `$internal`
	# (#179, #180).
	out$internal <- list(
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		calc_ses = res$calc_ses,
		variance_components = variance_components,
		# v1.13.3 (#174): see `fetwfe()` for rationale.
		first_year = first_year
	)
	# Validate constructed object's contracts (#85). Validator operates
	# on the list shape regardless of class, then class is assigned.
	.validate_twfeCovs(out)
	class(out) <- "twfeCovs"
	# Apply ci_type to the cohort-family bounds (#197). twfeCovs has no
	# `alpha` slot; its catt_df is built at alpha = 0.05 (the entry-point
	# default), so the finalizer must be passed alpha = 0.05 explicitly.
	# No-op unless ci_type == "simultaneous". Re-validates internally.
	out <- .finalize_ci_type(out, alpha = 0.05)
	return(out)
}


#' Run twfeCovs on Simulated Data
#'
#' @description
#' This function runs the bridge-penalized extended two-way fixed effects estimator (\code{twfeCovs()}) on
#' simulated data. It is simply a wrapper for \code{twfeCovs()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{twfeCovs()}. So the outputs match \code{twfeCovs()}, and the needed inputs
#' match their counterparts in \code{twfeCovs()}.
#'
#' @param simulated_obj An object of class \code{"FETWFE_simulated"} containing the simulated panel
#' data and design matrix.
#' @param verbose Logical; if TRUE, more details on the progress of the function will
#' be printed as the function executes. Default is FALSE.
#' @param alpha Numeric; function will calculate (1 - `alpha`) confidence intervals
#' for the cohort average treatment effects that will be returned in `catt_df`.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#' regularization to the (untransformed) coefficients to stabilize estimation.
#' Default is FALSE.
#' @param allow_no_never_treated (Optional.) Logical; if `TRUE` (default) and
#' the input panel contains no never-treated units, the panel is auto-truncated
#' by dropping time periods at and after the latest cohort's start time --- the
#' units in that latest cohort then serve as the never-treated comparison group
#' in the retained sub-panel --- with a warning naming the dropped periods. If
#' `FALSE`, the estimator stops with an error in this case (the package's
#' behavior prior to version 1.5.6). The argument has no effect when the input
#' already contains never-treated units. Default is `TRUE`.
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF) (asymptotically exact for the package's default cohort
#' sample-proportions estimator); `"conservative"` returns the
#' Cauchy-Schwarz upper bound from Theorem (c) (use only when the
#' propensity-score estimator violates (Psi-IF)); `"cluster"` is an
#' *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support (see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' Default is `"default"`. v1.12.0 introduced the tight Gaussian default;
#' versions <= 1.11.7 used the conservative Cauchy-Schwarz formula as
#' the default.
#' @param ci_type Character; one of `"simultaneous"` (default) or
#'   `"pointwise"`. Controls the confidence-interval bounds reported for the
#'   cohort-specific ATTs (in `catt_df`). `"simultaneous"` reports parametric
#'   simultaneous (family-wise, uniform) bands computed via [simultaneousCIs()]:
#'   the band covers all cohort effects jointly with probability `1 - alpha`,
#'   matching the default presentation of `did::aggte(cband = TRUE)`.
#'   `"pointwise"` reports per-effect Wald intervals (each covers its own
#'   effect with probability `1 - alpha`, no joint guarantee --- the behavior
#'   of versions <= 1.15.1). Only the interval bounds change; the standard
#'   errors (`se`) and per-cohort p-values (`p_value`) are identical under both
#'   settings. `twfeCovs` estimates a single pooled effect per cohort, so only
#'   the cohort family is affected (it has no event-study surface). When
#'   standard errors are unavailable (e.g., a rank-deficient design) the bounds
#'   are `NA` under both settings. Default is `"simultaneous"`.
#' @return An object of class \code{twfeCovs} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{A standard error for the ATT. If `indep_counts` was
#' provided, this standard error is asymptotically exact; otherwise, it is
#' asymptotically conservative. If the Gram matrix is not invertible, this
#' will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA`. Standard post-OLS interpretation; `twfeCovs` does
#' not perform selection.}
#' \item{catt_hats}{A named vector containing the estimated average treatment
#' effects for each cohort.} \item{catt_ses}{A named vector containing the
#' (asymptotically exact, non-conservative) standard errors for the estimated
#' average treatment effects within each cohort. If the Gram matrix is not
#' invertible, the entries are NA.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names
#' (`cohort`), average treatment effects (`estimate`), standard errors (`se`),
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), and
#' per-cohort p-values (`p_value`). No `selected` column; `twfeCovs` does
#' not perform selection. The `catt_df` S3 class makes `[[` / `$` / `[`
#' access on the pre-1.11.0 Title-Case column names (`Cohort`, `Estimated
#' TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`) `stop()` with a
#' migration message pointing to the new name. See `NEWS.md` for the
#' rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.}
#' \item{X_ints}{The design matrix created containing all
#' interactions, time and cohort dummies, etc.} \item{y}{The vector of
#' responses, containing `nrow(X_ints)` entries.} \item{X_final}{The design
#' matrix after applying the change in coordinates to fit the model and also
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{y_final}{The final response after
#' multiplying on the left by the square root inverse of the estimated
#' covariance matrix for each unit.} \item{N}{The final number of units that
#' were in the  data set used for estimation (after any units may have been
#' removed because they were treated in the first time period).} \item{T}{The
#' number of time periods in the final data set.} \item{G}{The final number of
#' treated cohorts that appear in the final data set.} \item{R}{Deprecated alias
#' for \code{G}, retained for backward compatibility; populated with the same
#' value. Use \code{G}. Will be removed in a future release.} \item{d}{The final
#' number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#' \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"`, `"conservative"`, or `"cluster"`).}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
#' \item{y_mean}{Numeric scalar; mean of the original (pre-centering) response.
#' Stored so downstream methods (`augment()`, `predict()`) can return fitted
#' values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the response column name in
#' the original `pdata`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the corresponding
#' arguments the user passed.}
#' \item{covs}{Character vector; the original `covs` argument (pre-factor-
#' expansion).}
#' \item{internal}{A list containing internal outputs that are typically
#'   not needed for interpretation, packaged here for parity with
#'   `fetwfe()` so downstream consumers can use a single canonical
#'   access path across all four estimator classes (#144). The first
#'   five sub-slots (`X_ints`, `y`, `X_final`, `y_final`, `calc_ses`)
#'   are also duplicated at top level for backward compat;
#'   `variance_components` and `first_year` live only under
#'   `$internal`:
#'   \describe{
#'     \item{X_ints}{The design matrix containing all interactions,
#'       time and cohort dummies, etc. Same value as top-level `X_ints`.}
#'     \item{y}{The vector of responses. Same as top-level `y`.}
#'     \item{X_final}{The design matrix after the change-of-coordinates
#'       step. Same as top-level `X_final`.}
#'     \item{y_final}{The transformed response vector. Same as top-level
#'       `y_final`.}
#'     \item{calc_ses}{Logical indicating whether standard errors were
#'       calculated. Same as top-level `calc_ses`.}
#'     \item{variance_components}{A list exposing the two variance pieces
#'       (`att_var_1`, `att_var_2`) plus paper-notation counterparts
#'       (`V_1`, `V_2`) and unit-scaled variance estimators
#'       (`tilde_v_N`, `hat_v_N`, `tilde_v_N_C`, `tilde_v_N_C_pi_hat`,
#'       `tilde_v_N_C_pi_hat_cons`, `tilde_v_N_cons`). The Wald CI is
#'       `[hat_T_N +- qnorm(1-alpha/2) * sqrt(tilde_v_N / N)]` (paper Eq.
#'       `conf.int.form`). New in v1.12.0 (issue #141 + #146).}
#'     \item{first_year}{Integer or numeric scalar; the first (earliest)
#'       `time_var` value in the panel after `idCohorts()` processing.
#'       Consumed by `eventStudy()` to map `cohort_probs`' cohort labels
#'       (treatment-start years) to 1-based panel-time-index offsets when
#'       the labels are integer-coercible. New in v1.13.3 (issue #174).}
#'   }
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- twfeCovsWithSimulatedData(sim_data)
#' }
#'
#' @export
twfeCovsWithSimulatedData <- function(
	simulated_obj,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default",
	ci_type = c("simultaneous", "pointwise")
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	ci_type <- match.arg(ci_type)

	if (!inherits(simulated_obj, "FETWFE_simulated")) {
		stop("simulated_obj must be an object of class 'FETWFE_simulated'")
	}

	pdata <- simulated_obj$pdata
	time_var <- simulated_obj$time_var
	unit_var <- simulated_obj$unit_var
	treatment <- simulated_obj$treatment
	response <- simulated_obj$response
	covs <- simulated_obj$covs
	sig_eps_sq <- simulated_obj$sig_eps_sq
	sig_eps_c_sq <- simulated_obj$sig_eps_c_sq
	indep_counts <- simulated_obj$indep_counts

	res <- twfeCovs(
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
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		se_type = se_type,
		ci_type = ci_type
	)

	return(res)
}


#' Core Estimation Logic for twfeCovs
#'
#' @description
#' This function implements the core estimation steps of the twfeCovs methodology.
#' It takes a pre-processed design matrix and response, handles variance components, performs
#' ordinary least squares regression, and calculates treatment effects and their standard errors.
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
#' @param verbose Logical; if `TRUE`, prints progress messages. Default is `FALSE`.
#' @param alpha Numeric; significance level for confidence intervals (e.g., 0.05 for
#'   95% CIs). Default is 0.05.
#' @param add_ridge (Optional.) Logical; if TRUE, adds a small amount of ridge
#'   regularization to the (untransformed) coefficients to stabilize estimation.
#'   Default is FALSE.
#' @param se_type Character; the standard-error type, one of "default"
#'   (tight Gaussian variance under (Psi-IF), Theorem (c')),
#'   "conservative" (Cauchy-Schwarz upper bound from Theorem (c) for
#'   non-(Psi-IF) propensity estimators), or "cluster" (experimental
#'   unit-clustered Liang-Zeger sandwich SE on the OLS-selected
#'   support). See the exported wrapper twfeCovs() for details. Default
#'   is "default".
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters via
#'     `check_etwfe_core_inputs()`.
#'   \item **Design preparation + GLS weighting:** Calls
#'     `prep_for_etwfe_regression(X_mod = X_ints, is_fetwfe = FALSE,
#'     is_twfe_covs = TRUE)`. This estimates the variance components via REML
#'     (`estOmegaSqrtInv()`) when `sig_eps_sq` or `sig_eps_c_sq` are `NA`, then
#'     GLS-weights the design and response by
#'     `sqrt(sig_eps_sq) * Omega_sqrt_inv` via a Kronecker product. The
#'     `is_twfe_covs = TRUE` flag collapses the per-period treatment columns
#'     into one column per cohort (so the design has `p_short = R + T - 1 + d
#'     + R` columns instead of FETWFE's `p`). No fusion transformation is
#'     applied. Also computes cohort membership probabilities from
#'     `in_sample_counts` and (if provided) `indep_counts`.
#'   \item **Optional ridge penalty:** If `add_ridge = TRUE`, `X_final_scaled`
#'     and `y_final` are augmented to add a small L2 penalty on the
#'     coefficients. The fitted coefficients are then multiplied by
#'     `(1 + lambda_ridge)` per the elastic-net adjustment.
#'   \item **OLS regression:** Fits `lm(y ~ . + 0, df)` on the GLS-weighted
#'     design. No bridge penalty, no `lambda` grid, no BIC selection step.
#'     The coefficient vector is in the original (untransformed) parameter
#'     space and is returned as `beta_hat`.
#'   \item **Cohort-specific treatment effects:** Calls
#'     `getCohortATTsFinal()` (with `fused = FALSE`, `include_selected =
#'     FALSE`) to compute per-cohort point estimates, standard errors (when
#'     the Gram matrix is invertible), and confidence intervals.
#'   \item **Overall ATT:** Calls `getTeResultsOLS()` for the overall ATT
#'     point estimate and SE, using both in-sample probabilities and
#'     independent probabilities when `indep_counts` is provided.
#' }
#' The standard errors for CATTs are asymptotically exact. For ATT, if
#' `indep_counts` are provided, the SE is asymptotically exact; otherwise, it
#' is asymptotically conservative.
#'
#' @return A list containing detailed estimation results:
#'   \item{in_sample_att_hat}{Estimated overall ATT using in-sample cohort probabilities.}
#'   \item{in_sample_att_se}{Standard error for `in_sample_att_hat`.}
#'   \item{in_sample_att_se_no_prob}{SE for `in_sample_att_hat` ignoring variability from estimating cohort probabilities.}
#'   \item{indep_att_hat}{Estimated overall ATT using `indep_counts` cohort probabilities (NA if `indep_counts` not provided).}
#'   \item{indep_att_se}{Standard error for `indep_att_hat` (NA if not applicable).}
#'   \item{catt_hats}{A named vector of estimated CATTs for each cohort.}
#'   \item{catt_ses}{A named vector of SEs for `catt_hats` (NA when the Gram matrix is not invertible).}
#'   \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) summarizing CATTs (`cohort`, `estimate`, `se`, `ci_low`, `ci_high`, `p_value`). The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names `stop()` with a migration message pointing to the new name.}
#'   \item{beta_hat}{The vector of estimated coefficients in the *original*
#'     space (no bridge fusion transformation; `twfeCovs` is pure OLS).}
#'   \item{treat_inds}{Indices in `beta_hat` corresponding to base treatment effects.}
#'   \item{treat_int_inds}{Indices in `beta_hat` corresponding to treatment-covariate interactions.}
#'   \item{cohort_probs}{Estimated cohort probabilities conditional on being treated, from `in_sample_counts`.}
#'   \item{indep_cohort_probs}{Estimated cohort probabilities from `indep_counts` (NA if not provided).}
#'   \item{sig_eps_sq}{The (possibly estimated) variance of observation-level noise.}
#'   \item{sig_eps_c_sq}{The (possibly estimated) variance of unit-level random effects.}
#'   \item{X_ints}{The original input design matrix from `prepXints`.}
#'   \item{y}{The original input centered response vector from `prepXints`.}
#'   \item{X_final}{The design matrix after GLS weighting (no fusion
#'     transformation for `twfeCovs_core`).}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, R, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
twfeCovs_core <- function(
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
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	se_type = "default"
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
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

	stopifnot(length(c_names) == R)

	res <- prep_for_etwfe_regression(
		verbose = verbose,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		y = y,
		X_ints = X_ints,
		X_mod = X_ints, # Don't transform matrix
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
		indep_counts = indep_counts,
		is_fetwfe = FALSE,
		is_twfe_covs = TRUE
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
	# OLS regression and fitted-coefficient extraction
	#
	#

	p_short <- R + T - 1 + d + R
	treat_inds_short <- (R + T - 1 + d + 1):p_short
	first_inds <- 1:R

	df <- data.frame(y = y_final, X_final_scaled)

	stopifnot(all(!is.na(df)))
	stopifnot("y" %in% colnames(df))
	stopifnot(ncol(df) == p_short + 1)

	# Response already centered; no intercept needed
	fit <- lm(y ~ . + 0, df)

	stopifnot(length(coef(fit)) == p_short)
	stopifnot(length(scale_scale) == p_short)

	beta_hat_slopes <- coef(fit) / scale_scale

	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
		beta_hat_slopes <- beta_hat_slopes * (1 + lambda_ridge)
	}

	stopifnot(length(beta_hat_slopes) == p_short)
	stopifnot(all(!is.na(beta_hat_slopes)))

	stopifnot(max(treat_inds_short) == p_short)

	treat_int_inds <- c()

	stopifnot(length(treat_inds_short) == R)

	# Get actual estimated treatment effects (in original, untransformed space)
	tes <- beta_hat_slopes[treat_inds_short]

	stopifnot(all(!is.na(tes)))

	stopifnot(length(tes) == R)

	stopifnot(length(first_inds) == R)
	stopifnot(max(first_inds) <= R)

	#
	#
	# Calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = NULL, # OLS path: no penalty selection occurred
		treat_inds = treat_inds_short, # Global indices for treatment effects
		num_treats = R,
		first_inds = first_inds,
		sel_treat_inds_shifted = seq_len(R),
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat_slopes[treat_inds])
		sig_eps_sq = sig_eps_sq,
		R = R,
		N = N,
		T = T,
		fused = FALSE,
		calc_ses = TRUE,
		include_selected = FALSE, # twfeCovs has no bridge selection
		alpha = alpha,
		se_type = se_type,
		y_final = y_final
	)

	cohort_te_df <- res$cohort_te_df
	cohort_tes <- res$cohort_tes
	cohort_te_ses <- res$cohort_te_ses
	psi_mat <- res$psi_mat
	gram_inv <- res$gram_inv
	calc_ses <- res$calc_ses
	sandwich_full <- res$sandwich_full
	treat_block_mask <- res$treat_block_mask

	rm(res)

	stopifnot(nrow(psi_mat) == R)
	stopifnot(ncol(psi_mat) == R)

	#
	#
	# Calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	stopifnot(length(tes) == R)
	stopifnot(nrow(psi_mat) == length(tes))

	in_sample_te_results <- getTeResultsOLS(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		R = R,
		num_treats = R,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_r | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		tes = tes, # Untransformed treatment effect estimates beta_hat[treat_inds]
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_r (unconditional on treated)
		first_inds = first_inds,
		calc_ses = calc_ses,
		indep_probs = FALSE,
		se_type = se_type,
		sandwich_full = sandwich_full,
		treat_block_mask = treat_block_mask
	)

	in_sample_att_hat <- in_sample_te_results$att_hat
	in_sample_att_se <- in_sample_te_results$att_te_se
	in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob
	in_sample_att_var_1 <- in_sample_te_results$att_var_1
	in_sample_att_var_2 <- in_sample_te_results$att_var_2

	if (indep_count_data_available) {
		indep_te_results <- getTeResultsOLS(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			R = R,
			num_treats = R,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_r | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			tes = tes,
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_r (unconditional)
			first_inds = first_inds,
			calc_ses = calc_ses,
			indep_probs = TRUE,
			se_type = se_type,
			sandwich_full = sandwich_full,
			treat_block_mask = treat_block_mask
		)
		indep_att_hat <- indep_te_results$att_hat
		indep_att_se <- indep_te_results$att_te_se
		indep_att_var_1 <- indep_te_results$att_var_1
		indep_att_var_2 <- indep_te_results$att_var_2
	} else {
		indep_att_hat <- NA
		indep_att_se <- NA
		indep_att_var_1 <- NA
		indep_att_var_2 <- NA
	}

	return(list(
		in_sample_att_hat = in_sample_att_hat,
		in_sample_att_se = in_sample_att_se,
		in_sample_att_se_no_prob = in_sample_att_se_no_prob,
		in_sample_att_var_1 = in_sample_att_var_1,
		in_sample_att_var_2 = in_sample_att_var_2,
		indep_att_hat = indep_att_hat,
		indep_att_se = indep_att_se,
		indep_att_var_1 = indep_att_var_1,
		indep_att_var_2 = indep_att_var_2,
		catt_hats = cohort_tes, # Already named if applicable from getCohortATTsFinal
		catt_ses = cohort_te_ses, # Already named if applicable
		catt_df = cohort_te_df,
		beta_hat = beta_hat_slopes, # Untransformed slopes
		treat_inds = treat_inds_short,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		X_ints = X_ints,
		y = y,
		X_final = X_final,
		y_final = y_final,
		N = N,
		T = T,
		R = R,
		d = d,
		p = p_short,
		calc_ses = calc_ses
	))
}
