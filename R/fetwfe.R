#' @import glmnet
#' @importFrom stats rbinom rmultinom rnorm runif sd lm coef pnorm

#' @title Fused extended two-way fixed effects
#'
#' @description Implementation of fused extended two-way fixed effects.
#' Estimates overall ATT as well as CATT (cohort average treatment effects on
#' the treated units).
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
#' @param covs (Optional.) Character; a vector containing the names of the
#' columns for covariates. All of these columns are expected to contain integer,
#' numeric, or factor values, and any categorical values will be automatically
#' encoded as binary indicators. If no covariates are provided, the treatment
#' effect estimation will proceed, but it will only be valid under unconditional
#' versions of the parallel trends and no anticipation assumptions. Default is c().
#' @param indep_counts (Optional.) Integer; a vector. If you have a sufficiently
#' large number of units, you can optionally randomly split your data set in
#' half (with `N` units in each data set). The data for half of the units should
#' go in the `pdata` argument provided above. For the other `N` units, simply
#' provide the counts for how many units appear in the untreated cohort plus
#' each of the other `R` cohorts in this argument `indep_counts`. The benefit
#' of doing this is that the standard error for the average treatment effect
#' will be (asymptotically) exact instead of conservative. The length of
#' `indep_counts` must equal 1 plus the number of treated cohorts in `pdata`.
#' All entries of `indep_counts` must be strictly positive (if you are concerned
#' that this might not work out, maybe your data set is on the small side and
#' it's best to just leave your full data set in `pdata`). The sum of all the
#' counts in `indep_counts` must match the total number of units in `pdata`.
#' Default is NA (in which case conservative standard errors will be calculated
#' if `q < 1`.)
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
#' @param lambda.max (Optional.) Numeric. A penalty parameter `lambda` will be
#' selected over a grid search by BIC in order to select a single model. The
#' largest `lambda` in the grid will be `lambda.max`. If no `lambda.max` is
#' provided, one will be selected automatically. When `q <= 1`, the model
#' will be sparse, and ideally all of the following are true at once: the
#' smallest model (the one corresponding to `lambda.max`) selects close to 0
#' features, the largest model (the one corresponding to `lambda.min`) selects
#' close to `p` features, `nlambda` is large enough so that models are
#' considered at every feasible model size, and `nlambda` is small enough so
#' that the computation doesn't become infeasible. You may
#' want to manually tweak `lambda.max`, `lambda.min`, and `nlambda` to try
#' to achieve these goals, particularly if the selected model size is very
#' close to the model corresponding to `lambda.max` or `lambda.min`, which could
#' indicate that the range of `lambda` values was too narrow or coarse. You can
#' use the function outputs `lambda.max_model_size`, `lambda.min_model_size`, and
#' `lambda_star_model_size` to try to assess this. Default is NA.
#' @param lambda.min (Optional.) Numeric. The smallest `lambda` penalty
#' parameter that will be considered. See the description of `lambda.max` for
#' details. Default is NA.
#' @param nlambda (Optional.) Integer. The total number of `lambda` penalty
#' parameters that will be considered. See the description of `lambda.max` for
#' details. Default is 100.
#' @param q (Optional.) Numeric; determines what `L_q` penalty is used for the
#' fusion regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is
#' possible to get standard errors and confidence intervals. `q` = 2 is ridge
#' regression. See Faletto (2025) for details. Default is 0.5.
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
#' @param se_type Character; one of `"default"` (the package's
#' Assumption-F1-based standard error from the paper) or `"cluster"`
#' (an *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' bridge-selected support; see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' `"cluster"` is only meaningful when `q < 1` (the bridge oracle property
#' is required); for `q >= 1` the SE will be `NA` regardless of `se_type`.
#' Default is `"default"`.
#' @return An object of class \code{fetwfe} containing the following elements:
#' \item{att_hat}{The estimated overall average treatment effect for a randomly selected treated unit.}
#' \item{att_se}{If `q < 1`, a standard error for the ATT. If `indep_counts` was provided, this standard error is asymptotically exact; if not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out fallback). See the package vignette section "Testing the zero-effect null" for interpretation guidance under selection consistency.}
#' \item{att_selected}{Logical scalar; `TRUE` if `att_hat` is not exactly zero (i.e., at least one cohort's bridge-penalized coefficient survived selection), `FALSE` otherwise. Under FETWFE Theorem 6.2 (restriction selection consistency), `att_selected = FALSE` is the asymptotic statement that the truth is zero. For ridge (`q = 2`) the bridge solver does not zero coefficients, so this will typically be `TRUE`.}
#' \item{catt_hats}{A named vector containing the estimated average treatment effects for each cohort.}
#' \item{catt_ses}{If `q < 1`, a named vector containing the (asymptotically exact, non-conservative) standard errors for the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each cohort conditional on being treated, which was used in calculating `att_hat`. If `indep_counts` was provided, `cohort_probs` was calculated from that; otherwise, it was calculated from the counts of units in each treated cohort in `pdata`.}
#' \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names (`cohort`), average treatment effects (`estimate`), standard errors (`se`), `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), per-cohort p-values (`p_value`), and a `selected` logical flag (`TRUE` when the bridge penalty left the cohort's CATT nonzero). For selected-out cohorts (`selected = FALSE`), `p_value` is `NA` --- the inferential content lives in `selected`. The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`) `stop()` with a migration message pointing to the new name. See `NEWS.md` for the rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the interactions between the treatment effects for each cohort at each time and the covariates.}
#' \item{sig_eps_sq}{Either the provided `sig_eps_sq` or the estimated one, if a value wasn't provided.}
#' \item{sig_eps_c_sq}{Either the provided `sig_eps_c_sq` or the estimated one, if a value wasn't provided.}
#' \item{lambda.max}{Either the provided `lambda.max` or the one that was used, if a value wasn't provided. (This is returned to help with getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding to `lambda.max` (for `q <= 1`, this will be the smallest model size). As mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one that was used, if a value wasn't provided.}
#' \item{lambda.min_model_size}{The size of the selected model corresponding to `lambda.min` (for `q <= 1`, this will be the largest model size). As mentioned above, for `q <= 1` ideally this value is close to `p`.}
#' \item{lambda_star}{The value of `lambda` chosen by BIC. If this value is close to `lambda.min` or `lambda.max`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If this value is close to `lambda.max_model_size` or `lambda.min_model_size`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{N}{The final number of units that were in the data set used for estimation (after any units may have been removed because they were treated in the first time period).}
#' \item{T}{The number of time periods in the final data set.}
#' \item{R}{The final number of treated cohorts that appear in the final data set.}
#' \item{d}{The final number of covariates that appear in the final data set (after any covariates may have been removed because they contained missing values or all contained the same value for every unit).}
#' \item{p}{The final number of columns in the full set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"` or `"cluster"`).}
#' \item{y_mean}{Numeric scalar; the mean of the original (pre-centering)
#'   response. Stored so downstream methods (`augment()`, `predict()`)
#'   can return fitted values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the name of the response
#'   column in the original `pdata`. Consumed by `augment.<class>()`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the
#'   `time_var` / `unit_var` / `treatment` arguments the user passed.
#'   Consumed by `augment.<class>()` when auto-aligning a user-supplied
#'   panel to the fitted design (e.g., dropping first-period-treated
#'   units the estimator removed internally, and sorting rows to match
#'   the design matrix's internal `(unit, time)` order).}
#' \item{covs}{Character vector; the original `covs` argument the user
#'   passed (before any factor expansion the estimator performed
#'   internally). Consumed by `augment.<class>()`.}
#' \item{internal}{A list containing internal outputs that are typically not needed for interpretation:
#'   \describe{
#'     \item{X_ints}{The design matrix created containing all interactions, time and cohort dummies, etc.}
#'     \item{y}{The vector of responses, containing `nrow(X_ints)` entries.}
#'     \item{X_final}{The design matrix after applying the change in coordinates to fit the model and also multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{y_final}{The final response after multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{theta_hat}{The vector of estimated coefficients in the transformed (fused) space, including the intercept as the first element.}
#'     \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#'   }
#' }
#'
#' The object has methods for \code{print()}, \code{summary()}, and \code{coef()}. By default, \code{print()} and \code{summary()} only show the essential outputs. To see internal details, use \code{print(x, show_internal = TRUE)} or \code{summary(x, show_internal = TRUE)}. The \code{coef()} method returns the vector of estimated coefficients (\code{beta_hat}).
#'
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
#' library(bacondecomp)
#'
#' data(castle)
#'
#' # Response: the log homicide rate. Treatment: `cdl` records the share of
#' # the year the castle-doctrine law was in effect, so `cdl > 0` gives the
#' # absorbing 0/1 treatment indicator `fetwfe()` requires.
#' castle$l_homicide <- log(castle$homicide)
#' castle$treated <- as.integer(castle$cdl > 0)
#'
#' # No `covs` here: castle's smallest adoption cohorts contain a single
#' # state, so the design is rank-deficient once any covariate is added.
#' res <- fetwfe(
#'     pdata = castle,
#'     time_var = "year",
#'     unit_var = "state",
#'     treatment = "treated",
#'     response = "l_homicide",
#'     verbose = TRUE)
#'
#' # Print results with internal details
#' print(res, max_cohorts = Inf)
#'
#' @export
fetwfe <- function(
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
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default"
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

	# Capture original user-supplied args so they can be stored on the output
	# for downstream methods (augment / predict) that need to re-prep `data`.
	# `covs` in particular gets reassigned later to its post-factor-expansion
	# form; we want the original on the output.
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
		lambda.max = lambda.max,
		lambda.min = lambda.min,
		q = q,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		estimator_type = "fetwfe"
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
	R <- prep$R
	indep_count_data_available <- prep$indep_count_data_available

	rm(prep)

	res <- fetwfe_core(
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
		lambda.max = lambda.max,
		lambda.min = lambda.min,
		nlambda = nlambda,
		q = q,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		se_type = se_type
	)

	att_branch <- .select_att_branch(
		res,
		indep_count_data_available = indep_count_data_available,
		q = q
	)
	att_hat <- att_branch$att_hat
	att_se <- att_branch$att_se
	cohort_probs <- att_branch$cohort_probs
	cohort_probs_overall <- att_branch$cohort_probs_overall

	att_p_value <- .compute_p_values(att_hat, att_se)
	att_selected <- att_hat != 0

	# Create the main output list with essential results
	out <- list(
		att_hat = att_hat,
		att_se = att_se,
		att_p_value = att_p_value,
		att_selected = att_selected,
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
		lambda.max = res$lambda.max,
		lambda.max_model_size = res$lambda.max_model_size,
		lambda.min = res$lambda.min,
		lambda.min_model_size = res$lambda.min_model_size,
		lambda_star = res$lambda_star,
		lambda_star_model_size = res$lambda_star_model_size,
		N = res$N,
		T = res$T,
		R = res$R,
		d = res$d,
		p = res$p,
		alpha = alpha,
		se_type = se_type,
		indep_counts_used = indep_count_data_available,
		y_mean = y_mean,
		response_col_name = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs_orig
	)

	# Add internal outputs in a separate list
	out$internal <- list(
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		theta_hat = res$theta_hat,
		calc_ses = res$calc_ses
	)

	# Validate constructed object's contracts (#85). Catches malformed
	# output at construction time rather than at downstream-method
	# confusion time.
	.validate_fetwfe(out)

	# Add the fetwfe class
	class(out) <- "fetwfe"

	return(out)
}


#' Run FETWFE on Simulated Data
#'
#' @description
#' This function runs the fused extended two-way fixed effects estimator (\code{fetwfe()}) on
#' simulated data. It is simply a wrapper for \code{fetwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{fetwfe()}. So the outputs match \code{fetwfe()}, and the needed inputs
#' match their counterparts in \code{fetwfe()}.
#'
#' @param simulated_obj An object of class \code{"FETWFE_simulated"} containing the simulated panel
#' data and design matrix.
#' @param lambda.max (Optional.) Numeric. A penalty parameter `lambda` will be
#' selected over a grid search by BIC in order to select a single model. The
#' largest `lambda` in the grid will be `lambda.max`. If no `lambda.max` is
#' provided, one will be selected automatically. For `lambda <= 1`, the model
#' will be sparse, and ideally all of the following are true at once: the
#' smallest model (the one corresponding to `lambda.max`) selects close to 0
#' features, the largest model (the one corresponding to `lambda.min`) selects
#' close to `p` features, `nlambda` is large enough so that models are
#' considered at every feasible model size, and `nlambda` is small enough so
#' that the computation doesn't become infeasible. You may
#' want to manually tweak `lambda.max`, `lambda.min`, and `nlambda` to try
#' to achieve these goals, particularly if the selected model size is very
#' close to the model corresponding to `lambda.max` or `lambda.min`, which could
#' indicate that the range of `lambda` values was too narrow. You can use the
#' function outputs `lambda.max_model_size`, `lambda.min_model_size`, and
#' `lambda_star_model_size` to try to assess this. Default is NA.
#' @param lambda.min (Optional.) Numeric. The smallest `lambda` penalty
#' parameter that will be considered. See the description of `lambda.max` for
#' details. Default is NA.
#' @param nlambda (Optional.) Integer. The total number of `lambda` penalty
#' parameters that will be considered. See the description of `lambda.max` for
#' details. Default is 100.
#' @param q (Optional.) Numeric; determines what `L_q` penalty is used for the
#' fusion regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is
#' possible to get standard errors and confidence intervals. `q` = 2 is ridge
#' regression. See Faletto (2025) for details. Default is 0.5.
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
#' @param se_type Character; one of `"default"` (the package's
#' Assumption-F1-based standard error from the paper) or `"cluster"`
#' (an *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' bridge-selected support; see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' `"cluster"` is only meaningful when `q < 1` (the bridge oracle property
#' is required); for `q >= 1` the SE will be `NA` regardless of `se_type`.
#' Default is `"default"`.
#' @return An object of class \code{fetwfe} containing the following elements:
#' \item{att_hat}{The estimated overall average treatment effect for a randomly selected treated unit.}
#' \item{att_se}{If `q < 1`, a standard error for the ATT. If `indep_counts` was provided, this standard error is asymptotically exact; if not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out fallback). See the package vignette section "Testing the zero-effect null" for interpretation guidance under selection consistency.}
#' \item{att_selected}{Logical scalar; `TRUE` if `att_hat` is not exactly zero (i.e., at least one cohort's bridge-penalized coefficient survived selection), `FALSE` otherwise. Under FETWFE Theorem 6.2 (restriction selection consistency), `att_selected = FALSE` is the asymptotic statement that the truth is zero. For ridge (`q = 2`) the bridge solver does not zero coefficients, so this will typically be `TRUE`.}
#' \item{catt_hats}{A named vector containing the estimated average treatment effects for each cohort.}
#' \item{catt_ses}{If `q < 1`, a named vector containing the (asymptotically exact, non-conservative) standard errors for the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each cohort conditional on being treated, which was used in calculating `att_hat`. If `indep_counts` was provided, `cohort_probs` was calculated from that; otherwise, it was calculated from the counts of units in each treated cohort in `pdata`.}
#' \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names (`cohort`), average treatment effects (`estimate`), standard errors (`se`), `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), per-cohort p-values (`p_value`), and a `selected` logical flag (`TRUE` when the bridge penalty left the cohort's CATT nonzero). For selected-out cohorts (`selected = FALSE`), `p_value` is `NA` --- the inferential content lives in `selected`. The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`) `stop()` with a migration message pointing to the new name. See `NEWS.md` for the rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the interactions between the treatment effects for each cohort at each time and the covariates.}
#' \item{sig_eps_sq}{Either the provided `sig_eps_sq` or the estimated one, if a value wasn't provided.}
#' \item{sig_eps_c_sq}{Either the provided `sig_eps_c_sq` or the estimated one, if a value wasn't provided.}
#' \item{lambda.max}{Either the provided `lambda.max` or the one that was used, if a value wasn't provided. (This is returned to help with getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding to `lambda.max` (for `q <= 1`, this will be the smallest model size). As mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one that was used, if a value wasn't provided.}
#' \item{lambda.min_model_size}{The size of the selected model corresponding to `lambda.min` (for `q <= 1`, this will be the largest model size). As mentioned above, for `q <= 1` ideally this value is close to `p`.}
#' \item{lambda_star}{The value of `lambda` chosen by BIC. If this value is close to `lambda.min` or `lambda.max`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If this value is close to `lambda.max_model_size` or `lambda.min_model_size`, that could suggest that the range of `lambda` values should be expanded.}
#' \item{N}{The final number of units that were in the data set used for estimation (after any units may have been removed because they were treated in the first time period).}
#' \item{T}{The number of time periods in the final data set.}
#' \item{R}{The final number of treated cohorts that appear in the final data set.}
#' \item{d}{The final number of covariates that appear in the final data set (after any covariates may have been removed because they contained missing values or all contained the same value for every unit).}
#' \item{p}{The final number of columns in the full set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"` or `"cluster"`).}
#' \item{y_mean}{Numeric scalar; the mean of the original (pre-centering)
#'   response. Stored so downstream methods (`augment()`, `predict()`)
#'   can return fitted values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the name of the response
#'   column in the original `pdata`. Consumed by `augment.<class>()`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the
#'   `time_var` / `unit_var` / `treatment` arguments the user passed.
#'   Consumed by `augment.<class>()` when auto-aligning a user-supplied
#'   panel to the fitted design (e.g., dropping first-period-treated
#'   units the estimator removed internally, and sorting rows to match
#'   the design matrix's internal `(unit, time)` order).}
#' \item{covs}{Character vector; the original `covs` argument the user
#'   passed (before any factor expansion the estimator performed
#'   internally). Consumed by `augment.<class>()`.}
#' \item{internal}{A list containing internal outputs that are typically not needed for interpretation:
#'   \describe{
#'     \item{X_ints}{The design matrix created containing all interactions, time and cohort dummies, etc.}
#'     \item{y}{The vector of responses, containing `nrow(X_ints)` entries.}
#'     \item{X_final}{The design matrix after applying the change in coordinates to fit the model and also multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{y_final}{The final response after multiplying on the left by the square root inverse of the estimated covariance matrix for each unit.}
#'     \item{theta_hat}{The vector of estimated coefficients in the transformed (fused) space, including the intercept as the first element.}
#'     \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#'   }
#' }
#'
#' The object has methods for \code{print()}, \code{summary()}, and \code{coef()}. By default, \code{print()} and \code{summary()} only show the essential outputs. To see internal details, use \code{print(x, show_internal = TRUE)} or \code{summary(x, show_internal = TRUE)}. The \code{coef()} method returns the vector of estimated coefficients (\code{beta_hat}).
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- fetwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
fetwfeWithSimulatedData <- function(
	simulated_obj,
	lambda.max = NA,
	lambda.min = NA,
	nlambda = 100,
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default"
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

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

	res <- fetwfe(
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
		nlambda = nlambda,
		q = q,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		se_type = se_type
	)

	return(res)
}


#' @title Extended two-way fixed effects
#'
#' @description Implementation of extended two-way fixed effects.
#' Estimates overall ATT as well as CATT (cohort average treatment effects on
#' the treated units).
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
#' @param covs (Optional.) Character; a vector containing the names of the
#' columns for covariates. All of these columns are expected to contain integer,
#' numeric, or factor values, and any categorical values will be automatically
#' encoded as binary indicators. If no covariates are provided, the treatment
#' effect estimation will proceed, but it will only be valid under unconditional
#' versions of the parallel trends and no anticipation assumptions. Default is c().
#' @param indep_counts (Optional.) Integer; a vector. If you have a sufficiently
#' large number of units, you can optionally randomly split your data set in
#' half (with `N` units in each data set). The data for half of the units should
#' go in the `pdata` argument provided above. For the other `N` units, simply
#' provide the counts for how many units appear in the untreated cohort plus
#' each of the other `R` cohorts in this argument `indep_counts`. The benefit
#' of doing this is that the standard error for the average treatment effect
#' will be (asymptotically) exact instead of conservative. The length of
#' `indep_counts` must equal 1 plus the number of treated cohorts in `pdata`.
#' All entries of `indep_counts` must be strictly positive (if you are concerned
#' that this might not work out, maybe your data set is on the small side and
#' it's best to just leave your full data set in `pdata`). The sum of all the
#' counts in `indep_counts` must match the total number of units in `pdata`.
#' Default is NA (in which case conservative standard errors will be calculated
#' if `q < 1`.)
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
#' @param se_type Character; one of `"default"` (the package's
#' Assumption-F1-based standard error from the paper) or `"cluster"`
#' (an *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support; see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' Default is `"default"`.
#' @return An object of class \code{etwfe} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{A standard error for the ATT. If the Gram matrix is not
#' invertible, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA`. Standard post-OLS interpretation; ETWFE does not
#' perform selection.}
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
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), and per-cohort
#' p-values (`p_value`). No `selected` column; ETWFE does not perform selection.
#' The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0
#' Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`,
#' `ConfIntHigh`, `P_value`) `stop()` with a migration message pointing to
#' the new name. See `NEWS.md` for the rename table.}
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
#' number of time periods in the final data set.} \item{R}{The final number of
#' treated cohorts that appear in the final data set.} \item{d}{The final number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"` or `"cluster"`).}
#' \item{y_mean}{Numeric scalar; the mean of the original (pre-centering)
#'   response. Stored so downstream methods (`augment()`, `predict()`)
#'   can return fitted values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the name of the response
#'   column in the original `pdata`. Consumed by `augment.<class>()`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the
#'   `time_var` / `unit_var` / `treatment` arguments the user passed.
#'   Consumed by `augment.<class>()` when auto-aligning a user-supplied
#'   panel to the fitted design.}
#' \item{covs}{Character vector; the original `covs` argument the user
#'   passed (before any factor expansion the estimator performed
#'   internally). Consumed by `augment.<class>()`.}
#' @author Gregory Faletto
#' @references
#' Wooldridge, J. M. (2021). Two-way fixed effects, the two-way mundlak
#' regression, and difference-in-differences estimators.
#' \emph{Available at SSRN 3906345}.
#' \doi{10.2139/ssrn.3906345}.
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
#' # No `covs` here: etwfe is pure OLS (no bridge penalty), and castle's
#' # smallest adoption cohorts contain a single state, so the design is
#' # rank-deficient once any covariate is added.
#' res <- etwfe(
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
etwfe <- function(
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
	se_type = "default"
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

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
	R <- prep$R
	indep_count_data_available <- prep$indep_count_data_available

	rm(prep)

	.check_cohort_rank_for_ols(
		in_sample_counts = in_sample_counts,
		R = R,
		d = d,
		add_ridge = add_ridge
	)

	res <- etwfe_core(
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

	# Create the main output list with essential results
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
		R = res$R,
		d = res$d,
		p = res$p,
		alpha = alpha,
		calc_ses = res$calc_ses,
		se_type = se_type,
		indep_counts_used = indep_count_data_available,
		y_mean = y_mean,
		response_col_name = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs = covs_orig
	)

	# Validate constructed object's contracts (#85).
	.validate_etwfe(out)

	# Add the etwfe class
	class(out) <- "etwfe"

	return(out)
}


#' Run ETWFE on Simulated Data
#'
#' @description
#' This function runs the extended two-way fixed effects estimator (\code{etwfe()}) on
#' simulated data. It is simply a wrapper for \code{etwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{etwfe()}. So the outputs match \code{etwfe()}, and the needed inputs
#' match their counterparts in \code{etwfe()}.
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
#' @param se_type Character; one of `"default"` (the package's
#' Assumption-F1-based standard error from the paper) or `"cluster"`
#' (an *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support; see the companion vignette `inference_vignette`
#' for the formula, the assumptions, and the theory-pending caveat).
#' Default is `"default"`.
#' @return An object of class \code{etwfe} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{A standard error for the ATT. If the Gram matrix is not
#' invertible, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA`. Standard post-OLS interpretation; ETWFE does not
#' perform selection.}
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
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), and per-cohort
#' p-values (`p_value`). No `selected` column; ETWFE does not perform selection.
#' The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0
#' Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`,
#' `ConfIntHigh`, `P_value`) `stop()` with a migration message pointing to
#' the new name. See `NEWS.md` for the rename table.}
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
#' number of time periods in the final data set.} \item{R}{The final number of
#' treated cohorts that appear in the final data set.} \item{d}{The final number
#' of covariates that appear in the final data set (after any covariates may
#' have been removed because they contained missing values or all contained the
#' same value for every unit).} \item{p}{The final number of columns in the full
#' set of covariates used to estimate the model.}
#' \item{alpha}{The alpha level used for confidence intervals.}
#' \item{calc_ses}{Logical indicating whether standard errors were calculated.}
#' \item{cohort_probs_overall}{A vector of the estimated cohort probabilities
#' on the overall sample (treated and untreated), used in computing the
#' variance of the overall ATT.}
#' \item{indep_counts_used}{Logical scalar; `TRUE` if a valid `indep_counts`
#' argument was provided and used for asymptotically-exact ATT inference,
#' `FALSE` otherwise.}
#' \item{se_type}{Character scalar; the `se_type` argument the user passed
#' (`"default"` or `"cluster"`).}
#' \item{y_mean}{Numeric scalar; the mean of the original (pre-centering)
#'   response. Stored so downstream methods (`augment()`, `predict()`) can
#'   return fitted values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the name of the response
#'   column in the original `pdata`. Consumed by `augment.<class>()`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the
#'   `time_var` / `unit_var` / `treatment` arguments the user passed.}
#' \item{covs}{Character vector; the original `covs` argument the user
#'   passed (before any factor expansion the estimator performed
#'   internally).}
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
#'
#'   result <- etwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
etwfeWithSimulatedData <- function(
	simulated_obj,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default"
) {
	se_type <- match.arg(se_type, c("default", "cluster"))

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

	res <- etwfe(
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
		se_type = se_type
	)

	return(res)
}
