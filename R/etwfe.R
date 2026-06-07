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
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF); this is asymptotically exact for the package's default
#' cohort sample-proportions estimator and for every standard
#' propensity-score estimator that satisfies (Psi-IF) (multinomial logit,
#' any GLM on `W | X`, kernel/series regression of `1{W = g}` on `X`).
#' `"conservative"` returns the Cauchy-Schwarz upper bound
#' `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
#' Theorem (c); use only if the propensity-score estimator violates
#' (Psi-IF) (e.g., a Robins-Rotnitzky-augmented doubly-robust estimator,
#' which the package does not currently implement). `"cluster"` is an
#' *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support (see the companion vignette
#' `inference_vignette` for the formula, the assumptions, and the
#' theory-pending caveat). The default value of `"default"` corresponds
#' to the new tight Gaussian default introduced in version 1.12.0;
#' previous versions used the conservative Cauchy-Schwarz formula as the
#' default. To recover the prior conservative default behavior, pass
#' `se_type = "conservative"`.
#' @param ci_type Character; one of `"simultaneous"` (default) or
#'   `"pointwise"`. Controls the confidence-interval bounds reported for the
#'   cohort-specific ATTs (in `catt_df`) and the event-study effects (from
#'   [eventStudy()], shown by `print` / `summary` / `plot`, and surfaced by
#'   [broom::tidy()] on the fitted object and on the [eventStudy()] /
#'   [cohortStudy()] outputs). `"simultaneous"` reports parametric simultaneous
#'   (family-wise, uniform) bands computed via [simultaneousCIs()]: each
#'   family's band covers all of its effects jointly with probability
#'   `1 - alpha`, matching the default presentation of
#'   `did::aggte(cband = TRUE)`. `"pointwise"` reports per-effect Wald intervals
#'   (each covers its own effect with probability `1 - alpha`, no joint
#'   guarantee --- the behavior of versions <= 1.15.1). Both the interval
#'   bounds and the per-cohort p-values (`p_value`) follow `ci_type` (max-T
#'   multiplicity-adjusted under `"simultaneous"`, per-cohort Wald under
#'   `"pointwise"`; #200); the standard errors (`se`) and the overall-ATT
#'   confidence interval (a single scalar) are identical under both settings.
#'   When standard errors are unavailable
#'   (e.g., a rank-deficient design) the bounds are `NA` under both settings.
#'   Default is `"simultaneous"`.
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
#' number of time periods in the final data set.} \item{G}{The final number of
#' treated cohorts that appear in the final data set.} \item{R}{Deprecated alias
#' for \code{G}, retained for backward compatibility; populated with the same
#' value. Use \code{G}. Will be removed in a future release.} \item{d}{The final
#' number
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
#' (`"default"`, `"conservative"`, or `"cluster"`).}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` / `eventStudy()` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
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

	.assemble_ols_estimator(
		prep = prep,
		core_fn = etwfe_core,
		validator_fn = .validate_etwfe,
		class_name = "etwfe",
		field_order = "etwfe",
		indep_counts = indep_counts,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		se_type = se_type,
		ci_type = ci_type,
		response = response,
		time_var = time_var,
		unit_var = unit_var,
		treatment = treatment,
		covs_orig = covs_orig
	)
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
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF); this is asymptotically exact for the package's default
#' cohort sample-proportions estimator and for every standard
#' propensity-score estimator that satisfies (Psi-IF) (multinomial logit,
#' any GLM on `W | X`, kernel/series regression of `1{W = g}` on `X`).
#' `"conservative"` returns the Cauchy-Schwarz upper bound
#' `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
#' Theorem (c); use only if the propensity-score estimator violates
#' (Psi-IF) (e.g., a Robins-Rotnitzky-augmented doubly-robust estimator,
#' which the package does not currently implement). `"cluster"` is an
#' *experimental* unit-clustered Liang-Zeger sandwich SE on the
#' OLS-selected support (see the companion vignette
#' `inference_vignette` for the formula, the assumptions, and the
#' theory-pending caveat). The default value of `"default"` corresponds
#' to the new tight Gaussian default introduced in version 1.12.0;
#' previous versions used the conservative Cauchy-Schwarz formula as the
#' default. To recover the prior conservative default behavior, pass
#' `se_type = "conservative"`.
#' @param ci_type Character; one of `"simultaneous"` (default) or
#'   `"pointwise"`. Controls the confidence-interval bounds reported for the
#'   cohort-specific ATTs (in `catt_df`) and the event-study effects (from
#'   [eventStudy()], shown by `print` / `summary` / `plot`, and surfaced by
#'   [broom::tidy()] on the fitted object and on the [eventStudy()] /
#'   [cohortStudy()] outputs). `"simultaneous"` reports parametric simultaneous
#'   (family-wise, uniform) bands computed via [simultaneousCIs()]: each
#'   family's band covers all of its effects jointly with probability
#'   `1 - alpha`, matching the default presentation of
#'   `did::aggte(cband = TRUE)`. `"pointwise"` reports per-effect Wald intervals
#'   (each covers its own effect with probability `1 - alpha`, no joint
#'   guarantee --- the behavior of versions <= 1.15.1). Both the interval
#'   bounds and the per-cohort p-values (`p_value`) follow `ci_type` (max-T
#'   multiplicity-adjusted under `"simultaneous"`, per-cohort Wald under
#'   `"pointwise"`; #200); the standard errors (`se`) and the overall-ATT
#'   confidence interval (a single scalar) are identical under both settings.
#'   When standard errors are unavailable
#'   (e.g., a rank-deficient design) the bounds are `NA` under both settings.
#'   Default is `"simultaneous"`.
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
#' number of time periods in the final data set.} \item{G}{The final number of
#' treated cohorts that appear in the final data set.} \item{R}{Deprecated alias
#' for \code{G}, retained for backward compatibility; populated with the same
#' value. Use \code{G}. Will be removed in a future release.} \item{d}{The final
#' number
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
#' (`"default"`, `"conservative"`, or `"cluster"`).}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` / `eventStudy()` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
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
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
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
	se_type = "default",
	ci_type = c("simultaneous", "pointwise")
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	ci_type <- match.arg(ci_type)

	sim <- .unpack_simulated_obj(simulated_obj)

	res <- etwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		indep_counts = sim$indep_counts,
		sig_eps_sq = sim$sig_eps_sq,
		sig_eps_c_sq = sim$sig_eps_c_sq,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		se_type = se_type,
		ci_type = ci_type
	)

	return(res)
}
