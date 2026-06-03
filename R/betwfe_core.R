#' @title Bridge-penalized extended two-way fixed effects
#'
#' @description Implementation of extended two-way fixed effects with a bridge
#' penalty. Estimates overall ATT as well as CATT (cohort average treatment
#' effects on the treated units).
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
#' regularization. `q` = 1 is the lasso, and for 0 < `q` < 1, it is
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
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF); this is asymptotically exact for the package's default
#' cohort sample-proportions estimator and for every standard
#' propensity-score estimator that satisfies (Psi-IF) (multinomial logit,
#' any GLM on `W | X`, kernel/series regression of `1{W = g}` on `X`).
#' `"conservative"` returns the Cauchy-Schwarz upper bound from Theorem
#' (c); use only if the propensity-score estimator violates (Psi-IF).
#' `"cluster"` is an *experimental* unit-clustered Liang-Zeger sandwich
#' SE on the bridge-selected support (see the companion vignette
#' `inference_vignette` for details). `"cluster"` is only meaningful
#' when `q < 1` (the bridge oracle property is required); for `q >= 1`
#' the SE will be `NA` regardless of `se_type`. The default `"default"`
#' was the conservative Cauchy-Schwarz formula in versions <= 1.11.7;
#' v1.12.0 switched the default to the tight Gaussian variance.
#' Default is `"default"`.
#' @param lambda_selection Character; method for selecting the bridge
#'   penalty parameter `lambda`. Either `"cv"` (10-fold cross-validation
#'   on `cv.grpreg`; the v1.13.0+ default) or `"bic"` (BIC over the
#'   `grpreg` lambda grid; the prior default for v1.12.0 and earlier).
#'   The default changed in v1.13.0 to address a finite-sample bias issue
#'   documented in simulation studies (see issue #164). Pass
#'   `lambda_selection = "bic"` to recover the prior behavior. See the
#'   inference vignette section "Choosing the bridge penalty parameter"
#'   for details.
#' @param cv_folds Integer; number of folds for the CV path. Ignored when
#'   `lambda_selection = "bic"`. Default is 10.
#' @param cv_seed Integer or `NULL`; the seed passed to `set.seed()`
#'   immediately before the `cv.grpreg()` call. If `NULL` (the default),
#'   the seed defaults internally to `as.integer(N * T)`. Ignored when
#'   `lambda_selection = "bic"`.
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
#'   bounds and the per-cohort p-values (`p_value`) follow `ci_type`: under
#'   `"simultaneous"` the `p_value` is the single-step max-T multiplicity-
#'   adjusted p-value matching the band, under `"pointwise"` the per-cohort
#'   Wald p-value (#200). The standard errors (`se`) and selection flags
#'   (`selected`) are identical under both settings, and the overall-ATT
#'   confidence interval (a single scalar) is
#'   unaffected. When standard errors are unavailable (`q >= 1`, or a
#'   rank-deficient design) the bounds are `NA` under both settings. Default
#'   is `"simultaneous"`.
#' @return An object of class \code{betwfe} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{If `q < 1`, a standard error for the ATT. If
#' `indep_counts` was provided, this standard error is asymptotically exact; if
#' not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out
#' fallback).}
#' \item{att_selected}{Logical scalar; `TRUE` if `att_hat` is not exactly zero,
#' `FALSE` otherwise. BETWFE uses bridge regression directly on the
#' coefficients (rather than on the fused restrictions used by FETWFE); under
#' the bridge oracle property of Kock (2013), `att_selected = FALSE` is an
#' analogous asymptotic statement that the truth is zero under a sparsity
#' assumption different from the one Theorem 6.2 establishes for FETWFE. For
#' ridge (`q = 2`) the bridge solver does not zero coefficients, so this will
#' typically be `TRUE`.}
#' \item{catt_hats}{A named vector containing the estimated average treatment
#' effects for each cohort.} \item{catt_ses}{If `q < 1`, a named vector
#' containing the (asymptotically exact, non-conservative) standard errors for
#' the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names
#' (`cohort`), average treatment effects (`estimate`), standard errors (`se`),
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), per-cohort
#' p-values (`p_value`), and a `selected` logical flag (`TRUE` when the
#' bridge penalty left the cohort's CATT nonzero). For selected-out cohorts
#' (`selected = FALSE`), `p_value` is `NA`. The `catt_df` S3 class makes
#' `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names
#' (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`)
#' `stop()` with a migration message pointing to the new name. See
#' `NEWS.md` for the rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.} \item{lambda.max}{Either the provided `lambda.max` or the one
#' that was used, if a value wasn't provided. (This is returned to help with
#' getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding
#' `lambda.max` (for `q <= 1`, this will be the smallest model size). As
#' mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one
#' that was used, if a value wasn't provided.} \item{lambda.min_model_size}{The
#' size of the selected model corresponding to `lambda.min` (for `q <= 1`, this
#' will be the largest model size). As mentioned above, for `q <= 1` ideally
#' this value is close to `p`.}\item{lambda_star}{The value of `lambda` chosen
#' by the method recorded in `lambda_selection`. If this value is close to
#' `lambda.min` or `lambda.max`, that could suggest that the range of
#' `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If
#' this value is close to `lambda.max_model_size` or `lambda.min_model_size`,
#' That could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_selection}{Character scalar; either `"cv"` or `"bic"`.
#' Mirrors the `lambda_selection` argument the user passed.}
#' \item{cv_folds}{Integer scalar; the `cv_folds` value used when
#' `lambda_selection = "cv"`, `NA_integer_` when `lambda_selection = "bic"`.}
#' \item{cv_seed}{Integer scalar; the seed actually fed to `set.seed()`
#' immediately before `cv.grpreg()` was called. Defaults to
#' `as.integer(N * T)` when the user did not pass a seed. `NA_integer_`
#' when `lambda_selection = "bic"`.}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` / `eventStudy()` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
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
#' \item{y_mean}{Numeric scalar; mean of the original (pre-centering) response.
#' Stored so downstream methods (`augment()`, `predict()`) can return fitted
#' values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the response column name in
#' the original `pdata`. Consumed by `augment.betwfe()`.}
#' \item{time_var, unit_var, treatment}{Character scalars; the corresponding
#' arguments the user passed. Consumed by `augment.betwfe()` when auto-aligning
#' a user-supplied panel to the fitted design.}
#' \item{covs}{Character vector; the original `covs` argument (pre-factor-
#' expansion). Consumed by `augment.betwfe()`.}
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
#' library(bacondecomp)
#'
#' data(castle)
#'
#' # Response: the log homicide rate. Treatment: `cdl` records the share of
#' # the year the castle-doctrine law was in effect, so `cdl > 0` gives the
#' # absorbing 0/1 treatment indicator. No `covs`: castle's smallest
#' # adoption cohorts contain a single state, so the design is
#' # rank-deficient once any covariate is added.
#' castle$l_homicide <- log(castle$homicide)
#' castle$treated <- as.integer(castle$cdl > 0)
#'
#' # On this panel betwfe's bridge penalty selects every cohort out, so the
#' # estimated ATT and cohort effects below are all zero.
#' res <- betwfe(
#'     pdata = castle,
#'     time_var = "year",
#'     unit_var = "state",
#'     treatment = "treated",
#'     response = "l_homicide",
#'     verbose = TRUE)
#'
#' # Average treatment effect on the treated units (in percentage point
#' # units)
#' 100 * res$att_hat
#'
#' # Conservative 95% confidence interval for ATT (in percentage point units)
#'
#' low_att <- 100 * (res$att_hat - qnorm(1 - 0.05 / 2) * res$att_se)
#' high_att <- 100 * (res$att_hat + qnorm(1 - 0.05 / 2) * res$att_se)
#'
#' c(low_att, high_att)
#'
#' # Cohort average treatment effects and confidence intervals (in percentage
#' # point units)
#'
#' catt_df_pct <- res$catt_df
#' catt_df_pct[["estimate"]] <- 100 * catt_df_pct[["estimate"]]
#' catt_df_pct[["se"]] <- 100 * catt_df_pct[["se"]]
#' catt_df_pct[["ci_low"]] <- 100 * catt_df_pct[["ci_low"]]
#' catt_df_pct[["ci_high"]] <- 100 * catt_df_pct[["ci_high"]]
#'
#' catt_df_pct
#' @export
betwfe <- function(
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
	se_type = "default",
	lambda_selection = "cv",
	cv_folds = 10L,
	cv_seed = NULL,
	ci_type = c("simultaneous", "pointwise")
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	ci_type <- match.arg(ci_type)
	# `lambda_selection` validated downstream by `checkFetwfeInputs()`
	# (collect-all-violations pattern).

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
		lambda.max = lambda.max,
		lambda.min = lambda.min,
		q = q,
		verbose = verbose,
		alpha = alpha,
		add_ridge = add_ridge,
		allow_no_never_treated = allow_no_never_treated,
		estimator_type = "fetwfe",
		lambda_selection = lambda_selection,
		cv_folds = cv_folds,
		cv_seed = cv_seed
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
	G <- prep$G
	indep_count_data_available <- prep$indep_count_data_available

	rm(prep)

	res <- betwfe_core(
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
		se_type = se_type,
		lambda_selection = lambda_selection,
		cv_folds = cv_folds,
		cv_seed = cv_seed
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
		# v1.13.0 (#164): lambda-selection method provenance. Mirrors
		# fetwfe()'s shape. `res$cv_seed_used` is already NA_integer_
		# under the BIC path (the core's dispatch sets it on that branch).
		lambda_selection = lambda_selection,
		cv_folds = if (lambda_selection == "cv") {
			as.integer(cv_folds)
		} else {
			NA_integer_
		},
		cv_seed = res$cv_seed_used,
		X_ints = res$X_ints,
		y = res$y,
		X_final = res$X_final,
		y_final = res$y_final,
		N = res$N,
		T = res$T,
		G = res$G,
		R = res$G,
		d = res$d,
		p = res$p,
		calc_ses = res$calc_ses,
		alpha = alpha,
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
	# Validate constructed object's contracts (#85).
	.validate_betwfe(out)
	class(out) <- "betwfe"
	# Apply ci_type to the cohort-family bounds (#197). No-op unless
	# ci_type == "simultaneous"; runs after classing. Re-validates internally.
	out <- .finalize_ci_type(out, alpha = alpha)
	return(out)
}

#' Run BETWFE on Simulated Data
#'
#' @description
#' This function runs the bridge-penalized extended two-way fixed effects estimator (\code{betwfe()}) on
#' simulated data. It is simply a wrapper for \code{betwfe()}: it accepts an object of class
#' \code{"FETWFE_simulated"} (produced by \code{simulateData()}) and unpacks the necessary
#' components to pass to \code{betwfe()}. So the outputs match \code{betwfe()}, and the needed inputs
#' match their counterparts in \code{betwfe()}.
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
#' @param se_type Character; one of `"default"`, `"conservative"`, or
#' `"cluster"`. `"default"` returns the tight Gaussian variance
#' `sqrt(att_var_1 + att_var_2)` from Theorem (c$'$) under Assumption
#' (Psi-IF); this is asymptotically exact for the package's default
#' cohort sample-proportions estimator and for every standard
#' propensity-score estimator that satisfies (Psi-IF) (multinomial logit,
#' any GLM on `W | X`, kernel/series regression of `1{W = g}` on `X`).
#' `"conservative"` returns the Cauchy-Schwarz upper bound from Theorem
#' (c); use only if the propensity-score estimator violates (Psi-IF).
#' `"cluster"` is an *experimental* unit-clustered Liang-Zeger sandwich
#' SE on the bridge-selected support (see the companion vignette
#' `inference_vignette` for details). `"cluster"` is only meaningful
#' when `q < 1` (the bridge oracle property is required); for `q >= 1`
#' the SE will be `NA` regardless of `se_type`. The default `"default"`
#' was the conservative Cauchy-Schwarz formula in versions <= 1.11.7;
#' v1.12.0 switched the default to the tight Gaussian variance.
#' Default is `"default"`.
#' @param lambda_selection Character; method for selecting the bridge
#'   penalty parameter `lambda`. Either `"cv"` (10-fold cross-validation
#'   on `cv.grpreg`; the v1.13.0+ default) or `"bic"` (BIC over the
#'   `grpreg` lambda grid; the prior default for v1.12.0 and earlier).
#'   The default changed in v1.13.0 to address a finite-sample bias issue
#'   documented in simulation studies (see issue #164). Pass
#'   `lambda_selection = "bic"` to recover the prior behavior. See the
#'   inference vignette section "Choosing the bridge penalty parameter"
#'   for details.
#' @param cv_folds Integer; number of folds for the CV path. Ignored when
#'   `lambda_selection = "bic"`. Default is 10.
#' @param cv_seed Integer or `NULL`; the seed passed to `set.seed()`
#'   immediately before the `cv.grpreg()` call. If `NULL` (the default),
#'   the seed defaults internally to `as.integer(N * T)`. Ignored when
#'   `lambda_selection = "bic"`.
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
#'   bounds and the per-cohort p-values (`p_value`) follow `ci_type`: under
#'   `"simultaneous"` the `p_value` is the single-step max-T multiplicity-
#'   adjusted p-value matching the band, under `"pointwise"` the per-cohort
#'   Wald p-value (#200). The standard errors (`se`) and selection flags
#'   (`selected`) are identical under both settings, and the overall-ATT
#'   confidence interval (a single scalar) is
#'   unaffected. When standard errors are unavailable (`q >= 1`, or a
#'   rank-deficient design) the bounds are `NA` under both settings. Default
#'   is `"simultaneous"`.
#' @return An object of class \code{betwfe} containing the following elements:
#' \item{att_hat}{The
#' estimated overall average treatment effect for a randomly selected treated
#' unit.} \item{att_se}{If `q < 1`, a standard error for the ATT. If
#' `indep_counts` was provided, this standard error is asymptotically exact; if
#' not, it is asymptotically conservative. If `q >= 1`, this will be NA.}
#' \item{att_p_value}{A two-sided p-value for the overall ATT against the
#' null `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
#' `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out
#' fallback).}
#' \item{att_selected}{Logical scalar; `TRUE` if `att_hat` is not exactly zero,
#' `FALSE` otherwise. BETWFE uses bridge regression directly on the
#' coefficients (rather than on the fused restrictions used by FETWFE); under
#' the bridge oracle property of Kock (2013), `att_selected = FALSE` is an
#' analogous asymptotic statement that the truth is zero under a sparsity
#' assumption different from the one Theorem 6.2 establishes for FETWFE. For
#' ridge (`q = 2`) the bridge solver does not zero coefficients, so this will
#' typically be `TRUE`.}
#' \item{catt_hats}{A named vector containing the estimated average treatment
#' effects for each cohort.} \item{catt_ses}{If `q < 1`, a named vector
#' containing the (asymptotically exact, non-conservative) standard errors for
#' the estimated average treatment effects within each cohort.}
#' \item{cohort_probs}{A vector of the estimated probabilities of being in each
#' cohort conditional on being treated, which was used in calculating `att_hat`.
#' If `indep_counts` was provided, `cohort_probs` was calculated from that;
#' otherwise, it was calculated from the counts of units in each treated
#' cohort in `pdata`.} \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) displaying the cohort names
#' (`cohort`), average treatment effects (`estimate`), standard errors (`se`),
#' `1 - alpha` confidence interval bounds (`ci_low`, `ci_high`), per-cohort
#' p-values (`p_value`), and a `selected` logical flag (`TRUE` when the
#' bridge penalty left the cohort's CATT nonzero). For selected-out cohorts
#' (`selected = FALSE`), `p_value` is `NA`. The `catt_df` S3 class makes
#' `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names
#' (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`)
#' `stop()` with a migration message pointing to the new name. See
#' `NEWS.md` for the rename table.}
#' \item{beta_hat}{The full vector of estimated coefficients.}
#' \item{treat_inds}{The indices of `beta_hat` corresponding to
#' the treatment effects for each cohort at each time.}
#' \item{treat_int_inds}{The indices of `beta_hat` corresponding to the
#' interactions between the treatment effects for each cohort at each time and
#' the covariates.} \item{sig_eps_sq}{Either the provided `sig_eps_sq` or
#' the estimated one, if a value wasn't provided.} \item{sig_eps_c_sq}{Either
#' the provided `sig_eps_c_sq` or the estimated one, if a value wasn't
#' provided.} \item{lambda.max}{Either the provided `lambda.max` or the one
#' that was used, if a value wasn't provided. (This is returned to help with
#' getting a reasonable range of `lambda` values for grid search.)}
#' \item{lambda.max_model_size}{The size of the selected model corresponding
#' `lambda.max` (for `q <= 1`, this will be the smallest model size). As
#' mentioned above, for `q <= 1` ideally this value is close to 0.}
#' \item{lambda.min}{Either the provided `lambda.min` or the one
#' that was used, if a value wasn't provided.} \item{lambda.min_model_size}{The
#' size of the selected model corresponding to `lambda.min` (for `q <= 1`, this
#' will be the largest model size). As mentioned above, for `q <= 1` ideally
#' this value is close to `p`.}\item{lambda_star}{The value of `lambda` chosen
#' by the method recorded in `lambda_selection`. If this value is close to
#' `lambda.min` or `lambda.max`, that could suggest that the range of
#' `lambda` values should be expanded.}
#' \item{lambda_star_model_size}{The size of the model that was selected. If
#' this value is close to `lambda.max_model_size` or `lambda.min_model_size`,
#' That could suggest that the range of `lambda` values should be expanded.}
#' \item{lambda_selection}{Character scalar; either `"cv"` or `"bic"`.
#' Mirrors the `lambda_selection` argument the user passed.}
#' \item{cv_folds}{Integer scalar; the `cv_folds` value used when
#' `lambda_selection = "cv"`, `NA_integer_` when `lambda_selection = "bic"`.}
#' \item{cv_seed}{Integer scalar; the seed actually fed to `set.seed()`
#' immediately before `cv.grpreg()` was called. Defaults to
#' `as.integer(N * T)` when the user did not pass a seed. `NA_integer_`
#' when `lambda_selection = "bic"`.}
#' \item{ci_type}{Character scalar; the `ci_type` argument the user passed (`"simultaneous"` or `"pointwise"`), controlling whether the reported `catt_df` / `eventStudy()` confidence-interval bounds are simultaneous (family-wise) or pointwise.}
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
#' \item{y_mean}{Numeric scalar; mean of the original (pre-centering) response.
#' Stored so downstream methods (`augment()`, `predict()`) can return fitted
#' values on the original-response scale.}
#' \item{response_col_name}{Character scalar; the response column name in
#' the original `pdata`. Consumed by `augment.betwfe()`.}
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
#'   result <- betwfeWithSimulatedData(sim_data)
#' }
#'
#' @export
betwfeWithSimulatedData <- function(
	simulated_obj,
	lambda.max = NA,
	lambda.min = NA,
	nlambda = 100,
	q = 0.5,
	verbose = FALSE,
	alpha = 0.05,
	add_ridge = FALSE,
	allow_no_never_treated = TRUE,
	se_type = "default",
	lambda_selection = "cv",
	cv_folds = 10L,
	cv_seed = NULL,
	ci_type = c("simultaneous", "pointwise")
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	ci_type <- match.arg(ci_type)
	# `lambda_selection` validated downstream by `checkFetwfeInputs()`
	# (collect-all-violations pattern).

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

	res <- betwfe(
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
		se_type = se_type,
		lambda_selection = lambda_selection,
		cv_folds = cv_folds,
		cv_seed = cv_seed,
		ci_type = ci_type
	)

	return(res)
}


#' Core Estimation Logic for Bridge-Penalized Extended Two-Way Fixed Effects
#'
#' @description
#' This function implements the core estimation steps of the BETWFE methodology.
#' It takes a pre-processed design matrix and response, handles variance
#' components, performs bridge regression, selects the optimal penalty via BIC,
#' and calculates treatment effects and their standard errors.
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
#'   appear in the untreated cohort plus each of the other `G` cohorts, derived
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
#'   the coefficients to stabilize estimation. Default is `FALSE`.
#' @param se_type Character; the standard-error type, one of `"default"`
#'   (tight Gaussian variance under (Psi-IF), Theorem (c$'$)),
#'   `"conservative"` (Cauchy-Schwarz upper bound from Theorem (c) for
#'   non-(Psi-IF) propensity estimators), or `"cluster"` (experimental
#'   unit-clustered Liang-Zeger sandwich SE on the bridge-selected
#'   support). See the exported wrapper `betwfe()` for details. Default
#'   is `"default"`.
#'
#' @details
#' The function executes the following main steps:
#' \enumerate{
#'   \item **Input Checks:** Validates the provided parameters.
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
#'     penalty. The augmentation uses the identity basis (via
#'     `prep_for_etwfe_regression(..., is_fetwfe = FALSE)`) because BETWFE's
#'     design is untransformed --- same as ETWFE and `twfeCovs`, and unlike FETWFE
#'     which uses the inverse fusion-transform basis.
#'   \item **Cohort Probabilities:** Calculates cohort membership probabilities
#'     conditional on being treated, using `in_sample_counts` and `indep_counts`
#'     if available.
#'   \item **Bridge Regression:** Fits a bridge regression model using
#'     `grpreg::gBridge` on `X_final_scaled` and `y_final` with the specified `q`
#'     and lambda sequence.
#'   \item **Coefficient Selection (BIC):** Calls `getBetaBIC` to select the
#'     optimal `lambda` using BIC and retrieve the corresponding estimated
#'     coefficients.
#'   \item **Handle Zero-Feature Case:** If BIC selects a model with zero features,
#'     treatment effects are set to zero.
#'   \item **Treatment Effect Calculation:**
#'     \itemize{
#'       \item Extracts cohort-specific average treatment effects (CATTs) from
#'         `beta_hat`.
#'       \item Calls `getCohortATTsFinal` to calculate CATT point estimates,
#'         standard errors (if `q < 1`), and confidence intervals. This involves
#'         computing the Gram matrix and related quantities.
#'     }
#'   \item **Overall ATT Calculation:** Calls `getTeResultsOLS` to calculate the
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
#'   \item{catt_df}{A data frame (with S3 class `c("catt_df", "data.frame")`) summarizing CATTs (`cohort`, `estimate`, `se`, `ci_low`, `ci_high`, `p_value`) and a `selected` logical flag (`TRUE` when the bridge penalty left the cohort's CATT nonzero). The `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0 Title-Case column names `stop()` with a migration message pointing to the new name.}
#'   \item{beta_hat}{The vector of estimated coefficients.}
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
#'   \item{X_final}{The design matrix after GLS weighting.}
#'   \item{y_final}{The response vector after GLS weighting.}
#'   \item{N, T, G, d, p}{Dimensions used in estimation.}
#' @keywords internal
#' @noRd
betwfe_core <- function(
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
	add_ridge = FALSE,
	se_type = "default",
	lambda_selection = "cv",
	cv_folds = 10L,
	cv_seed = NULL
) {
	se_type <- match.arg(
		se_type,
		c("default", "conservative", "cluster")
	)
	lambda_selection <- match.arg(lambda_selection, c("cv", "bic"))
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

	G <- ret$G
	c_names <- ret$c_names
	indep_count_data_available <- ret$indep_count_data_available

	rm(ret)

	res <- prep_for_etwfe_regression(
		verbose = verbose,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		y = y,
		X_ints = X_ints,
		X_mod = X_ints, # Don't transform matrix
		N = N,
		T = T,
		G = G,
		d = d,
		p = p,
		num_treats = num_treats,
		add_ridge = add_ridge,
		first_inds = first_inds,
		in_sample_counts = in_sample_counts,
		indep_count_data_available = indep_count_data_available,
		indep_counts = indep_counts,
		is_fetwfe = FALSE
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

	# Dispatch on `lambda_selection` via the shared CV/BIC helper. BETWFE
	# uses the untransformed design (X_ints) for the BIC SSE computation,
	# unlike FETWFE which passes the fusion-transformed `X_mod`.
	bridge_sel <- .dispatch_bridge_selection(
		lambda_selection = lambda_selection,
		X_final_scaled = X_final_scaled,
		y_final = y_final,
		q = q,
		lambda.max = lambda.max,
		lambda.min = lambda.min,
		nlambda = nlambda,
		cv_folds = cv_folds,
		cv_seed = cv_seed,
		N = N,
		T = T,
		p = p,
		X_mod_bic = X_ints,
		y_bic = y,
		scale_center = scale_center,
		scale_scale = scale_scale,
		verbose = verbose
	)
	beta_hat <- bridge_sel$theta_hat
	lambda_star_ind <- bridge_sel$lambda_star_ind
	lambda_star_model_size <- bridge_sel$lambda_star_model_size
	fit <- bridge_sel$fit
	lambda.max <- bridge_sel$lambda.max
	lambda.min <- bridge_sel$lambda.min
	lambda.max_model_size <- bridge_sel$lambda.max_model_size
	lambda.min_model_size <- bridge_sel$lambda.min_model_size
	cv_seed_used <- bridge_sel$cv_seed_used

	lambda_star <- fit$lambda[lambda_star_ind]

	# c_names <- names(in_sample_counts)[2:(G + 1)] # Moved definition up
	stopifnot(length(c_names) == G)

	# Indices corresponding to base treatment effects
	ti <- .compute_treat_inds(
		G = G,
		T = T,
		d = d,
		num_treats = num_treats,
		p = p
	)
	treat_inds <- ti$treat_inds
	treat_int_inds <- ti$treat_int_inds

	# Handle edge case where no features are selected (model_size includes intercept)
	if (lambda_star_model_size <= 1 && all(beta_hat[2:(p + 1)] == 0)) {
		# Only the intercept might be non-zero. Delegate to the shared
		# helper (`.build_selected_out_result()` in `R/core_funcs.R`) that
		# also serves the no-treatment branch below and the two FETWFE
		# early-exits. BETWFE blocks omit `theta_hat` from the return.
		return(.build_selected_out_result(
			message_text = "No features selected (or only intercept); all treatment effects estimated to be 0.",
			verbose = verbose,
			G = G,
			c_names = c_names,
			q = q,
			beta_hat = beta_hat[2:(p + 1)], # Slopes are all zero
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			cohort_probs_overall = cohort_probs_overall,
			indep_cohort_probs_overall = indep_cohort_probs_overall,
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
			d = d,
			p = p,
			cv_seed_used = cv_seed_used
		))
	}

	# estimated coefficients
	beta_hat_slopes <- beta_hat[2:(p + 1)]

	# Indices of selected features
	sel_feat_inds <- which(beta_hat_slopes != 0)

	sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

	stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
	stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
	stopifnot(length(sel_treat_inds) <= length(treat_inds))
	stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

	# 1. Get beta_hat_slopes[treat_inds] -> these are the treatment coefficients
	# 2. Find which of these are non-zero: `which(beta_hat_slopes[treat_inds] != 0)` -> these are
	# indices *within* the treat_inds block. Let's call this `sel_treat_inds_relative_to_block`.
	# `sel_treat_inds` itself contains the global indices of selected treatment features.
	# So `beta_hat_slopes[sel_treat_inds]` are the non-zero treatment coefs.

	beta_hat_treat_block = beta_hat_slopes[treat_inds]
	sel_treat_inds_shifted <- which(beta_hat_treat_block != 0) # these are 1 to num_treats

	stopifnot(all(sel_treat_inds_shifted >= 1))
	stopifnot(all(sel_treat_inds_shifted <= num_treats))

	# Handle edge case where no treatment features selected
	if (length(sel_treat_inds_shifted) == 0) {
		# Apply ridge adjustment locally before early-exit return; doesn't
		# affect later code paths (they re-compute `beta_hat` from
		# `beta_hat_slopes` separately and run the non-early-exit
		# `add_ridge` scaling at line ~1108). Plan D3.
		if (add_ridge) {
			lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)
			beta_hat_slopes <- beta_hat_slopes * (1 + lambda_ridge)
		}

		return(.build_selected_out_result(
			message_text = "No treatment features selected; all treatment effects estimated to be 0.",
			verbose = verbose,
			G = G,
			c_names = c_names,
			q = q,
			beta_hat = beta_hat_slopes, # Untransformed slopes
			treat_inds = treat_inds,
			treat_int_inds = treat_int_inds,
			cohort_probs = cohort_probs,
			indep_cohort_probs = indep_cohort_probs,
			cohort_probs_overall = cohort_probs_overall,
			indep_cohort_probs_overall = indep_cohort_probs_overall,
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
			d = d,
			p = p,
			cv_seed_used = cv_seed_used
		))
	}

	beta_hat <- beta_hat_slopes

	# If using ridge regularization, multiply the "naive" estimated coefficients
	# by 1 + lambda_ridge, similar to suggestion in original elastic net paper.
	if (add_ridge) {
		lambda_ridge <- ifelse(is.na(lambda_ridge), 0, lambda_ridge)

		beta_hat <- beta_hat * (1 + lambda_ridge)
	}

	# Get actual estimated treatment effects
	tes <- beta_hat[treat_inds]

	stopifnot(length(tes) <= num_treats)

	stopifnot(all(beta_hat_slopes[treat_inds][sel_treat_inds_shifted] != 0))
	stopifnot(all(
		beta_hat_slopes[treat_inds][setdiff(
			1:num_treats,
			sel_treat_inds_shifted
		)] ==
			0
	))

	stopifnot(length(first_inds) == G)
	stopifnot(max(first_inds) <= num_treats)

	stopifnot(length(sel_feat_inds) > 0) # sel_feat_inds are indices in beta_hat_slopes
	stopifnot(length(sel_treat_inds_shifted) > 0) # sel_treat_inds_shifted are indices within the treat_inds block of beta_hat_slopes

	#
	#
	# Step 6: calculate cohort-specific treatment effects and standard
	# errors
	#
	#

	res <- getCohortATTsFinal(
		X_final = X_final, # This is X_mod * GLS_transform_matrix
		sel_feat_inds = sel_feat_inds, # Indices of non-zero elements in beta_hat_slopes
		treat_inds = treat_inds, # Global indices for treatment effects
		num_treats = num_treats,
		first_inds = first_inds,
		sel_treat_inds_shifted = sel_treat_inds_shifted, # Indices (1 to num_treats) of non-zero treat. coefs.
		c_names = c_names,
		tes = tes, # Treatment effect estimates (beta_hat[treat_inds])
		sig_eps_sq = sig_eps_sq,
		G = G,
		N = N,
		T = T,
		fused = FALSE,
		calc_ses = q < 1,
		include_selected = TRUE,
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

	#
	#
	# Step 7: calculate overall average treatment effect on treated units
	#
	#

	# Get overal estimated ATT!
	in_sample_te_results <- getTeResultsOLS(
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T,
		G = G,
		num_treats = num_treats,
		cohort_tes = cohort_tes, # CATTs (point estimates)
		cohort_probs = cohort_probs, # In-sample pi_g | treated
		psi_mat = psi_mat,
		gram_inv = gram_inv,
		tes = tes[sel_treat_inds_shifted], # Untransformed treatment effect estimates beta_hat[treat_inds]
		cohort_probs_overall = cohort_probs_overall, # In-sample pi_g (unconditional on treated)
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

	if ((q < 1) & calc_ses) {
		stopifnot(!is.na(in_sample_att_se))
	}

	if (indep_count_data_available) {
		stopifnot(nrow(psi_mat) == length(tes[sel_treat_inds_shifted]))

		indep_te_results <- getTeResultsOLS(
			sig_eps_sq = sig_eps_sq,
			N = N,
			T = T,
			G = G,
			num_treats = num_treats,
			cohort_tes = cohort_tes,
			cohort_probs = indep_cohort_probs, # indep pi_g | treated
			psi_mat = psi_mat,
			gram_inv = gram_inv,
			tes = tes[sel_treat_inds_shifted],
			cohort_probs_overall = indep_cohort_probs_overall, # indep pi_g (unconditional)
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
		beta_hat = beta_hat, # Untransformed slopes
		treat_inds = treat_inds,
		treat_int_inds = treat_int_inds,
		cohort_probs = cohort_probs,
		indep_cohort_probs = indep_cohort_probs,
		cohort_probs_overall = cohort_probs_overall,
		indep_cohort_probs_overall = indep_cohort_probs_overall,
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
		G = G,
		d = d,
		p = p,
		calc_ses = calc_ses,
		# v1.13.0 (#164): CV-path provenance.
		cv_seed_used = cv_seed_used
	))
}
