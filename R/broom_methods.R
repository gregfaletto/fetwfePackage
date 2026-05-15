#-------------------------------------------------------------------------------
# broom-package S3 methods for fetwfe / etwfe / betwfe / event_study / getTes.
#
# These methods plug fetwfe-family outputs into the broom ecosystem so users can
# pipe results into ggplot2 / modelsummary / gt without knowing the package's
# internal slot names. The generics come from the `generics` package (which
# `broom` re-exports); declaring them via @importFrom generics ... + plain
# @export on each method is the parsnip / estimatr idiom.
#-------------------------------------------------------------------------------

#' @importFrom generics tidy glance augment
NULL

#-------------------------------------------------------------------------------
# Internal helpers
#-------------------------------------------------------------------------------

#' @title Assemble a broom-style tidy data frame from an estimator output
#'
#' @description
#' Shared implementation behind `tidy.fetwfe()` / `tidy.etwfe()` /
#' `tidy.betwfe()`. Builds a data frame with row 1 = overall ATT and rows
#' 2..(R+1) = per-cohort CATTs, sorted by ascending cohort label. Columns
#' follow broom conventions: `term`, `estimate`, `std.error`, `statistic`,
#' `p.value`, and (when `conf.int = TRUE`) `conf.low` / `conf.high`. When
#' `include_selected = TRUE` (FETWFE / BETWFE), a logical `selected` column
#' is appended carrying the per-row selection status.
#'
#' @param x A fitted object of class `"fetwfe"`, `"etwfe"`, or `"betwfe"`.
#' @param conf.int Logical; include `conf.low` / `conf.high` columns.
#' @param conf.level Numeric in (0, 1); confidence level for the CI columns.
#' @param include_selected Logical; include the `selected` column.
#' @return A data frame with `R + 1` rows.
#' @keywords internal
#' @noRd
.tidy_estimator_output <- function(
	x,
	conf.int = TRUE,
	conf.level = 1 - x$alpha,
	include_selected = inherits(x, c("fetwfe", "betwfe"))
) {
	stopifnot(conf.level > 0, conf.level < 1)

	z <- stats::qnorm(1 - (1 - conf.level) / 2)

	att_term <- "ATT"
	att_estimate <- x$att_hat
	att_se <- x$att_se
	att_statistic <- if (is.na(att_se) || att_se == 0) {
		NA_real_
	} else {
		att_estimate / att_se
	}
	att_pvalue <- x$att_p_value

	catt <- x$catt_df
	catt <- catt[order(catt$Cohort), , drop = FALSE]
	cohort_terms <- paste0("Cohort ", catt$Cohort)
	cohort_estimates <- catt[["Estimated TE"]]
	cohort_ses <- catt$SE
	cohort_statistics <- ifelse(
		!is.na(cohort_ses) & cohort_ses > 0,
		cohort_estimates / cohort_ses,
		NA_real_
	)
	cohort_pvalues <- catt$P_value

	out <- data.frame(
		term = c(att_term, cohort_terms),
		estimate = c(att_estimate, cohort_estimates),
		std.error = c(att_se, cohort_ses),
		statistic = c(att_statistic, cohort_statistics),
		p.value = c(att_pvalue, cohort_pvalues),
		stringsAsFactors = FALSE
	)

	if (conf.int) {
		att_lo <- att_estimate - z * att_se
		att_hi <- att_estimate + z * att_se
		cohort_lo <- cohort_estimates - z * cohort_ses
		cohort_hi <- cohort_estimates + z * cohort_ses
		out$conf.low <- c(att_lo, cohort_lo)
		out$conf.high <- c(att_hi, cohort_hi)
	}

	if (include_selected) {
		att_selected <- x$att_selected
		cohort_selected <- if ("selected" %in% names(catt)) {
			catt$selected
		} else {
			rep(NA, nrow(catt))
		}
		out$selected <- c(att_selected, cohort_selected)
	}

	out
}

#' @title Glance helper for FETWFE / BETWFE (includes lambda_star columns)
#' @keywords internal
#' @noRd
.glance_fetwfe_betwfe <- function(x) {
	data.frame(
		nobs = x$N * x$T,
		n_units = x$N,
		n_periods = x$T,
		n_cohorts = x$R,
		n_covs = x$d,
		n_features = x$p,
		lambda_star = x$lambda_star,
		lambda_star_model_size = x$lambda_star_model_size,
		sig_eps_sq = x$sig_eps_sq,
		sig_eps_c_sq = x$sig_eps_c_sq,
		alpha = x$alpha,
		se_type = x$se_type,
		indep_counts_used = x$indep_counts_used,
		stringsAsFactors = FALSE
	)
}

#' @title Glance helper for ETWFE (omits lambda_star columns)
#' @keywords internal
#' @noRd
.glance_etwfe <- function(x) {
	data.frame(
		nobs = x$N * x$T,
		n_units = x$N,
		n_periods = x$T,
		n_cohorts = x$R,
		n_covs = x$d,
		n_features = x$p,
		sig_eps_sq = x$sig_eps_sq,
		sig_eps_c_sq = x$sig_eps_c_sq,
		alpha = x$alpha,
		se_type = x$se_type,
		indep_counts_used = x$indep_counts_used,
		stringsAsFactors = FALSE
	)
}

#' @title Class-dispatched access to the fitted design matrix
#'
#' @description
#' `fetwfe` wraps its design matrix in a nested `$internal` sublist
#' (`x$internal$X_ints`), while `etwfe` and `betwfe` put it at the top level
#' (`x$X_ints`). This helper hides that asymmetry from `augment()`.
#' @keywords internal
#' @noRd
.get_X_ints <- function(x) {
	if (inherits(x, "fetwfe")) {
		x$internal$X_ints
	} else {
		x$X_ints
	}
}

#' @title Shared augment implementation for fetwfe / etwfe / betwfe
#'
#' @description
#' Computes `.fitted = X %*% beta_hat + mean(data[[response]])` and
#' `.resid = data[[response]] - .fitted`, then column-binds to `data`. The
#' mean-add is needed because the estimator centers `y` during fitting
#' (`R/etwfe_core.R:1706`) and does not store the mean.
#' @keywords internal
#' @noRd
.augment_estimator_output <- function(x, data, response, ...) {
	if (missing(data)) {
		stop(
			"augment(): the fitted object does not store the original panel ",
			"data; supply `data = <your post-idCohorts() panel>` and ",
			"`response = <name of original response column>`."
		)
	}
	if (missing(response)) {
		stop(
			"augment(): supply `response = <name of original response ",
			"column>` so fitted values can be returned on the original-",
			"response scale (the estimator centers `y` internally during ",
			"fitting and does not store the response mean)."
		)
	}
	if (!is.character(response) || length(response) != 1L) {
		stop("augment(): `response` must be a single column name.")
	}
	if (!(response %in% names(data))) {
		stop(
			"augment(): column '",
			response,
			"' not found in `data`."
		)
	}

	X <- .get_X_ints(x)
	beta <- x$beta_hat

	if (nrow(data) != nrow(X)) {
		stop(
			"augment(): `data` has ",
			nrow(data),
			" rows but the fitted design has ",
			nrow(X),
			" (= x$N * x$T). `data` should match the panel the estimator ",
			"actually fit on; the estimator drops any units treated in the ",
			"first time period before fitting because their treatment effect ",
			"cannot be identified, so any such units must be removed from ",
			"`data` before calling augment()."
		)
	}

	y_obs <- data[[response]]
	if (!is.numeric(y_obs)) {
		stop("augment(): column '", response, "' must be numeric.")
	}
	y_mean <- mean(y_obs)

	fitted_centered <- as.vector(X %*% beta)
	.fitted <- fitted_centered + y_mean
	.resid <- y_obs - .fitted

	cbind(data, .fitted = .fitted, .resid = .resid)
}

#-------------------------------------------------------------------------------
# tidy() methods
#-------------------------------------------------------------------------------

#' Tidy an `fetwfe` fitted object
#'
#' Returns a `broom`-style tidy data frame for an object of class `"fetwfe"`.
#' Row 1 is the overall ATT (`term = "ATT"`); subsequent rows are the
#' cohort-specific ATTs (`term = "Cohort <adoption-time>"`), one per
#' treated cohort, sorted by ascending cohort label. Standard error,
#' z-statistic, and p-value reflect the value of `se_type` used at fit time
#' (model-based by default, cluster-robust under `se_type = "cluster"`).
#' Cohorts that the bridge penalty zeroed out (`selected = FALSE`) carry
#' `NA` for `std.error` / `statistic` / `p.value`.
#'
#' @param x An object of class `"fetwfe"` returned by [fetwfe()].
#' @param conf.int Logical. If `TRUE` (default), `conf.low` and `conf.high`
#'   columns are included.
#' @param conf.level Numeric in (0, 1). Confidence level for the CI columns.
#'   Defaults to `1 - x$alpha` (faithful to the alpha used at fit time;
#'   deviates from `broom::tidy.lm`'s `0.95` default by design).
#' @param ... Unused; present for S3 compatibility.
#' @return A data frame with `R + 1` rows and columns `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, optionally `conf.low` /
#'   `conf.high`, and `selected` (logical).
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(res)
#' }
#' @export
tidy.fetwfe <- function(
	x,
	conf.int = TRUE,
	conf.level = 1 - x$alpha,
	...
) {
	.tidy_estimator_output(x, conf.int, conf.level, include_selected = TRUE)
}

#' Tidy an `etwfe` fitted object
#'
#' Like [tidy.fetwfe()] but for an ETWFE fit. Has no `selected` column
#' (ETWFE does no regularized selection).
#'
#' @param x An object of class `"etwfe"` returned by [etwfe()].
#' @param conf.int Logical; include CI columns.
#' @param conf.level Numeric in (0, 1); defaults to `1 - x$alpha`.
#' @param ... Unused.
#' @return A data frame with `R + 1` rows.
#' @examples
#' \dontrun{
#'   res <- etwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(res)
#' }
#' @export
tidy.etwfe <- function(
	x,
	conf.int = TRUE,
	conf.level = 1 - x$alpha,
	...
) {
	.tidy_estimator_output(x, conf.int, conf.level, include_selected = FALSE)
}

#' Tidy a `betwfe` fitted object
#'
#' Like [tidy.fetwfe()] but for a BETWFE fit. Includes the `selected`
#' column reflecting BETWFE's bridge-penalized selection.
#'
#' @param x An object of class `"betwfe"` returned by [betwfe()].
#' @param conf.int Logical; include CI columns.
#' @param conf.level Numeric in (0, 1); defaults to `1 - x$alpha`.
#' @param ... Unused.
#' @return A data frame with `R + 1` rows.
#' @examples
#' \dontrun{
#'   res <- betwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(res)
#' }
#' @export
tidy.betwfe <- function(
	x,
	conf.int = TRUE,
	conf.level = 1 - x$alpha,
	...
) {
	.tidy_estimator_output(x, conf.int, conf.level, include_selected = TRUE)
}

#' Tidy a `fetwfe_event_study` object
#'
#' Returns a `broom`-style tidy data frame for the output of
#' [event_study()]. Renames the existing columns to broom conventions
#' (`se` → `std.error`, `ci_low` / `ci_high` → `conf.low` / `conf.high`,
#' `p_value` → `p.value`) and adds a `term` column (`"e<event_time>"`)
#' plus a `statistic` column (`estimate / std.error`) so the schema
#' matches `tidy.<estimator>()` for downstream `bind_rows()` consumers.
#'
#' @param x An object of class `"fetwfe_event_study"` returned by
#'   [event_study()].
#' @param ... Unused.
#' @return A data frame with one row per event-time and columns `term`,
#'   `event_time`, `n_cohorts`, `estimate`, `std.error`, `statistic`,
#'   `p.value`, `conf.low`, `conf.high`.
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(event_study(res))
#' }
#' @export
tidy.fetwfe_event_study <- function(x, ...) {
	statistic <- ifelse(
		!is.na(x$se) & x$se > 0,
		x$estimate / x$se,
		NA_real_
	)
	data.frame(
		term = paste0("e", x$event_time),
		event_time = x$event_time,
		n_cohorts = x$n_cohorts,
		estimate = x$estimate,
		std.error = x$se,
		statistic = statistic,
		p.value = x$p_value,
		conf.low = x$ci_low,
		conf.high = x$ci_high,
		stringsAsFactors = FALSE
	)
}

#' Tidy a `FETWFE_tes` simulation truth object
#'
#' Returns a `broom`-style tidy data frame for the population-truth
#' object returned by [getTes()]. Row 1 is the overall true ATT
#' (`term = "ATT_true"`); subsequent rows are the true cohort ATTs
#' (`term = "Cohort <integer>"`, matching `print.FETWFE_tes`'s
#' rendering, since `actual_cohort_tes` is an unnamed numeric vector).
#' Standard error / statistic / p-value / CI columns are all
#' `NA_real_` — there is no sampling distribution for a population
#' truth.
#'
#' @param x An object of class `"FETWFE_tes"` returned by [getTes()].
#' @param ... Unused.
#' @return A data frame with `R + 1` rows.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   broom::tidy(getTes(coefs))
#' }
#' @export
tidy.FETWFE_tes <- function(x, ...) {
	R <- length(x$actual_cohort_tes)
	terms <- c("ATT_true", paste0("Cohort ", seq_len(R)))
	estimates <- c(x$att_true, x$actual_cohort_tes)
	data.frame(
		term = terms,
		estimate = estimates,
		std.error = rep(NA_real_, R + 1L),
		statistic = rep(NA_real_, R + 1L),
		p.value = rep(NA_real_, R + 1L),
		conf.low = rep(NA_real_, R + 1L),
		conf.high = rep(NA_real_, R + 1L),
		stringsAsFactors = FALSE
	)
}

#-------------------------------------------------------------------------------
# glance() methods
#-------------------------------------------------------------------------------

#' Glance an `fetwfe` fitted object
#'
#' Returns a one-row `broom`-style summary data frame with model-level
#' scalars: panel-shape counts (`nobs`, `n_units`, `n_periods`,
#' `n_cohorts`, `n_covs`, `n_features`), bridge-regression tuning
#' (`lambda_star`, `lambda_star_model_size`), variance components
#' (`sig_eps_sq`, `sig_eps_c_sq`), and inference settings (`alpha`,
#' `se_type`, `indep_counts_used`).
#'
#' @param x An object of class `"fetwfe"`.
#' @param ... Unused.
#' @return A one-row data frame with 13 columns.
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::glance(res)
#' }
#' @export
glance.fetwfe <- function(x, ...) {
	.glance_fetwfe_betwfe(x)
}

#' Glance an `etwfe` fitted object
#'
#' Like [glance.fetwfe()] but omits the `lambda_star` /
#' `lambda_star_model_size` columns — ETWFE has no regularization.
#'
#' @param x An object of class `"etwfe"`.
#' @param ... Unused.
#' @return A one-row data frame with 11 columns.
#' @examples
#' \dontrun{
#'   res <- etwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::glance(res)
#' }
#' @export
glance.etwfe <- function(x, ...) {
	.glance_etwfe(x)
}

#' Glance a `betwfe` fitted object
#'
#' Same schema as [glance.fetwfe()] (BETWFE also has regularization).
#'
#' @param x An object of class `"betwfe"`.
#' @param ... Unused.
#' @return A one-row data frame with 13 columns.
#' @examples
#' \dontrun{
#'   res <- betwfeWithSimulatedData(
#'     simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::glance(res)
#' }
#' @export
glance.betwfe <- function(x, ...) {
	.glance_fetwfe_betwfe(x)
}

#-------------------------------------------------------------------------------
# augment() methods
#-------------------------------------------------------------------------------

#' Augment user-supplied data with fitted values and residuals from a fetwfe fit
#'
#' Computes `.fitted = X %*% beta_hat + mean(data[[response]])` and
#' `.resid = data[[response]] - .fitted`, then column-binds those two
#' columns onto `data`. The mean-add is necessary because the
#' estimator centers `y` during fitting (`R/etwfe_core.R:1706`) and does
#' not store the response mean.
#'
#' Both arguments are required: `data` because the package does not
#' store the original panel on the fitted object, and `response`
#' because the package does not store the response column name. `data`
#' must already have first-period-treated units dropped (the same
#' transformation `idCohorts()` applies internally during fitting);
#' the method errors with a clear message if the row counts disagree.
#'
#' @param x An object of class `"fetwfe"`.
#' @param data A balanced panel `data.frame` with one row per
#'   unit-period, already preprocessed by `idCohorts()` so the row
#'   count matches the fitted design.
#' @param response Character; name of the original (uncentered)
#'   response column in `data`.
#' @param ... Unused.
#' @return A copy of `data` with two extra numeric columns: `.fitted`
#'   and `.resid`.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- fetwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata, response = "y")
#' }
#' @export
augment.fetwfe <- function(x, data, response, ...) {
	.augment_estimator_output(x, data, response, ...)
}

#' Augment user-supplied data with fitted values and residuals from an etwfe fit
#'
#' Same shape as [augment.fetwfe()], dispatched on class `"etwfe"`.
#'
#' @param x An object of class `"etwfe"`.
#' @param data Post-`idCohorts()` panel.
#' @param response Character column name.
#' @param ... Unused.
#' @return `data` with `.fitted` and `.resid` columns appended.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- etwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata, response = "y")
#' }
#' @export
augment.etwfe <- function(x, data, response, ...) {
	.augment_estimator_output(x, data, response, ...)
}

#' Augment user-supplied data with fitted values and residuals from a betwfe fit
#'
#' Same shape as [augment.fetwfe()], dispatched on class `"betwfe"`.
#'
#' @param x An object of class `"betwfe"`.
#' @param data Post-`idCohorts()` panel.
#' @param response Character column name.
#' @param ... Unused.
#' @return `data` with `.fitted` and `.resid` columns appended.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(R = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- betwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata, response = "y")
#' }
#' @export
augment.betwfe <- function(x, data, response, ...) {
	.augment_estimator_output(x, data, response, ...)
}
