#-------------------------------------------------------------------------------
# broom-package S3 methods for fetwfe / etwfe / betwfe / eventStudy / getTes.
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
#' is appended carrying the per-row selection status. The cohort-row CI
#' columns pass through the fit's `catt_df` bounds (reflecting `ci_type`,
#' #197); the overall-ATT row's CI is the scalar Wald interval at `conf.level`.
#' The cohort-row `p.value` likewise follows `ci_type` (the single-step max-T
#' multiplicity-adjusted p under `"simultaneous"`, the Wald p under
#' `"pointwise"`; #200), while `statistic` stays the raw Wald z
#' (`estimate / std.error`) --- so under the simultaneous default `p.value`
#' is the adjusted dual of the band rather than `2 * pnorm(-|statistic|)`.
#'
#' @param x A fitted object of class `"fetwfe"`, `"etwfe"`, or `"betwfe"`.
#' @param conf.int Logical; include `conf.low` / `conf.high` columns.
#' @param conf.level Numeric in (0, 1); confidence level for the overall-ATT
#'   row's CI only (the cohort rows pass through the fit-time `catt_df`
#'   bounds).
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
	.check_for_tidy(x)
	stopifnot(conf.level > 0, conf.level < 1)

	z <- stats::qnorm(1 - (1 - conf.level) / 2)

	att_term <- "ATT"
	att_estimate <- x$att_hat
	att_se <- x$att_se
	att_statistic <- if (!is.na(att_se) && att_se > 0) {
		att_estimate / att_se
	} else {
		NA_real_
	}
	att_pvalue <- x$att_p_value

	catt <- x$catt_df
	# Composite sort key: numeric-cohort-time first, character tiebreak.
	# See R/class_helpers.R::.truncate_catt for the rationale.
	catt <- catt[
		order(
			suppressWarnings(as.numeric(catt$cohort)),
			catt$cohort
		),
		,
		drop = FALSE
	]
	cohort_terms <- paste0("Cohort ", catt$cohort)
	cohort_estimates <- catt[["estimate"]]
	cohort_ses <- catt$se
	cohort_statistics <- ifelse(
		!is.na(cohort_ses) & cohort_ses > 0,
		cohort_estimates / cohort_ses,
		NA_real_
	)
	cohort_pvalues <- catt$p_value

	out <- data.frame(
		term = c(att_term, cohort_terms),
		estimate = c(att_estimate, cohort_estimates),
		std.error = c(att_se, cohort_ses),
		statistic = c(att_statistic, cohort_statistics),
		p.value = c(att_pvalue, cohort_pvalues),
		stringsAsFactors = FALSE
	)

	if (conf.int) {
		# Overall-ATT row (K = 1) stays the scalar Wald interval at
		# `conf.level` (pointwise == simultaneous for a scalar).
		att_lo <- att_estimate - z * att_se
		att_hi <- att_estimate + z * att_se
		# Cohort rows pass through `catt_df`'s stored bounds (#197, Option B),
		# so they reflect the fit's `ci_type` (simultaneous by default). `catt`
		# is the SAME post-sort frame `cohort_estimates` / `cohort_ses` came
		# from, so `catt$ci_low` / `catt$ci_high` are already in `cohort_terms`
		# order. This is a pass-through SIMPLIFICATION, not a pointwise
		# recompute: `conf.level` no longer governs the cohort-row CI columns
		# (they carry the fit-time alpha + ci_type); see the roxygen note.
		cohort_lo <- catt$ci_low
		cohort_hi <- catt$ci_high
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
	.check_for_glance(x)
	data.frame(
		nobs = x$N * x$T,
		n_units = x$N,
		n_periods = x$T,
		n_cohorts = x$R,
		n_covs = x$d,
		n_features = x$p,
		lambda_star = x$lambda_star,
		lambda_star_model_size = x$lambda_star_model_size,
		# v1.13.0 (#164): lambda-selection method.
		lambda_selection = x$lambda_selection,
		cv_folds = x$cv_folds,
		cv_seed = x$cv_seed,
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
	.check_for_glance(x)
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
#' Computes `.fitted = X %*% beta_hat + x$y_mean` and
#' `.resid = data[[x$response_col_name]] - .fitted`, then column-binds to
#' the (possibly auto-trimmed) panel. If `nrow(data) != nrow(X_ints)`, the
#' method calls the package's internal cohort-identification routine on
#' `data` using the `time_var` / `unit_var` / `treatment` / `covs` slots
#' stashed at fit time, which drops units treated in the first time
#' period (their treatment effect is unidentifiable). The returned data
#' frame therefore corresponds to the panel the estimator actually fit
#' on — broom's `augment.lm` convention of returning the model frame
#' rather than the user's input.
#' @keywords internal
#' @noRd
.augment_estimator_output <- function(x, data, ...) {
	.check_for_augment(x)
	if (missing(data)) {
		stop(
			"augment(): the fitted object does not store the original panel ",
			"data; supply `data = <pdata you originally passed to the ",
			"estimator>`."
		)
	}

	response <- x$response_col_name
	if (is.null(response)) {
		stop(
			"augment(): the fitted object does not carry `response_col_name` ",
			"(expected on objects produced by fetwfe / etwfe / betwfe / ",
			"twfeCovs from package version 1.9.0 and later). This looks ",
			"like a fitted object from an earlier dev build."
		)
	}
	if (!(response %in% names(data))) {
		stop(
			"augment(): the response column '",
			response,
			"' (from `x$response_col_name`) is not in `data`. Rename or ",
			"restore that column before calling augment()."
		)
	}

	X <- .get_X_ints(x)

	# The four metadata slots are needed both for auto-trim (idCohorts) and
	# for sorting `data` into X_ints' row order. Validate them up front.
	needed_slots <- c("time_var", "unit_var", "treatment", "covs")
	missing_slots <- needed_slots[
		vapply(needed_slots, function(s) is.null(x[[s]]), logical(1))
	]
	if (length(missing_slots) > 0) {
		stop(
			"augment(): the fitted object is missing slots ",
			paste(missing_slots, collapse = ", "),
			" needed to align `data` with the fitted design. This looks ",
			"like a fitted object from an earlier dev build."
		)
	}

	# Auto-align row counts via idCohorts() when needed (the same drop the
	# estimator applied internally during fitting). Skip when the user has
	# already pre-trimmed, in which case `data` is used as-is.
	if (nrow(data) != nrow(X)) {
		ret <- idCohorts(
			df = data,
			time_var = x$time_var,
			unit_var = x$unit_var,
			treat_var = x$treatment
		)
		# idCohorts() filters first-period-treated units at the UNIT level
		# (R/utility.R:126) and drops the treat_var column. Augment needs the
		# unit-level filter but wants to preserve every column the user
		# supplied (treatment, plus any extras like state names or panel
		# IDs). So we use ret$df only to identify the surviving unit set,
		# and filter the user's original `data` directly. Invariant relied
		# on: idCohorts has no partial-time row drops. If that ever changes,
		# the unit-set filter here will under-trim and the nrow() check
		# below will fire.
		kept_units <- unique(ret$df[[x$unit_var]])
		data <- data[
			data[[x$unit_var]] %in% kept_units,
			,
			drop = FALSE
		]
		if (nrow(data) != nrow(X)) {
			stop(
				"augment(): even after auto-trimming first-period-treated ",
				"units via idCohorts(), `data` has ",
				nrow(data),
				" rows but the fitted design has ",
				nrow(X),
				". The panel structure does not match what the estimator was ",
				"fit on."
			)
		}
	}

	# Sort data into the same (unit_var, time_var) ascending order that
	# prepXints used when building X_ints (see `prepXints()` and `addDummies()`
	# in R/design_matrix.R). Without this, .fitted[i] would correspond to a
	# different (unit, time) than data[i, ] when the user's panel comes in a
	# different order from what the estimator saw internally.
	data <- data[
		order(data[[x$unit_var]], data[[x$time_var]], decreasing = FALSE),
		,
		drop = FALSE
	]

	y_obs <- data[[response]]
	if (!is.numeric(y_obs)) {
		stop("augment(): column '", response, "' must be numeric.")
	}

	beta <- x$beta_hat

	y_mean <- x$y_mean
	if (is.null(y_mean)) {
		stop(
			"augment(): the fitted object does not carry `y_mean`. This ",
			"looks like a fitted object from an earlier dev build."
		)
	}

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
#' The cohort-row `conf.low` / `conf.high` columns pass through the fit's
#' stored `catt_df` bounds, so they reflect the fit's `ci_type` (#197):
#' simultaneous (family-wise) by default, or pointwise when the fit used
#' `ci_type = "pointwise"`. They are NOT recomputed from `conf.level` (see the
#' `conf.level` note). The overall-ATT row (row 1) is a scalar, so its CI is
#' the pointwise Wald interval at `conf.level` (pointwise == simultaneous for a
#' single effect).
#'
#' @param x An object of class `"fetwfe"` returned by [fetwfe()].
#' @param conf.int Logical. If `TRUE` (default), `conf.low` and `conf.high`
#'   columns are included.
#' @param conf.level Numeric in (0, 1). Applies only to the overall-ATT row
#'   (row 1), whose CI is recomputed at this level; defaults to `1 - x$alpha`
#'   (faithful to the alpha used at fit time; deviates from `broom::tidy.lm`'s
#'   `0.95` default by design). The cohort rows pass through the fit-time
#'   `catt_df` bounds (reflecting the fit's `ci_type` and the fit-time alpha)
#'   and are NOT recomputed at `conf.level` (#197) --- mirroring
#'   [tidy.cohortStudy()]'s stored-bounds behavior.
#' @param ... Unused; present for S3 compatibility.
#' @return A data frame with `R + 1` rows and columns `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, optionally `conf.low` /
#'   `conf.high`, and `selected` (logical).
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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
#' @param conf.level Numeric in (0, 1); defaults to `1 - x$alpha`. Applies only
#'   to the overall-ATT row; the cohort rows pass through the fit-time
#'   `catt_df` bounds (reflecting the fit's `ci_type`) and are not recomputed
#'   at `conf.level` (#197). See [tidy.fetwfe()].
#' @param ... Unused.
#' @return A data frame with `R + 1` rows.
#' @examples
#' \dontrun{
#'   res <- etwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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
#' @param conf.level Numeric in (0, 1); defaults to `1 - x$alpha`. Applies only
#'   to the overall-ATT row; the cohort rows pass through the fit-time
#'   `catt_df` bounds (reflecting the fit's `ci_type`) and are not recomputed
#'   at `conf.level` (#197). See [tidy.fetwfe()].
#' @param ... Unused.
#' @return A data frame with `R + 1` rows.
#' @examples
#' \dontrun{
#'   res <- betwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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

#' Tidy an `eventStudy` object
#'
#' Returns a `broom`-style tidy data frame for the output of
#' [eventStudy()]. Renames existing columns to broom conventions
#' (`se` \eqn{\to} `std.error`, `p_value` \eqn{\to} `p.value`) and adds a `term`
#' column (`"e<event_time>"`) plus a `statistic` column
#' (`estimate / std.error`) so the schema matches `tidy.<estimator>()`
#' for downstream `bind_rows()` consumers.
#'
#' The `eventStudy()` output stores its confidence-interval bounds (`ci_low` /
#' `ci_high`), which reflect the fit's `ci_type` (#197): simultaneous
#' (family-wise, uniform) by default, or pointwise when the fit used
#' `ci_type = "pointwise"`. When `conf.int = TRUE` (the default), `conf.low` /
#' `conf.high` PASS THROUGH those stored bounds rather than recomputing from
#' `estimate +/- z * se` --- so the tidied event-study CIs agree with
#' `print` / `summary` / `plot` and with [simultaneousCIs()] under the default.
#' When `conf.int = FALSE`, the CI columns are omitted. (Degenerate event times
#' carry `NA` bounds under both `ci_type` settings.)
#'
#' @param x An object of class `"eventStudy"` returned by
#'   [eventStudy()].
#' @param conf.int Logical; include `conf.low` / `conf.high` columns.
#' @param conf.level Numeric in (0, 1). Retained for `broom`-convention parity
#'   (default `0.95`) but no longer recomputes the event-study CIs: as of
#'   version 1.16.0 (#197) the `conf.low` / `conf.high` columns pass through the
#'   `eventStudy` object's stored `ci_low` / `ci_high` (reflecting the fit's
#'   `ci_type`). To change the confidence level, recompute `eventStudy(fit,
#'   alpha = ...)` at the desired alpha first.
#' @param ... Unused.
#' @return A data frame with one row per event-time and columns `term`,
#'   `event_time`, `n_cohorts`, `estimate`, `std.error`, `statistic`,
#'   `p.value`, and (when `conf.int = TRUE`) `conf.low` / `conf.high`.
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(eventStudy(res))
#' }
#' @export
tidy.eventStudy <- function(
	x,
	conf.int = TRUE,
	conf.level = 0.95,
	...
) {
	stopifnot(conf.level > 0, conf.level < 1)
	# Guard against a user-mutated `eventStudy` frame that's missing a
	# required column: localizes the error message rather than the cryptic
	# "arguments imply differing number of rows" from `data.frame()` when
	# `x$<col>` returns NULL on an absent column. Parallels the guard added
	# to `tidy.cohortStudy()` in PR #150. The CI columns `ci_low` / `ci_high`
	# ARE required now (#197, Option B): the method passes them through
	# (reflecting the fit's `ci_type`) rather than recomputing from
	# `estimate +/- z * se`.
	required <- c(
		"event_time",
		"n_cohorts",
		"estimate",
		"se",
		"ci_low",
		"ci_high",
		"p_value"
	)
	missing_cols <- setdiff(required, names(x))
	if (length(missing_cols) > 0L) {
		stop(
			"tidy.eventStudy(): input is missing required columns: ",
			paste(missing_cols, collapse = ", "),
			". If you have mutated the `eventStudy` object, restore the ",
			"original data frame; otherwise please file an issue.",
			call. = FALSE
		)
	}
	statistic <- ifelse(
		!is.na(x$se) & x$se > 0,
		x$estimate / x$se,
		NA_real_
	)
	out <- data.frame(
		term = paste0("e", x$event_time),
		event_time = x$event_time,
		n_cohorts = x$n_cohorts,
		estimate = x$estimate,
		std.error = x$se,
		statistic = statistic,
		p.value = x$p_value,
		stringsAsFactors = FALSE
	)
	if (conf.int) {
		# Pass through the eventStudy object's stored bounds (#197, Option B),
		# which are simultaneous-or-pointwise per the fit's `ci_type`. This is
		# a SIMPLIFICATION of the prior `estimate +/- z * se` recompute;
		# `conf.level` no longer governs these columns (see the roxygen note).
		out$conf.low <- x$ci_low
		out$conf.high <- x$ci_high
	}
	out
}

#' Tidy a `cohortStudy` object
#'
#' Returns a `broom`-style tidy data frame for the output of
#' [cohortStudy()]. Renames the snake_case columns of `catt_df` to broom
#' conventions (`se` -> `std.error`, `p_value` -> `p.value`,
#' `ci_low` / `ci_high` -> `conf.low` / `conf.high`) and adds a `term`
#' column (`"cohort_<cohort label>"`) plus a `statistic` column
#' (`estimate / std.error`) so the schema matches `tidy.eventStudy()`
#' for downstream `bind_rows()` consumers. When the input carries a
#' `selected` column (`fetwfe` / `betwfe`), it is passed through as the
#' final column.
#'
#' Confidence intervals come from the cohort fit's stored bounds (which
#' encode the alpha passed at fit time); unlike [tidy.eventStudy()], this
#' method does not recompute the CIs at a custom `conf.level` because the
#' standard errors in `catt_df` are already paired with the fit-time
#' bounds (`ci_low` / `ci_high`), so re-emitting those is the
#' minimum-surprise behavior.
#'
#' @param x A `cohortStudy` object returned by [cohortStudy()].
#' @param ... Unused; present for S3 compatibility.
#' @return A data frame with one row per treated cohort and columns
#'   `term`, `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
#'   `conf.high`, and (if present in the input) `selected`.
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
#'                  N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   )
#'   broom::tidy(cohortStudy(res))
#' }
#' @export
tidy.cohortStudy <- function(x, ...) {
	# Guard against a user-mutated `catt_df` slot that's missing required
	# columns: localizes the error to a clear message rather than the
	# cryptic "arguments imply differing number of rows" from `data.frame()`
	# when `.subset2()` returns NULL on an absent column.
	required <- c("cohort", "estimate", "se", "ci_low", "ci_high", "p_value")
	missing_cols <- setdiff(required, names(x))
	if (length(missing_cols) > 0L) {
		stop(
			"tidy.cohortStudy(): input is missing required columns: ",
			paste(missing_cols, collapse = ", "),
			". If you have mutated the `catt_df` slot, restore the original ",
			"data frame; otherwise please file an issue.",
			call. = FALSE
		)
	}
	# Read the source columns via `.subset2()` to bypass the `catt_df`
	# helpful-error S3 layer on `$` / `[[` (which is keyed on old-name
	# access; new-name access falls through harmlessly, but `.subset2()`
	# is the most direct route and documents intent).
	estimate <- .subset2(x, "estimate")
	se <- .subset2(x, "se")
	statistic <- ifelse(
		!is.na(se) & se > 0,
		estimate / se,
		NA_real_
	)
	out <- data.frame(
		term = paste0("cohort_", .subset2(x, "cohort")),
		estimate = estimate,
		std.error = se,
		statistic = statistic,
		p.value = .subset2(x, "p_value"),
		conf.low = .subset2(x, "ci_low"),
		conf.high = .subset2(x, "ci_high"),
		stringsAsFactors = FALSE
	)
	if ("selected" %in% names(x)) {
		out$selected <- .subset2(x, "selected")
	}
	out
}

#' Tidy a `FETWFE_tes` simulation truth object
#'
#' Returns a `broom`-style tidy data frame for the population-truth
#' object returned by [getTes()]. Row 1 is the overall true ATT
#' (`term = "ATT_true"`); subsequent rows are the true cohort ATTs
#' (`term = "Cohort <adoption-time>"`, using the simulator's
#' convention that cohort `r` adopts at calendar time `r + 1`, so
#' the labels match what `tidy.<estimator>` uses on a fitted panel
#' generated from the same `FETWFE_coefs`). Standard error /
#' statistic / p-value columns are always `NA_real_` — there is no
#' sampling distribution for a population truth. When
#' `conf.int = TRUE` (default, matching the sibling tidy methods),
#' `conf.low` / `conf.high` columns are included and also set to
#' `NA_real_`. When `conf.int = FALSE`, those columns are omitted.
#'
#' @param x An object of class `"FETWFE_tes"` returned by [getTes()].
#' @param conf.int Logical; include `conf.low` / `conf.high` columns.
#'   Defaults to `TRUE` to match the sibling tidy methods and preserve
#'   pre-#84 backward compatibility. Population-truth objects have no
#'   sampling distribution, so the CI columns are always filled with
#'   `NA_real_` when included.
#' @param conf.level Numeric in (0, 1). Accepted for broom-convention
#'   parity but unused (no CIs to compute for a population truth);
#'   validated regardless. Defaults to `0.95` (`FETWFE_tes` objects do
#'   not carry an alpha slot, so there is no fitted-object value to
#'   default to).
#' @param ... Unused.
#' @return A data frame with `R + 1` rows and columns `term`,
#'   `estimate`, `std.error`, `statistic`, `p.value`, and (when
#'   `conf.int = TRUE`) `conf.low` / `conf.high`.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   broom::tidy(getTes(coefs))
#' }
#' @export
tidy.FETWFE_tes <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
	stopifnot(conf.level > 0, conf.level < 1)
	R <- length(x$actual_cohort_tes)
	cohort_times <- if (!is.null(x$cohort_times)) {
		x$cohort_times
	} else {
		# Fallback for pre-1.9.0 dev objects without the slot: use the
		# simulator's convention that cohort r adopts at calendar time r + 1.
		as.integer(seq_len(R) + 1L)
	}
	terms <- c("ATT_true", paste0("Cohort ", cohort_times))
	estimates <- c(x$att_true, x$actual_cohort_tes)
	out <- data.frame(
		term = terms,
		estimate = estimates,
		std.error = rep(NA_real_, R + 1L),
		statistic = rep(NA_real_, R + 1L),
		p.value = rep(NA_real_, R + 1L),
		stringsAsFactors = FALSE
	)
	if (conf.int) {
		out$conf.low <- rep(NA_real_, R + 1L)
		out$conf.high <- rep(NA_real_, R + 1L)
	}
	out
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
#' @return A one-row data frame with 16 columns.
#' @examples
#' \dontrun{
#'   res <- fetwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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
#' @return A one-row data frame with 16 columns.
#' @examples
#' \dontrun{
#'   res <- betwfeWithSimulatedData(
#'     simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
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
#' Computes `.fitted = X %*% beta_hat + x$y_mean` and
#' `.resid = data[[x$response_col_name]] - .fitted`, then column-binds those
#' two columns onto `data`. The response mean and column name are stored on
#' the fitted object during fitting (the estimator internally centers `y`
#' before solving), so fitted values come back on the original-response
#' scale without the caller having to remember either.
#'
#' `data` is auto-handled to match the fitted design: rows are auto-sorted
#' by `(unit, time)`, and any first-period-treated units (whose treatment
#' effect cannot be identified by the estimator) are auto-trimmed via
#' `idCohorts()`. So you can pass the same raw `pdata` you handed to
#' `fetwfe()` — the method takes care of alignment. The only hard
#' requirement is that `data` contains the response column under its
#' original name.
#'
#' @param x An object of class `"fetwfe"`.
#' @param data A panel `data.frame` with one row per unit-period (any
#'   sort order — augment auto-sorts), containing the response column
#'   under the same name used at fit time (see `x$response_col_name`).
#'   First-period-treated units, if present, are auto-trimmed.
#' @param ... Unused.
#' @return A copy of `data` with two extra numeric columns: `.fitted`
#'   and `.resid`.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- fetwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata)
#' }
#' @export
augment.fetwfe <- function(x, data, ...) {
	.augment_estimator_output(x, data, ...)
}

#' Augment user-supplied data with fitted values and residuals from an etwfe fit
#'
#' Same shape as [augment.fetwfe()], dispatched on class `"etwfe"`. `data`
#' is auto-sorted by `(unit, time)` and any first-period-treated units
#' are auto-trimmed; pass the same raw `pdata` you handed to `etwfe()`.
#'
#' @param x An object of class `"etwfe"`.
#' @param data A panel `data.frame` with the response column under
#'   `x$response_col_name`. Any sort order; first-period-treated units
#'   are auto-trimmed.
#' @param ... Unused.
#' @return `data` with `.fitted` and `.resid` columns appended.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- etwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata)
#' }
#' @export
augment.etwfe <- function(x, data, ...) {
	.augment_estimator_output(x, data, ...)
}

#' Augment user-supplied data with fitted values and residuals from a betwfe fit
#'
#' Same shape as [augment.fetwfe()], dispatched on class `"betwfe"`. `data`
#' is auto-sorted by `(unit, time)` and any first-period-treated units
#' are auto-trimmed; pass the same raw `pdata` you handed to `betwfe()`.
#'
#' @param x An object of class `"betwfe"`.
#' @param data A panel `data.frame` with the response column under
#'   `x$response_col_name`. Any sort order; first-period-treated units
#'   are auto-trimmed.
#' @param ... Unused.
#' @return `data` with `.fitted` and `.resid` columns appended.
#' @examples
#' \dontrun{
#'   sim <- simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5,
#'                                eff_size = 2),
#'                       N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- betwfeWithSimulatedData(sim)
#'   broom::augment(res, data = sim$pdata)
#' }
#' @export
augment.betwfe <- function(x, data, ...) {
	.augment_estimator_output(x, data, ...)
}
