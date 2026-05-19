#-------------------------------------------------------------------------------
# Shared helpers for fetwfe / etwfe / betwfe S3 class methods
#-------------------------------------------------------------------------------

#' @title Reorder and truncate a CATT data frame for printing
#'
#' @description
#' Internal helper used by `print.fetwfe()`, `print.etwfe()`, `print.betwfe()`
#' and their `summary.*` companions to render the cohort-average-treatment-
#' effect table. Reorders `catt_df` by one of five orderings and truncates to
#' the first `max_cohorts` rows. When truncation occurs, attaches the
#' attributes `truncated = TRUE` and `n_discarded = <number dropped>` so the
#' caller can print a "... and N more cohorts" suffix.
#'
#' Every call site passes `max_cohorts` explicitly, so the parameter is
#' positional-required. The per-class option keys
#' (`fetwfe.max_cohorts` / `etwfe.max_cohorts` / `betwfe.max_cohorts`) are
#' resolved on each `print.<class>()` method's own `max_cohorts =`
#' argument, not here.
#'
#' @param catt_df A data frame with at least the columns `Cohort` and
#'   `Estimated TE`, plus `P_value` if `order_by = "pvalue"` is used.
#' @param max_cohorts Integer; maximum number of cohort rows to retain after
#'   sorting. Required positional argument.
#' @param order_by Character; one of `"cohort"`, `"estimate"`, `"abs_estimate"`,
#'   `"pvalue"`, `"none"`. Defaults to `"cohort"`.
#' @return A data frame, possibly truncated, with attributes `truncated`
#'   (logical) and (if truncated) `n_discarded` (integer).
#' @keywords internal
#' @noRd
.truncate_catt <- function(
	catt_df,
	max_cohorts,
	order_by = c("cohort", "estimate", "abs_estimate", "pvalue", "none")
) {
	order_by <- match.arg(order_by)

	idx <- switch(
		order_by,
		# Composite sort key: numeric-cohort-time first, character as
		# tiebreak. The character fallback only matters if some labels
		# aren't numerically coercible (hypothetical -- catt_df$Cohort is
		# always as.character(integer time) given the package's
		# `is.integer(time_var)` input validation).
		cohort = order(
			suppressWarnings(as.numeric(catt_df$Cohort)),
			catt_df$Cohort
		),
		estimate = order(catt_df$`Estimated TE`),
		abs_estimate = order(abs(catt_df$`Estimated TE`), decreasing = TRUE),
		pvalue = order(catt_df$P_value),
		none = seq_len(nrow(catt_df))
	)
	catt_df <- catt_df[idx, , drop = FALSE]

	if (nrow(catt_df) > max_cohorts) {
		attr(catt_df, "truncated") <- TRUE
		attr(catt_df, "n_discarded") <- nrow(catt_df) - max_cohorts
		catt_df <- catt_df[seq_len(max_cohorts), , drop = FALSE]
	} else {
		attr(catt_df, "truncated") <- FALSE
	}
	catt_df
}

#' @title Print a CATT data frame with the package's standard formatting
#'
#' @description
#' Internal one-line helper that wraps `print()` with `row.names = FALSE` and
#' `right = TRUE`, used by every `print.fetwfe()` / `print.etwfe()` /
#' `print.betwfe()` and their `print.summary.*` companions to render the
#' cohort table consistently. Returns `invisible()` so it can be chained
#' inside a print-method body without producing a second visible printout.
#'
#' @param df A data frame to print.
#' @return `invisible(NULL)`.
#' @keywords internal
#' @noRd
.print_catt_tbl <- function(df) {
	print(df, row.names = FALSE, right = TRUE)
	invisible()
}

#-------------------------------------------------------------------------------
# Constructor-validator helpers (issue #85).
#
# Shared cross-slot contract assertions used by .validate_fetwfe /
# .validate_etwfe / .validate_betwfe / .validate_twfeCovs. Each helper
# checks one contract family and stops with a structured error message
# on violation. Contract violations are programmer bugs (the estimator
# constructed a malformed object) -- they should never reach a user under
# correct usage. The error message includes (a) which class's validator
# fired, (b) which contract failed, (c) actual vs expected when useful.
#-------------------------------------------------------------------------------

#' @title Stop with a structured error if any expected slot is missing
#' @keywords internal
#' @noRd
.stop_if_missing_slots <- function(x, expected, cls, where = NULL) {
	missing_slots <- setdiff(expected, names(x))
	if (length(missing_slots) > 0) {
		prefix <- if (is.null(where)) "" else paste0(" (in `", where, "`)")
		stop(
			".validate_",
			cls,
			"(): malformed `",
			cls,
			"` object",
			prefix,
			". Missing slot(s): ",
			paste(missing_slots, collapse = ", "),
			". This is a programmer-side contract violation -- please ",
			"report at https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
}

#' @title Assert a single contract; stop with a structured message on violation
#' @keywords internal
#' @noRd
.assert_contract <- function(cond, contract_name, cls, detail = "") {
	if (!isTRUE(cond)) {
		stop(
			".validate_",
			cls,
			"(): contract violation `",
			contract_name,
			"`. ",
			detail,
			"\nThis is a programmer-side contract violation -- please ",
			"report at https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
}

#' @title C1 -- SE consistency
#' @description
#' When `calc_ses` is FALSE, the estimator cannot produce valid standard
#' errors, so all SE-bearing slots must be NA. This is the contract whose
#' violation produced issue #73 (`event_study()` reporting finite SEs for
#' q >= 1). `calc_ses_path` is the string path to read `calc_ses` from
#' (e.g., `"calc_ses"` for top-level, `"internal$calc_ses"` for fetwfe).
#' @keywords internal
#' @noRd
.check_se_consistency <- function(x, calc_ses_path, cls) {
	# Resolve the calc_ses value via path-string evaluation (covers both
	# top-level and nested-under-$internal cases).
	calc_ses <- eval(
		parse(text = paste0("x$", calc_ses_path)),
		envir = environment()
	)
	if (!isTRUE(calc_ses)) {
		.assert_contract(
			is.na(x$att_se),
			"C1 SE consistency: att_se must be NA when calc_ses is FALSE",
			cls,
			detail = paste0(
				"att_se = ",
				format(x$att_se),
				", ",
				calc_ses_path,
				" = ",
				format(calc_ses)
			)
		)
		.assert_contract(
			all(is.na(x$catt_ses)),
			"C1 SE consistency: all catt_ses must be NA when calc_ses is FALSE",
			cls
		)
	}
}

#' @title C2 -- Selection consistency (fetwfe / betwfe only)
#' @description
#' `att_selected` is TRUE iff at least one cohort's coefficient was selected
#' (i.e., `att_hat != 0`). When the all-zero fallback fires, `att_hat == 0`
#' AND `att_selected == FALSE`.
#' @keywords internal
#' @noRd
.check_selection_consistency <- function(x, cls) {
	.assert_contract(
		isTRUE(x$att_selected) == isTRUE(x$att_hat != 0),
		"C2 selection consistency: att_selected == (att_hat != 0)",
		cls,
		detail = paste0(
			"att_selected = ",
			format(x$att_selected),
			", att_hat = ",
			format(x$att_hat)
		)
	)
}

#' @title C3 -- p-value NA-derivation consistency
#' @keywords internal
#' @noRd
.check_p_value_na <- function(x, cls) {
	se_is_zero_or_na <- is.na(x$att_se) || isTRUE(x$att_se == 0)
	.assert_contract(
		is.na(x$att_p_value) == se_is_zero_or_na,
		"C3 p-value NA: is.na(att_p_value) iff is.na(att_se) || att_se == 0",
		cls,
		detail = paste0(
			"att_se = ",
			format(x$att_se),
			", att_p_value = ",
			format(x$att_p_value)
		)
	)
}

#' @title C4 -- catt_df shape
#' @keywords internal
#' @noRd
.check_catt_df_shape <- function(x, cls) {
	.assert_contract(
		is.data.frame(x$catt_df),
		"C4 catt_df: must be a data.frame",
		cls
	)
	.assert_contract(
		nrow(x$catt_df) == x$R,
		"C4 catt_df: nrow == R",
		cls,
		detail = paste0(
			"nrow(catt_df) = ",
			nrow(x$catt_df),
			", R = ",
			x$R
		)
	)
	.assert_contract(
		"Cohort" %in% names(x$catt_df) && is.character(x$catt_df$Cohort),
		"C4 catt_df: Cohort column present and character",
		cls
	)
}

#' @title C5 -- Cohort-probability structural sanity
#' @keywords internal
#' @noRd
.check_cohort_probs <- function(x, cls) {
	.assert_contract(
		length(x$cohort_probs) == x$R,
		"C5 cohort_probs length == R",
		cls,
		detail = paste0(
			"length(cohort_probs) = ",
			length(x$cohort_probs),
			", R = ",
			x$R
		)
	)
	.assert_contract(
		length(x$cohort_probs_overall) == x$R,
		"C5 cohort_probs_overall length == R",
		cls
	)
	.assert_contract(
		all(x$cohort_probs >= 0 & x$cohort_probs <= 1 + 1e-10),
		"C5 cohort_probs in [0, 1]",
		cls
	)
	.assert_contract(
		all(x$cohort_probs_overall >= 0 & x$cohort_probs_overall <= 1 + 1e-10),
		"C5 cohort_probs_overall in [0, 1]",
		cls
	)
	# Strict < 1 - 1e-6 matches the existing guard at R/fetwfe_core.R:2025.
	.assert_contract(
		sum(x$cohort_probs_overall) < 1 - 1e-6,
		"C5 sum(cohort_probs_overall) < 1 - 1e-6",
		cls,
		detail = paste0(
			"sum = ",
			format(sum(x$cohort_probs_overall), digits = 10)
		)
	)
}

#' @title C7 -- Lambda monotonicity (fetwfe / betwfe only)
#' @description
#' Larger lambda -> smaller model, so model-size direction is REVERSED from
#' lambda direction: `lambda.max_model_size <= lambda_star_model_size <=
#' lambda.min_model_size`. The upper bound for `lambda.min_model_size` is
#' `p + 1` (intercept counted in grpreg's beta-path).
#' @keywords internal
#' @noRd
.check_lambda_monotonicity <- function(x, cls) {
	.assert_contract(
		is.finite(x$lambda.max) && x$lambda.max > 0,
		"C7 lambda.max finite and > 0",
		cls
	)
	.assert_contract(
		is.finite(x$lambda.min) && x$lambda.min >= 0,
		"C7 lambda.min finite and >= 0",
		cls
	)
	.assert_contract(
		x$lambda.max >= x$lambda.min,
		"C7 lambda.max >= lambda.min",
		cls
	)
	.assert_contract(
		x$lambda_star >= x$lambda.min && x$lambda_star <= x$lambda.max,
		"C7 lambda_star in [lambda.min, lambda.max]",
		cls
	)
	# Model size REVERSED from lambda direction.
	.assert_contract(
		x$lambda.max_model_size <= x$lambda_star_model_size,
		"C7 lambda.max_model_size <= lambda_star_model_size",
		cls
	)
	.assert_contract(
		x$lambda_star_model_size <= x$lambda.min_model_size,
		"C7 lambda_star_model_size <= lambda.min_model_size",
		cls
	)
	.assert_contract(
		x$lambda.min_model_size <= x$p + 1L,
		"C7 lambda.min_model_size <= p + 1",
		cls
	)
}

#' @title C8 -- Type sanity (universal across all 4 classes)
#' @description
#' `has_alpha` is TRUE for fetwfe/etwfe/betwfe (which all carry an `alpha`
#' slot); FALSE for twfeCovs (which omits alpha because it has no inference
#' output). `has_att_selected` is TRUE for fetwfe/betwfe.
#' @keywords internal
#' @noRd
.check_type_sanity <- function(
	x,
	cls,
	has_alpha = TRUE,
	has_att_selected = FALSE
) {
	for (slot in c("N", "T", "R", "p")) {
		val <- x[[slot]]
		.assert_contract(
			is.numeric(val) &&
				length(val) == 1L &&
				val > 0L &&
				val == round(val),
			paste0("C8 ", slot, " is positive integer-valued scalar"),
			cls,
			detail = paste0(slot, " = ", format(val))
		)
	}
	# d may be 0 (no covariates).
	.assert_contract(
		is.numeric(x$d) && length(x$d) == 1L && x$d >= 0L && x$d == round(x$d),
		"C8 d is non-negative integer-valued scalar",
		cls,
		detail = paste0("d = ", format(x$d))
	)
	for (slot in c("time_var", "unit_var", "treatment", "response_col_name")) {
		val <- x[[slot]]
		.assert_contract(
			is.character(val) && length(val) == 1L,
			paste0("C8 ", slot, " is length-1 character"),
			cls
		)
	}
	# covs is NULL when d == 0, character otherwise.
	.assert_contract(
		is.null(x$covs) || is.character(x$covs),
		"C8 covs is NULL or character",
		cls
	)
	.assert_contract(
		x$se_type %in% c("default", "cluster"),
		"C8 se_type in c('default', 'cluster')",
		cls,
		detail = paste0("se_type = ", format(x$se_type))
	)
	.assert_contract(
		is.logical(x$indep_counts_used) && length(x$indep_counts_used) == 1L,
		"C8 indep_counts_used is length-1 logical",
		cls
	)
	.assert_contract(
		is.numeric(x$y_mean) && length(x$y_mean) == 1L,
		"C8 y_mean is length-1 numeric",
		cls
	)
	.assert_contract(
		is.numeric(x$treat_inds) && length(x$treat_inds) > 0L,
		"C8 treat_inds is non-empty numeric/integer",
		cls
	)
	# treat_int_inds is NULL when d == 0.
	.assert_contract(
		is.null(x$treat_int_inds) || is.numeric(x$treat_int_inds),
		"C8 treat_int_inds is NULL or numeric/integer",
		cls
	)
	if (has_alpha) {
		.assert_contract(
			is.numeric(x$alpha) &&
				length(x$alpha) == 1L &&
				x$alpha > 0 &&
				x$alpha < 1,
			"C8 alpha in (0, 1)",
			cls,
			detail = paste0("alpha = ", format(x$alpha))
		)
	}
	if (has_att_selected) {
		.assert_contract(
			is.logical(x$att_selected) && length(x$att_selected) == 1L,
			"C8 att_selected is length-1 logical",
			cls
		)
	}
}

#-------------------------------------------------------------------------------
# Method-entry preconditions (issue #86).
#
# Each downstream method that reads from a fitted estimator object calls
# a .check_for_<method>(x) helper at its top. The helper:
#   (1) Re-runs the constructor validator from #85 (defense-in-depth;
#       the object may have been hand-modified between construction and
#       method call).
#   (2) Returns a small named list of method-relevant invariants the
#       method can use rather than re-deriving them.
#
# This is what fixes #73: .check_for_event_study(x) returns
# `has_valid_ses` derived from `x$internal$calc_ses` (fetwfe) or
# `x$calc_ses` (etwfe/betwfe), and the event_study dispatchers
# AND-gate their SE-computation branch on it.
#-------------------------------------------------------------------------------

#' @title Universal dispatcher: validate any estimator-class object
#' @description Dispatches via `inherits()` to the appropriate
#' `.validate_<class>` helper from #85. As of #76, `twfeCovs` is now
#' a classed list and its `.validate_twfeCovs` branch is wired up.
#' @keywords internal
#' @noRd
.assert_estimator_object <- function(x) {
	if (inherits(x, "fetwfe")) {
		.validate_fetwfe(x)
	} else if (inherits(x, "etwfe")) {
		.validate_etwfe(x)
	} else if (inherits(x, "betwfe")) {
		.validate_betwfe(x)
	} else if (inherits(x, "twfeCovs")) {
		.validate_twfeCovs(x)
	} else {
		stop(
			"Expected a `fetwfe`, `etwfe`, `betwfe`, or `twfeCovs` object; got class(es): ",
			paste(class(x), collapse = ", "),
			call. = FALSE
		)
	}
	invisible(x)
}

#' @title Method precondition: event_study
#' @description Validates the input + derives `has_valid_ses` (the
#' contract gate that fixes #73). `calc_ses` lives in different paths
#' across classes: nested under `$internal` for fetwfe; top-level for
#' etwfe/betwfe.
#' @return list(has_valid_ses = logical).
#' @keywords internal
#' @noRd
.check_for_event_study <- function(x) {
	.assert_estimator_object(x)
	calc_ses <- if (inherits(x, "fetwfe")) {
		x$internal$calc_ses
	} else {
		x$calc_ses
	}
	list(has_valid_ses = isTRUE(calc_ses))
}

#' @title Method precondition: augment
#' @keywords internal
#' @noRd
.check_for_augment <- function(x) {
	.assert_estimator_object(x)
	invisible(x)
}

#' @title Method precondition: tidy
#' @keywords internal
#' @noRd
.check_for_tidy <- function(x) {
	.assert_estimator_object(x)
	invisible(x)
}

#' @title Method precondition: glance
#' @keywords internal
#' @noRd
.check_for_glance <- function(x) {
	.assert_estimator_object(x)
	invisible(x)
}

#' @title Method precondition: plot
#' @keywords internal
#' @noRd
.check_for_plot <- function(x) {
	.assert_estimator_object(x)
	invisible(x)
}

#' @title Method precondition: coef
#' @keywords internal
#' @noRd
.check_for_coef <- function(x) {
	.assert_estimator_object(x)
	invisible(x)
}

#-------------------------------------------------------------------------------
# Consolidated print/summary/print.summary bodies for the three S3 classes
# (`fetwfe`, `etwfe`, `betwfe`). Each per-class method in R/<class>_class.R
# is now a thin wrapper that pre-resolves the per-class options
# (header text, gating flags, path-extractor functions, max_cohorts /
# order_by) and delegates here. Issue #77 step 2 / PR following PR #90.
#
# Byte-identity of the output is enforced by the six snapshots at
# tests/testthat/_snaps/print-method-snapshot.md (PR #90).
#-------------------------------------------------------------------------------

#' @title Consolidated body of `print.<class>` for the three estimator classes
#'
#' @description
#' Internal helper called by `print.fetwfe()`, `print.etwfe()`, and
#' `print.betwfe()`. The three per-class methods differ only in (a) the
#' header banner, (b) whether the `Selected:` row appears under the
#' Overall-ATT block, (c) whether the `Selected size` / `Lambda*` rows
#' appear in Model Details, and (d) where `X_ints`, `y`, and `calc_ses`
#' live on the fit object (top-level for `etwfe` / `betwfe`, nested under
#' `$internal` for `fetwfe`). All four are passed in by the caller; the
#' helper renders the full output.
#'
#' @param x The fit object (one of `"fetwfe"` / `"etwfe"` / `"betwfe"`).
#' @param header Single string containing the title line plus its `====`
#'   divider plus two trailing newlines (e.g.
#'   `"Foo Results\n===\n\n"`). Passed verbatim to `cat()`. The
#'   title-versus-divider length relationship is inconsistent across the
#'   three classes, so the helper does not compute the divider --
#'   each per-class wrapper passes its full banner as a literal string.
#' @param show_att_selected Logical; if `TRUE`, render the
#'   `Selected:` row under the Overall ATT block. `fetwfe` / `betwfe`
#'   pass `TRUE`; `etwfe` passes `FALSE` (etwfe is pure OLS with no
#'   selection).
#' @param show_lambda Logical; if `TRUE`, render `Selected size:` and
#'   `Lambda*:` rows under Model Details. Same pattern as
#'   `show_att_selected`.
#' @param X_ints_path,y_path,calc_ses_path Path-extractor functions that
#'   read `X_ints`, `y`, and `calc_ses` off `x`. For `fetwfe` these
#'   return `x$internal$X_ints` / `x$internal$y` / `x$internal$calc_ses`;
#'   for `etwfe` / `betwfe` they return the top-level fields. Functions
#'   are passed (not path strings) for readability.
#' @param max_cohorts Integer; pre-resolved
#'   `getOption("<class>.max_cohorts", 10)` value. The wrapper resolves
#'   the per-class option key (`fetwfe.max_cohorts` etc.) and passes
#'   the integer.
#' @param order_by Character; pre-validated (via `match.arg()` in the
#'   wrapper) one of `"cohort"`, `"estimate"`, `"abs_estimate"`,
#'   `"pvalue"`, `"none"`.
#' @param show_internal Logical; if `TRUE`, render the `Internal Details:`
#'   block (X dims / y length / SEs computed). Passed through from the
#'   wrapper's `show_internal` argument.
#' @param ... Currently ignored (accepted for S3-method-arg compatibility).
#' @return `invisible(x)`.
#' @keywords internal
#' @noRd
.print_estimator_output <- function(
	x,
	header,
	show_att_selected,
	show_lambda,
	X_ints_path,
	y_path,
	calc_ses_path,
	max_cohorts,
	order_by,
	show_internal,
	...
) {
	cat(header)

	## Overall ATT
	ci_pct <- 100 * (1 - x$alpha)
	ci_low <- x$att_hat - qnorm(1 - x$alpha / 2) * x$att_se
	ci_high <- x$att_hat + qnorm(1 - x$alpha / 2) * x$att_se
	cat(sprintf(
		"Overall Average Treatment Effect (ATT):\n  Estimate:   %.4f\n",
		x$att_hat
	))
	if (identical(x$se_type, "cluster")) {
		cat(sprintf(
			"  Std. Error (cluster-robust): %.4f\n",
			x$att_se
		))
	} else {
		cat(sprintf("  Std. Error: %.4f\n", x$att_se))
	}
	if (!is.null(x$att_p_value) && !is.na(x$att_p_value)) {
		cat(sprintf("  P-value:    %.4g\n", x$att_p_value))
	} else {
		cat("  P-value:    NA\n")
	}
	if (show_att_selected) {
		cat(sprintf("  Selected:   %s\n", x$att_selected))
	}
	cat(sprintf(
		"  %.0f%% CI:    [%.4f, %.4f]\n\n",
		ci_pct,
		ci_low,
		ci_high
	))

	## Cohort effects
	catt_df <- .truncate_catt(x$catt_df, max_cohorts, order_by)
	cat("Cohort Average Treatment Effects (CATT):\n")
	.print_catt_tbl(catt_df)
	if (isTRUE(attr(catt_df, "truncated"))) {
		cat(sprintf(
			"  ... and %d more cohorts.\n",
			attr(catt_df, "n_discarded")
		))
	}
	cat("\n")

	## Model info
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", x$N))
	cat(sprintf("  Time periods (T)    : %d\n", x$T))
	cat(sprintf("  Treated cohorts (R) : %d\n", x$R))
	cat(sprintf("  Covariates (d)      : %d\n", x$d))
	cat(sprintf("  Features (p)        : %d\n", x$p))
	if (show_lambda) {
		cat(sprintf("  Selected size       : %d\n", x$lambda_star_model_size))
		cat(sprintf("  Lambda*             : %.4f\n", x$lambda_star))
	}

	if (show_internal) {
		cat("\nInternal Details:\n")
		cat(
			"  X dims     :",
			paste(dim(X_ints_path(x)), collapse = " x "),
			"\n"
		)
		cat("  y length   :", length(y_path(x)), "\n")
		cat("  SEs computed:", calc_ses_path(x), "\n")
	}

	invisible(x)
}

#' @title Consolidated body of `summary.<class>` for the three estimator classes
#'
#' @description
#' Internal helper called by `summary.fetwfe()`, `summary.etwfe()`, and
#' `summary.betwfe()`. Builds the summary list whose `print` method is
#' dispatched on `output_class`. The three per-class summaries differ in
#' (a) the assigned class name (`"summary.fetwfe"` / `"summary.etwfe"` /
#' `"summary.betwfe"`), (b) whether `att_selected` appears between
#' `att` and `catt` (fetwfe / betwfe only), and (c) whether
#' `lambda_star` / `model_size` appear inside `model_info` between
#' `p` and `sig_eps_sq` (fetwfe / betwfe only).
#'
#' **List-key ordering is contractual.** `att_selected` (when present)
#' must come between `att` and `catt`; `lambda_star` / `model_size`
#' (when present) must come between `p` and `sig_eps_sq` inside
#' `model_info`. The snapshot tests do not catch ordering drift
#' (they access by name, not by position), but downstream users who
#' deparse the structure or write order-sensitive tests would see it.
#' The implementation builds the full union list and reorders via
#' a name-vector selection to keep a single source of truth.
#'
#' @param object The fit object.
#' @param output_class Character; one of `"summary.fetwfe"`,
#'   `"summary.etwfe"`, `"summary.betwfe"`. Set as the class of
#'   the returned list.
#' @param include_att_selected Logical; if `TRUE`, include
#'   `att_selected = object$att_selected` between `att` and
#'   `catt`. fetwfe / betwfe pass `TRUE`; etwfe passes `FALSE`.
#' @param include_lambda Logical; if `TRUE`, include
#'   `lambda_star` and `model_size` between `p` and `sig_eps_sq`
#'   inside `model_info`. Same gating as `include_att_selected`.
#' @param full_catt Logical; passed through from the wrapper's
#'   `full_catt` argument. If `TRUE`, the un-truncated `catt_df`
#'   appears; otherwise the helper truncates to the first 20
#'   rows.
#' @return A list with class set to `output_class`.
#' @keywords internal
#' @noRd
.summary_estimator_output <- function(
	object,
	output_class,
	include_att_selected,
	include_lambda,
	full_catt
) {
	# Build the FULL union list with every possible field, then drop
	# the fields that don't belong to this class via name-vector
	# selection. Single source of truth for field contents; the
	# `keep` vector encodes the contractual ordering.
	#
	# The per-summary preview truncation defaults to `max_cohorts = 20`
	# (kept verbatim from the pre-refactor `summary.<class>` bodies);
	# the print-side default is 10. Summary keeps a larger preview
	# because it is a one-shot interactive inspection rather than a
	# screen-formatted layout.
	out <- list(
		att = c(
			estimate = object$att_hat,
			se = object$att_se,
			p_value = object$att_p_value
		),
		att_selected = object$att_selected,
		catt = if (full_catt) {
			object$catt_df
		} else {
			.truncate_catt(object$catt_df, max_cohorts = 20)
		},
		model_info = list(
			N = object$N,
			T = object$T,
			R = object$R,
			d = object$d,
			p = object$p,
			lambda_star = object$lambda_star,
			model_size = object$lambda_star_model_size,
			sig_eps_sq = object$sig_eps_sq,
			sig_eps_c_sq = object$sig_eps_c_sq
		),
		alpha = object$alpha,
		se_type = object$se_type
	)

	keep <- c(
		"att",
		if (include_att_selected) "att_selected",
		"catt",
		"model_info",
		"alpha",
		"se_type"
	)
	out <- out[keep]

	out$model_info <- out$model_info[c(
		"N",
		"T",
		"R",
		"d",
		"p",
		if (include_lambda) c("lambda_star", "model_size"),
		"sig_eps_sq",
		"sig_eps_c_sq"
	)]

	structure(out, class = output_class)
}

#' @title Consolidated body of `print.summary.<class>` for the three classes
#'
#' @description
#' Internal helper called by `print.summary.fetwfe()`,
#' `print.summary.etwfe()`, and `print.summary.betwfe()`. The three
#' per-class methods differ only in (a) the header banner, (b) whether
#' the `Selected:` line appears under the one-line ATT, and (c) whether
#' the `Selected size` / `Lambda*` rows appear under Model Details.
#' All three are passed in by the caller.
#'
#' Reads `lambda_star` / `model_size` from `x$model_info$...` (the
#' nested summary-list shape produced by
#' `.summary_estimator_output()`), NOT from the top-level
#' `lambda_star_model_size` / `lambda_star` slots that
#' `.print_estimator_output()` reads.
#'
#' @param x A summary object (one of `"summary.fetwfe"` /
#'   `"summary.etwfe"` / `"summary.betwfe"`).
#' @param header Single string containing the title plus its
#'   `===` divider plus two trailing newlines, verbatim. Same
#'   shape as in `.print_estimator_output()`.
#' @param show_att_selected Logical; if `TRUE`, render
#'   `Selected: <value>` after the one-line Overall ATT. fetwfe /
#'   betwfe pass `TRUE`; etwfe passes `FALSE`. Spacing matches the
#'   existing snapshots: one blank line after the Overall ATT line
#'   in all three cases.
#' @param show_lambda Logical; if `TRUE`, render `Selected size:` and
#'   `Lambda*:` rows under Model Details.
#' @return `invisible(x)`.
#' @keywords internal
#' @noRd
.print_summary_estimator_output <- function(
	x,
	header,
	show_att_selected,
	show_lambda
) {
	cat(header)

	ci_pct <- 100 * (1 - x$alpha)
	ci_low <- x$att["estimate"] - qnorm(1 - x$alpha / 2) * x$att["se"]
	ci_high <- x$att["estimate"] + qnorm(1 - x$alpha / 2) * x$att["se"]
	p_val <- x$att["p_value"]
	p_str <- if (is.na(p_val)) "NA" else sprintf("%.4g", p_val)
	se_label <- if (identical(x$se_type, "cluster")) {
		"SE (cluster-robust)"
	} else {
		"SE"
	}
	cat(sprintf(
		"Overall ATT: %.4f  (%s = %.4f, p = %s, %.0f%% CI = [%.4f, %.4f])\n",
		x$att["estimate"],
		se_label,
		x$att["se"],
		p_str,
		ci_pct,
		ci_low,
		ci_high
	))
	if (show_att_selected) {
		cat(sprintf("Selected: %s\n\n", x$att_selected))
	} else {
		cat("\n")
	}

	cat("CATT (preview):\n")
	.print_catt_tbl(x$catt)
	if (isTRUE(attr(x$catt, "truncated"))) {
		cat(sprintf("  ... + %d more cohorts.\n", attr(x$catt, "n_discarded")))
	}
	cat("\n")

	## Model info
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", x$model_info$N))
	cat(sprintf("  Time periods (T)    : %d\n", x$model_info$T))
	cat(sprintf("  Treated cohorts (R) : %d\n", x$model_info$R))
	cat(sprintf("  Covariates (d)      : %d\n", x$model_info$d))
	cat(sprintf("  Features (p)        : %d\n", x$model_info$p))
	if (show_lambda) {
		cat(sprintf("  Selected size       : %d\n", x$model_info$model_size))
		cat(sprintf("  Lambda*             : %.4f\n", x$model_info$lambda_star))
	}

	invisible(x)
}
