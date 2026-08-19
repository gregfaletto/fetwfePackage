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
#' @param catt_df A data frame with at least the columns `cohort` and
#'   `estimate`, plus `p_value` if `order_by = "pvalue"` is used.
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
		# aren't numerically coercible (hypothetical -- catt_df$cohort is
		# always as.character(integer time) given the package's
		# `is.integer(time_var)` input validation).
		cohort = order(
			suppressWarnings(as.numeric(catt_df$cohort)),
			catt_df$cohort
		),
		estimate = order(catt_df$estimate),
		abs_estimate = order(abs(catt_df$estimate), decreasing = TRUE),
		pvalue = order(catt_df$p_value),
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

#' @title Truncate an eventStudy frame to a max number of event times
#'
#' @description
#' Internal helper used by `.print_estimator_output()` and
#' `.summary_estimator_output()` to render a preview of event-study
#' effects matching the existing CATT-preview style. Mirrors
#' `.truncate_catt()` in shape (preserves class via `head()`-style
#' subsetting, sets `truncated` / `n_discarded` attributes), but
#' simpler -- `eventStudy()` already returns rows in `event_time`
#' order so no sort step is needed.
#'
#' @param es An `eventStudy`-classed data frame.
#' @param max_event_times Integer; the maximum number of rows to
#'   preserve. If `nrow(es) > max_event_times`, the head is kept and
#'   the rest discarded; otherwise `es` is returned unchanged.
#' @return A data frame, possibly truncated, with attributes
#'   `truncated` (logical) and (if truncated) `n_discarded` (integer).
#' @keywords internal
#' @noRd
.truncate_event_study <- function(es, max_event_times) {
	if (nrow(es) > max_event_times) {
		attr_truncated <- TRUE
		n_discarded <- nrow(es) - max_event_times
		es <- es[seq_len(max_event_times), , drop = FALSE]
	} else {
		attr_truncated <- FALSE
		n_discarded <- 0L
	}
	attr(es, "truncated") <- attr_truncated
	if (attr_truncated) {
		attr(es, "n_discarded") <- n_discarded
	}
	es
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
#' violation produced issue #73 (`eventStudy()` reporting finite SEs for
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
		nrow(x$catt_df) == x$G,
		"C4 catt_df: nrow == G",
		cls,
		detail = paste0(
			"nrow(catt_df) = ",
			nrow(x$catt_df),
			", G = ",
			x$G
		)
	)
	.assert_contract(
		"cohort" %in% names(x$catt_df) && is.character(x$catt_df$cohort),
		"C4 catt_df: cohort column present and character",
		cls
	)
	.assert_contract(
		inherits(x$catt_df, "catt_df"),
		"C4 catt_df: object carries the `catt_df` S3 class",
		cls
	)
}

#' @title Resolve the alpha a fit's catt_df was built at
#' @description Returns `x$alpha` for all four classes
#'   (fetwfe/etwfe/betwfe/twfeCovs); every fit now carries an `alpha`
#'   slot recording the level its catt_df band was built at (#204).
#' @keywords internal
#' @noRd
.alpha_of <- function(x) {
	x$alpha
}

#' @title C10 -- Simultaneous band is never narrower than the pointwise band
#' @description
#' Encodes the core invariant of the `ci_type = "simultaneous"` default (#197):
#' a family-wise (uniform) band is never narrower than the per-effect pointwise
#' Wald band. Gated to run only when meaningful: `ci_type == "simultaneous"`,
#' standard errors available, at least two cohorts (no widening when `K = 1`),
#' and all `se` finite and strictly positive (degenerate / selected-out rows
#' are exempt). For `ci_type == "pointwise"` the check is skipped entirely, and
#' the gate is satisfied-by-equality on pre-overwrite pointwise bounds, so the
#' double-validation (pre-classing pointwise + post-finalizer simultaneous)
#' never produces a false positive.
#' @keywords internal
#' @noRd
.check_ci_band_width <- function(x, cls) {
	if (!identical(x$ci_type, "simultaneous")) {
		return(invisible(NULL))
	}
	cd <- x$catt_df
	se <- cd$se
	if (
		isTRUE(nrow(cd) >= 2L) &&
			all(is.finite(se)) &&
			all(se > 0)
	) {
		z <- stats::qnorm(1 - .alpha_of(x) / 2)
		pw_width <- 2 * z * se
		sim_width <- cd$ci_high - cd$ci_low
		.assert_contract(
			all(sim_width >= pw_width - 1e-6 * pmax(pw_width, 1)),
			"C10 simultaneous band >= pointwise band",
			cls
		)
	}
	invisible(NULL)
}

#' @title C5 -- Cohort-probability structural sanity
#' @keywords internal
#' @noRd
.check_cohort_probs <- function(x, cls) {
	.assert_contract(
		length(x$cohort_probs) == x$G,
		"C5 cohort_probs length == G",
		cls,
		detail = paste0(
			"length(cohort_probs) = ",
			length(x$cohort_probs),
			", G = ",
			x$G
		)
	)
	.assert_contract(
		length(x$cohort_probs_overall) == x$G,
		"C5 cohort_probs_overall length == G",
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
	# Strict < 1 - 1e-6 matches the existing guard inside
	# `getSecondVarTermDataApp()` (R/variance_machinery.R).
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

#' @title Top-level C6 dimension + C8 calc_ses contract checks (class-invariant)
#' @description The C6 (length / nrow) and C8 (`calc_ses`) top-level contract
#'   checks shared byte-identically by the etwfe / betwfe / twfeCovs validators.
#'   Single-sourced (#401) so a future C6 invariant added to one class cannot
#'   silently skip the others. fetwfe's variant is nested and adds C11, so it
#'   stays separate.
#' @param x A fitted object.
#' @param cls The object's class label, used in the assertion messages.
#' @keywords internal
#' @noRd
.check_c6_dims_toplevel <- function(x, cls) {
	# #431 item 6: read every slot with `[[`, never `$`. `$` partial-matches on
	# lists, so a slot that is absent can silently resolve to a longer-named one
	# -- the trap tests/testthat/test-internal-slot-parity.R:88-90 documents and
	# mandates `[[` for. All seven reads are converted, not just the four the
	# issue named, so no `$` remains in this function.
	#
	# Be precise about the hazard, because the obvious example does not apply
	# here: on the three classes this function serves, an absent `y` does NOT
	# resolve to `y_mean`, because `y_final` is also present and the prefix is
	# therefore ambiguous -- `$` returns NULL just as `[[` does. The reads that
	# can diverge are the ones with exactly one longer-named sibling, which is
	# why the tests build one fixture per read (test-small-fixes-431.R item 6).
	# So this is behavior-neutral for well-formed etwfe / betwfe / twfeCovs
	# objects, which carry all seven slots exactly; it matters for a malformed
	# object and for the class-agnostic future caller this helper's `_toplevel`
	# name anticipates.
	#
	# NOT the only contract-checker in this file, and the others still use `$`
	# with live prefix collisions (cohort_probs -> cohort_probs_overall on all
	# four classes, att_se -> att_selected, three lambda.* pairs). Those are
	# latent because .stop_if_missing_slots() runs an exact setdiff first on
	# every validator path; converting them is its own change, filed separately.
	.assert_contract(
		length(x[["beta_hat"]]) == x[["p"]],
		"C6 length(beta_hat) == p",
		cls
	)
	.assert_contract(
		length(x[["y"]]) == x[["N"]] * x[["T"]],
		"C6 length(y) == N * T",
		cls
	)
	.assert_contract(
		nrow(x[["X_ints"]]) == x[["N"]] * x[["T"]],
		"C6 nrow(X_ints) == N * T",
		cls
	)
	.assert_contract(
		is.logical(x[["calc_ses"]]) && length(x[["calc_ses"]]) == 1L,
		"C8 calc_ses is length-1 logical",
		cls
	)
}

#' @title C7 -- Lambda monotonicity (fetwfe / betwfe only)
#' @description
#' Larger lambda -> smaller model, so model-size direction is REVERSED from
#' lambda direction: `lambda.max_model_size <= lambda_star_model_size <=
#' lambda.min_model_size`. The upper bound for `lambda.min_model_size` is
#' `p`, since the reported sizes count selected features and exclude the
#' always-present intercept (#269).
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
		x$lambda.min_model_size <= x$p,
		"C7 lambda.min_model_size <= p",
		cls
	)
}

#' @title C8 -- Type sanity (universal across all 4 classes)
#' @description
#' `has_alpha` is TRUE for all four classes (fetwfe/etwfe/betwfe/twfeCovs all
#' carry an `alpha` slot; twfeCovs gained one in #204). `has_att_selected` is
#' TRUE for fetwfe/betwfe.
#' @keywords internal
#' @noRd
.check_type_sanity <- function(
	x,
	cls,
	has_alpha = TRUE,
	has_att_selected = FALSE
) {
	for (slot in c("N", "T", "G", "p")) {
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
	# C8 -- the canonical cohort-count field `G` is dual-populated with the
	# same value as the deprecated `R` alias (#41).
	.assert_contract(
		identical(as.numeric(x$G), as.numeric(x$R)),
		"C8 G equals R (dual-populated cohort count)",
		cls,
		detail = paste0("G = ", format(x$G), ", R = ", format(x$R))
	)
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
	# A d == 0 fit from the shipped entry points stores character(0), NOT NULL;
	# NULL is allowed here only for hand-built or legacy objects. (#403)
	.assert_contract(
		is.null(x$covs) || is.character(x$covs),
		"C8 covs is NULL or character",
		cls
	)
	.assert_contract(
		x$se_type %in% c("default", "conservative", "cluster"),
		"C8 se_type in c('default', 'conservative', 'cluster')",
		cls,
		detail = paste0("se_type = ", format(x$se_type))
	)
	# C9 (#197): ci_type slot sanity. Universal across all four classes
	# (every fit produced by version >= 1.16.0 carries a `ci_type` slot).
	.assert_contract(
		length(x$ci_type) == 1L &&
			x$ci_type %in% c("simultaneous", "pointwise"),
		"C9 ci_type in c('simultaneous', 'pointwise')",
		cls,
		detail = paste0("ci_type = ", format(x$ci_type))
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
# `x$calc_ses` (etwfe/betwfe), and the .event_study_* dispatchers
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

#' @title Method precondition: eventStudy
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

#' @title Method precondition: simultaneousCIs
#' @description Validates the input + derives `has_valid_ses` (the contract
#' gate). Mirrors `.check_for_event_study(x)`. `calc_ses` lives in different
#' paths across classes: nested under `$internal` for fetwfe; top-level for
#' etwfe/betwfe/twfeCovs.
#' @return list(has_valid_ses = logical).
#' @keywords internal
#' @noRd
.check_for_simultaneous_cis <- function(x) {
	.assert_estimator_object(x)
	calc_ses <- if (inherits(x, "fetwfe")) {
		x$internal$calc_ses
	} else {
		x$calc_ses
	}
	list(has_valid_ses = isTRUE(calc_ses))
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
# Byte-identity of the output is enforced by the snapshots at
# tests/testthat/_snaps/print-method-snapshot.md (PR #90, extended by #439).
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Shared rendering pieces (#439, issue #366 part A2).
#
# `print()` and `print(summary())` are meant to report the same facts in the same
# layout. Before #439 each of the five pieces below was written twice, once per
# function, kept in agreement only by snapshot tests that lock each copy
# separately -- so an edit to one could land without the other and the package
# would disagree with itself about, say, whether its intervals are simultaneous
# or pointwise. Each is now emitted from one place.
#
# The callers keep the event-study GATE (`print` gates on a live `eventStudy(x)`
# call, `summary` on the cached `x$event_study`); only the rendering skeleton is
# shared. Absorbing the gate would collapse the two independent facts that
# test-event-study-present-in-print-summary-174.R asserts into one.
#-------------------------------------------------------------------------------

#' @title Parenthetical naming the standard-error flavor
#'
#' @description
#' Returns the qualifier the print and summary paths append after the words
#' "Std. Error" / "SE" -- `" (cluster-robust)"`, `" (conservative)"`, or `""`
#' for the default. The helper deliberately stops at the qualifier rather than
#' returning a finished label, because **the two callers' literals differ**:
#' `.print_estimator_output()` renders `"  Std. Error (cluster-robust): %.4f\n"`
#' while `.print_summary_estimator_output()` renders `"SE (cluster-robust)"`.
#' One shared finished label could not serve both without changing user-visible
#' output, so each caller composes its own surrounding text.
#'
#' @param se_type Character scalar; the fit's `se_type` slot, one of
#'   `"default"`, `"conservative"`, `"cluster"`.
#' @return A length-1 character: the parenthetical including its leading space,
#'   or `""` for anything other than `"cluster"` / `"conservative"`.
#' @keywords internal
#' @noRd
.se_qualifier <- function(se_type) {
	if (identical(se_type, "cluster")) {
		" (cluster-robust)"
	} else if (identical(se_type, "conservative")) {
		" (conservative)"
	} else {
		""
	}
}

#' @title Band-type label for a confidence-interval table header
#'
#' @description
#' Returns `"simultaneous"` or `"pointwise"` for the `[<label> NN% CI]` suffix on
#' the CATT and event-study preview headers (#197).
#'
#' Takes the **scalar** `ci_type`, never the object: the print path holds a
#' `fetwfe`-family fit and the summary path a `summary.fetwfe`-family list, so a
#' single object parameter would carry two unrelated shapes and could validate
#' neither.
#'
#' @details
#' `NULL` -- a pre-1.16.0 fit or summary carrying no `ci_type` slot -- returns
#' `"pointwise"`. That is a **compatibility behavior, not a default**: those
#' objects' stored bounds really are pointwise, so the label is accurate.
#'
#' The same compatibility rule is encoded independently by
#' `.resolve_event_study_ci_type()` in `R/event_study.R`, with different control
#' flow and a divergence on invalid input (this helper falls through to
#' `"pointwise"`; that one raises a `match.arg()` error). If the rule is ever
#' revisited, both must change.
#'
#' Not to be confused with the `identical(x$ci_type, "simultaneous")` tests in
#' `.check_ci_band_width()` (this file), in `.highdim_postselection_band_notice()`
#' (this file, just below), in `R/event_study.R`, and in `.finalize_ci_type()`
#' (`R/simultaneous_cis.R`). Those four are **gates**, not labels; routing a gate
#' through this helper would make a control-flow decision depend on a display
#' string. (The fourth arrived with #433 and is the closest call of the four,
#' since it gates a *string* the way this helper returns one -- but it decides
#' whether a caveat is rendered at all, which is control flow.)
#'
#' @param ci_type Character scalar or `NULL`; the object's `ci_type` slot.
#' @return A length-1 character, `"simultaneous"` or `"pointwise"`.
#' @keywords internal
#' @noRd
.band_label <- function(ci_type) {
	if (identical(ci_type, "simultaneous")) {
		"simultaneous"
	} else {
		"pointwise"
	}
}

#' @title Rendered caveat for a high-dimensional post-selection fallback band
#'
#' @description
#' Returns the two-line notice `print()` and `summary()` render under the CATT
#' preview when the fit's stored band is the `p >= NT` post-selection fallback,
#' or `NULL` when it is not. The interactive user reads the band on screen and
#' never sees the `warning()` [simultaneousCIs()] and [eventStudy()] emit for
#' the same band (#433); the two reach disjoint users and neither substitutes
#' for the other.
#'
#' Takes **scalars, never the object**, for the reason `.band_label()` states
#' just above: the print path holds a `fetwfe`-family fit and the summary path a
#' `summary.fetwfe`-family list, which carries no top-level `p` / `N` / `T` and
#' is not `inherits(., "fetwfe")`. An object parameter here aborts every
#' `print(summary(fit))` of every class at every dimension -- `x$p >= x$N * x$T`
#' is `logical(0)` on a summary and `&&` raises on it. So the summary side
#' resolves these values while it still holds the fit
#' (`.summary_estimator_output()`) and stores the finished string; its print
#' method renders whatever is stored and gates on nothing.
#'
#' @details
#' **A fit-time band is always the analytic one.** `.fit_band_for_family()`
#' passes no `method`, so it takes the `"analytic"` default -- which means a
#' `p >= NT` fit carrying `ci_type = "simultaneous"` has a post-selection
#' fallback band stored in `catt_df` whatever its class, and both `fetwfe` and
#' non-`fetwfe` fits need the notice (with different remedies, from the shared
#' `.highdim_band_remedy()`).
#'
#' **`isTRUE(calc_ses)` is load-bearing, not defensive dressing.** A plain
#' `fetwfe(gls = FALSE)` fit can be `p >= NT` with `ci_type = "simultaneous"`
#' and `catt_df` bounds all `NA`: the fit-time band came back `NULL` and
#' `.finalize_ci_type()` left `ci_type` alone. Without this conjunct the notice
#' asserts a full sentence about the properties of a band that does not exist.
#' It is necessary but **not** sufficient -- a singular selected-support Gram
#' can leave a `calc_ses = TRUE` object with no applied band, and the object
#' carries no positive signal that the band was applied (a width comparison does
#' not work: under the `K = 1` bypass the applied band's width equals the
#' pointwise width). That residual is the pre-existing
#' `ci_type`-over-`NA`-bounds mislabelling, filed separately.
#'
#' **Degenerate fits render nothing, deliberately.** A `p >= NT` fit whose
#' bridge zeroed every treatment effect carries the all-zero degenerate band,
#' not a post-selection fallback, and its own #304 condition already covers it
#' -- the same carve-out the warning makes. In practice such a fit reaches this
#' helper only if it also satisfies the gate; the omission is recorded here so
#' it reads as a decision rather than an oversight.
#'
#' **No precondition runs before these reads.** The print / summary /
#' print.summary family calls none (#447; `.workflow/PROFILE.md` section 9), so
#' the three guards below are the only thing between a malformed object and a
#' rendered result. Each is written as its own `if`, wrapping the numeric
#' comparison in `isTRUE()`, so that no conjunct ordering is load-bearing and a
#' `NULL` or zero-length slot returns `NULL` rather than raising.
#'
#' @param ci_type Character scalar or `NULL`; the object's `ci_type` slot.
#' @param p,N,T_ Numeric scalars; the design's column count, unit count, and
#'   time-period count. `T_` is spelled with the underscore because `T` is
#'   `TRUE`'s alias.
#' @param calc_ses Logical scalar; whether the fit computed analytic standard
#'   errors.
#' @param is_fetwfe Logical scalar; `inherits(x, "fetwfe")` **on the fit**, not
#'   on a summary object (a `summary.fetwfe` does not inherit `fetwfe`).
#' @return A length-1 character ready to `cat()` -- two lines plus a trailing
#'   blank line -- or `NULL` when no notice applies.
#' @keywords internal
#' @noRd
.highdim_postselection_band_notice <- function(
	ci_type,
	p,
	N,
	T_,
	calc_ses,
	is_fetwfe
) {
	if (!identical(ci_type, "simultaneous")) {
		return(NULL)
	}
	if (!isTRUE(calc_ses)) {
		return(NULL)
	}
	if (!isTRUE(p >= N * T_)) {
		return(NULL)
	}
	paste0(
		"Note: p >= N*T and this band is the post-selection fallback, which ",
		"under-covers.\n",
		"  ",
		.highdim_band_remedy(is_fetwfe),
		"\n\n"
	)
}

#' @title Two-sided Wald interval for a scalar estimate
#'
#' @description
#' Returns `estimate +/- qnorm(1 - alpha / 2) * se` as an unnamed length-2
#' numeric, lower bound first. Value-returning rather than `cat()`-emitting
#' because one of its three callers, `.tidy_estimator_output()` in
#' `R/broom_methods.R`, writes the bounds into a data frame and never prints.
#'
#' @details
#' **The `alpha` parameter carries two different notions of "the level."** The
#' two rendering callers pass the fit-time `x$alpha`, the level the fit's own
#' bands were built at; the `broom` caller passes `1 - conf.level` from the
#' user's argument to `tidy()`. Both are correct for their site and the
#' arithmetic is the same, but the join is worth naming.
#'
#' **The callers also disagree about what has been validated first.**
#' `.tidy_estimator_output()` calls `.check_for_tidy(x)` before reaching here;
#' the two rendering callers call no method-entry precondition, because the
#' print / summary family has none (see `.workflow/PROFILE.md`). In particular
#' nothing establishes `has_valid_ses`, so on a `q >= 1` fit `att_se` is `NA`,
#' this helper propagates it, and the rendered interval is `[NA, NA]`. That is
#' the status quo rather than a regression, but a reader of a shared helper is
#' entitled to know its precondition differs by caller.
#'
#' `unname()` makes the documented return type true regardless of the caller:
#' at the summary site the inputs are the *named* elements `x$att["estimate"]`
#' and `x$att["se"]`, so without it this helper would return a named pair there
#' and an unnamed one at the other two sites.
#'
#' @param estimate Numeric scalar; the point estimate.
#' @param se Numeric scalar; its standard error. May be `NA`.
#' @param alpha Numeric scalar in (0, 1); the two-sided level.
#' @return An unnamed numeric of length 2: `c(lower, upper)`.
#' @keywords internal
#' @noRd
.att_wald_ci <- function(estimate, se, alpha) {
	z <- stats::qnorm(1 - alpha / 2)
	unname(c(estimate - z * se, estimate + z * se))
}

#' @title Render one preview block: header, table, truncation footer, blank line
#'
#' @description
#' The rendering skeleton shared by **four** sites: the CATT preview and the
#' event-study preview, in each of `.print_estimator_output()` and
#' `.print_summary_estimator_output()`. Its real value is that the
#' `truncated` / `n_discarded` attribute protocol -- the contract between
#' `.truncate_catt()` / `.truncate_event_study()` and the renderers -- had four
#' independent readers and now has one.
#'
#' @details
#' `header` arrives fully formatted, so this helper does **not** single-source
#' the header text; `.band_label()` covers the part of it that can drift
#' semantically.
#'
#' `more_fmt` has **four distinct correct values** across the four sites,
#' varying on two axes at once -- `and` vs `+` (print vs summary) and `cohorts`
#' vs `event times` (CATT vs event study):
#'
#' \preformatted{
#'   print   / CATT         "  ... and %d more cohorts.\n"
#'   print   / event study  "  ... and %d more event times.\n"
#'   summary / CATT         "  ... + %d more cohorts.\n"
#'   summary / event study  "  ... + %d more event times.\n"
#' }
#'
#' Getting one wrong renders plausible-looking output at the wrong site, which
#' is why all four are pinned in
#' `tests/testthat/test-print-summary-single-source-439.R`.
#'
#' The `n_discarded` guard is not decorative: `sprintf("%d", NULL)` returns
#' `character(0)` and `cat(character(0))` emits nothing, so a frame flagged
#' `truncated` with no count would silently emit no footer at all -- in the
#' helper whose whole purpose is that footer.
#'
#' Reachability of that error, so it is not re-derived: both producers
#' (`.truncate_catt()`, `.truncate_event_study()`) set `truncated` and
#' `n_discarded` together, and every call site but one routes its frame through
#' one of them, which overwrites any stale `truncated` attribute. The exception
#' is `summary(object, full_catt = TRUE)`, which passes `object$catt_df`
#' straight through un-truncated -- so a `catt_df` carrying a hand-set
#' `truncated = TRUE` and no count would reach here. That still requires a
#' mutated fit object; it is not reachable from a well-formed one.
#'
#' @param df A data frame to render, optionally carrying the attributes
#'   `truncated` (logical) and `n_discarded` (integer).
#' @param header Single pre-formatted string, passed verbatim to `cat()`.
#' @param more_fmt A `sprintf` format taking one `%d`, used only when `df` is
#'   flagged truncated.
#' @return `invisible(NULL)`.
#' @keywords internal
#' @noRd
.cat_preview_block <- function(df, header, more_fmt) {
	cat(header)
	.print_catt_tbl(df)
	if (isTRUE(attr(df, "truncated"))) {
		n_discarded <- attr(df, "n_discarded")
		if (length(n_discarded) != 1L) {
			stop(
				".cat_preview_block(): `truncated` is TRUE but `n_discarded` ",
				"is missing or not a scalar. This is a programmer-side ",
				"contract violation -- please report at ",
				"https://github.com/gregfaletto/fetwfePackage/issues.",
				call. = FALSE
			)
		}
		cat(sprintf(more_fmt, n_discarded))
	}
	cat("\n")
	invisible(NULL)
}

#' @title Render the Model Details block
#'
#' @description
#' The header plus five dimension rows, and -- when `show_lambda` is `TRUE` --
#' the two bridge-selection rows. Shared by `.print_estimator_output()` and
#' `.print_summary_estimator_output()`, which rendered byte-identical text from
#' differently-named sources before #439.
#'
#' @details
#' **`info` is a named list, never positional scalars.** Every dimension row is
#' `%d` applied to a same-typed integer, so a positional `(N, T, G, d, p)`
#' signature would accept any permutation of its arguments silently -- and every
#' snapshot fixture predating #439 had `G == 2` and `d == 2`, so a G/d
#' transposition would have rendered identically and been invisible to the whole
#' suite. Passing names makes the transposition impossible to write.
#'
#' The summary caller passes `object$model_info` unmodified; the extra fields it
#' carries (`R`, `sig_eps_sq`, `sig_eps_c_sq`) are ignored. The print caller
#' builds the list from the fit's own slots, where the one vocabulary difference
#' lives: the fit's slot is `lambda_star_model_size`, the field here is
#' `model_size`.
#'
#' The required-field check exists because `sprintf("%d", NULL)` returns
#' `character(0)` and `cat(character(0))` emits nothing, so a mis-named *or*
#' `NULL`-valued field would print this block silently missing a whole row. It
#' tests `length(info[[nm]]) == 1L` rather than name presence, because a
#' names-only test would be inert at the `print()` call site: that caller builds
#' the list from fit slots, so a renamed slot arrives as a present name bound to
#' `NULL`. The required set is
#' derived from `show_lambda` so the two optional fields are demanded exactly
#' when they are read. Reads use `[[` rather than `$` for exact matching -- `$`
#' partial-matches on lists, and this helper's contract tolerates extra names,
#' which is the situation where a partial match can appear later. The two are
#' independent improvements: `[[` on a list returns `NULL` for a missing name
#' and never errors.
#'
#' Note the deliberate asymmetry with `.band_label()`, two helpers above, which
#' *tolerates* a missing `ci_type` on a pre-1.16.0 object. A missing `ci_type`
#' has a correct fallback -- those bounds really are pointwise -- whereas a
#' missing `N` has none: there is no value to print and no way to guess one.
#' Reachable in practice only from a `summary.<class>` object deserialized from
#' a build predating the #222 `R` -> `G` rename, since `model_info` is built in
#' exactly one place from an already-validated fit.
#'
#' @param info Named list with `N`, `T`, `G`, `d`, `p`, plus `model_size` and
#'   `lambda_star` when `show_lambda` is `TRUE`. Extra names are ignored.
#' @param show_lambda Logical; render the `Selected size` / `Lambda*` rows.
#' @return `invisible(NULL)`.
#' @keywords internal
#' @noRd
.cat_model_details <- function(info, show_lambda) {
	required <- c("N", "T", "G", "d", "p")
	if (show_lambda) {
		required <- c(required, "model_size", "lambda_star")
	}
	# Checks LENGTH, not just presence of the name. A names-only check
	# (`setdiff(required, names(info))`) is inert at the `print()` call site:
	# that caller builds `info` as `list(N = x$N, ...)`, so a renamed or dropped
	# slot arrives as a present name bound to `NULL` -- and `list(N = NULL)`
	# retains the name. `setdiff` then returns nothing, `sprintf("%d", NULL)`
	# returns `character(0)`, `cat(character(0))` emits nothing, and the row
	# vanishes silently. `length(info[[nm]]) == 1L` catches the absent name and
	# the NULL-valued name alike, so both callers are protected.
	is_usable <- vapply(
		required,
		function(nm) length(info[[nm]]) == 1L,
		logical(1)
	)
	if (!all(is_usable)) {
		stop(
			".cat_model_details(): missing or non-scalar field(s): ",
			paste(required[!is_usable], collapse = ", "),
			". This is a programmer-side contract violation -- please report ",
			"at https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
	cat("Model Details:\n")
	cat(sprintf("  Units (N)           : %d\n", info[["N"]]))
	cat(sprintf("  Time periods (T)    : %d\n", info[["T"]]))
	cat(sprintf("  Treated cohorts (G) : %d\n", info[["G"]]))
	cat(sprintf("  Covariates (d)      : %d\n", info[["d"]]))
	cat(sprintf("  Features (p)        : %d\n", info[["p"]]))
	if (show_lambda) {
		cat(sprintf("  Selected size       : %d\n", info[["model_size"]]))
		cat(sprintf("  Lambda*             : %.4f\n", info[["lambda_star"]]))
	}
	invisible(NULL)
}

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
	max_event_times,
	include_event_study = TRUE,
	...
) {
	cat(header)

	## Overall ATT
	ci_pct <- 100 * (1 - x$alpha)
	att_ci <- .att_wald_ci(x$att_hat, x$att_se, x$alpha)
	ci_low <- att_ci[1]
	ci_high <- att_ci[2]
	cat(sprintf(
		"Overall Average Treatment Effect (ATT):\n  Estimate:   %.4f\n",
		x$att_hat
	))
	cat(sprintf(
		"  Std. Error%s: %.4f\n",
		.se_qualifier(x$se_type),
		x$att_se
	))
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
	# Band-type label (#197): the reported CI bounds are simultaneous
	# (family-wise) by default, or pointwise when the fit used
	# ci_type = "pointwise". Older fits (pre-1.16.0) carry no ci_type slot
	# and are labeled pointwise (their bounds are pointwise).
	band_label <- .band_label(x$ci_type)
	catt_df <- .truncate_catt(x$catt_df, max_cohorts, order_by)
	.cat_preview_block(
		catt_df,
		sprintf(
			"Cohort Average Treatment Effects (CATT) [%s %.0f%% CI]:\n",
			band_label,
			ci_pct
		),
		"  ... and %d more cohorts.\n"
	)

	## High-dimensional post-selection fallback caveat (#433). `x` is the fit,
	## so every ingredient is read straight off it. Rendered here, under the
	## CATT preview, rather than once per preview: the two previews display the
	## same band under the same `[<label> NN% CI]` header, and one caveat per
	## object is the point.
	highdim_notice <- .highdim_postselection_band_notice(
		ci_type = x$ci_type,
		p = x$p,
		N = x$N,
		T_ = x$T,
		calc_ses = x$calc_ses,
		is_fetwfe = inherits(x, "fetwfe")
	)
	if (!is.null(highdim_notice)) {
		cat(highdim_notice)
	}

	## Event-study effects (#174 / #138). Strict policy: `eventStudy()`'s
	## contract is "succeeds on any fit produced by `fetwfe()` /
	## `betwfe()` / `etwfe()`". If it fails on a valid fit, that is
	## always a bug — let the error propagate so it surfaces to the user
	## and the developer rather than producing a print/summary that
	## silently omits the section. The prior tryCatch swallow (added
	## defensively for hypothetical future twfeCovs reuse, #138) masked
	## the scattered-cohort bug on `bacondecomp::divorce` for two
	## releases; #174 is the cure.
	# `include_event_study = FALSE` (twfeCovs, #58) skips the event-study
	# section: twfeCovs estimates one pooled effect per cohort, so there is no
	# per-(cohort,time) / event-study basis and `eventStudy()` rejects it.
	es <- if (isTRUE(include_event_study)) eventStudy(x) else NULL
	if (!is.null(es) && nrow(es) > 0L) {
		es_preview <- .truncate_event_study(es, max_event_times)
		.cat_preview_block(
			es_preview,
			sprintf(
				"Event-Study Average Treatment Effects (per event time) [%s %.0f%% CI]:\n",
				band_label,
				ci_pct
			),
			"  ... and %d more event times.\n"
		)
	}

	## Model info
	# Named fields, not positional scalars: every row is `%d` on a same-typed
	# integer, so a G/d transposition would render identically. `model_size` is
	# this caller's name for the fit's `lambda_star_model_size` slot.
	.cat_model_details(
		list(
			N = x$N,
			T = x$T,
			G = x$G,
			d = x$d,
			p = x$p,
			model_size = x$lambda_star_model_size,
			lambda_star = x$lambda_star
		),
		show_lambda
	)

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
	full_catt,
	max_event_times = 20L,
	include_event_study = TRUE
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
	# Compute event-study on demand (#174 / #138). Strict policy: see the
	# matching block in `.print_estimator_output()` for rationale. Any
	# `eventStudy()` failure on a valid fit is a bug and surfaces here
	# rather than being swallowed.
	# `include_event_study = FALSE` (twfeCovs, #58): skip the event-study
	# preview (twfeCovs has no per-(cohort,time) basis; see .print_estimator_output).
	es <- if (isTRUE(include_event_study)) eventStudy(object) else NULL
	event_study <- if (is.null(es) || nrow(es) == 0L) {
		NULL
	} else {
		.truncate_event_study(es, max_event_times = max_event_times)
	}

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
		event_study = event_study,
		model_info = list(
			N = object$N,
			T = object$T,
			G = object$G,
			R = object$G,
			d = object$d,
			p = object$p,
			lambda_star = object$lambda_star,
			model_size = object$lambda_star_model_size,
			sig_eps_sq = object$sig_eps_sq,
			sig_eps_c_sq = object$sig_eps_c_sq
		),
		alpha = object$alpha,
		se_type = object$se_type,
		# #197: carry ci_type so print.summary.<class> can label the CATT /
		# event-study previews (simultaneous vs pointwise). Pre-1.16.0 fits
		# have no slot -> NULL -> labeled pointwise downstream.
		ci_type = object$ci_type,
		# #433: the high-dimensional post-selection-fallback caveat, resolved
		# HERE because this is the last point on the summary path that holds the
		# fit. `print.summary.<class>` receives a `summary.<class>` list with no
		# top-level p / N / T and no `calc_ses` at all, so it could neither gate
		# nor build this. A length-1 character, or NULL when no notice applies.
		highdim_band_notice = .highdim_postselection_band_notice(
			ci_type = object$ci_type,
			p = object$p,
			N = object$N,
			T_ = object$T,
			calc_ses = object$calc_ses,
			is_fetwfe = inherits(object, "fetwfe")
		)
	)

	keep <- c(
		"att",
		if (include_att_selected) "att_selected",
		"catt",
		# `event_study` field is always present in the list slot when
		# `eventStudy()` succeeded (NULL otherwise). Placed between
		# `catt` and `model_info` to mirror the print-order convention.
		"event_study",
		"model_info",
		"alpha",
		"se_type",
		# #197: ci_type appears last, after se_type (the band-type label
		# source for the CATT / event-study previews).
		"ci_type",
		# #433: grouped with `ci_type`, the other band-related field. This
		# vector is a WHITELIST -- `out <- out[keep]` below silently drops any
		# name not listed here, with no error -- so a field added to `out` and
		# not to `keep` never reaches the summary object and never renders.
		"highdim_band_notice"
	)
	out <- out[keep]

	out$model_info <- out$model_info[c(
		"N",
		"T",
		"G",
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
	# Band-type label (#197): simultaneous (family-wise) by default, pointwise
	# under ci_type = "pointwise" or for pre-1.16.0 summaries with no slot.
	band_label <- .band_label(x$ci_type)
	att_ci <- .att_wald_ci(x$att["estimate"], x$att["se"], x$alpha)
	ci_low <- att_ci[1]
	ci_high <- att_ci[2]
	p_val <- x$att["p_value"]
	p_str <- if (is.na(p_val)) "NA" else sprintf("%.4g", p_val)
	se_label <- paste0("SE", .se_qualifier(x$se_type))
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

	.cat_preview_block(
		x$catt,
		sprintf(
			"CATT (preview) [%s %.0f%% CI]:\n",
			band_label,
			ci_pct
		),
		"  ... + %d more cohorts.\n"
	)

	## High-dimensional post-selection fallback caveat (#433). Renders the
	## string `.summary_estimator_output()` resolved while it still held the
	## fit; this function performs NO gating, because it does not hold the
	## ingredients (no top-level p / N / T, no `calc_ses`) and must not pretend
	## to. The `is.null()` guard is required, not stylistic: `cat(NULL, "\n\n")`
	## emits two blank lines, which would appear in every `summary()` of every
	## class at every dimension and drift both snapshot goldens. It also keeps a
	## pre-#433 serialized summary object (no such field) rendering correctly.
	if (!is.null(x$highdim_band_notice)) {
		cat(x$highdim_band_notice)
	}

	## Event-study preview (#174 / #138). Reads from the cached field set
	## by `.summary_estimator_output()`; no recompute. NULL means
	## `eventStudy()` returned an empty data frame (e.g., G = 0); a true
	## eventStudy() failure now propagates rather than being silently
	## swallowed (strict policy, #174).
	if (!is.null(x$event_study)) {
		.cat_preview_block(
			x$event_study,
			sprintf(
				"Event Study (preview) [%s %.0f%% CI]:\n",
				band_label,
				ci_pct
			),
			"  ... + %d more event times.\n"
		)
	}

	## Model info
	# `model_info` already carries exactly the names .cat_model_details() wants,
	# so it is passed with no adapter; its extra fields (R, sig_eps_sq,
	# sig_eps_c_sq) are ignored.
	.cat_model_details(x$model_info, show_lambda)

	invisible(x)
}
