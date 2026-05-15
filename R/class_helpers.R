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
		cohort = order(catt_df$Cohort),
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
