# Cohort-study accessor: a function-style wrapper around `result$catt_df`,
# parallel to `eventStudy()` in shape and discoverability. See issue #137
# and `.plans/cohort-study-wrapper-137/PLAN.md` for the design rationale.

#' Per-cohort average treatment effects
#'
#' @description
#' Extracts per-cohort ATT estimates from a fitted FETWFE / ETWFE / BETWFE /
#' twfeCovs object as a tidy data frame. Parallel to [eventStudy()] for
#' event-time aggregation: the per-cohort information is already available
#' via `result$catt_df`, and `cohortStudy()` surfaces it through a
#' discoverable function with its own help page (so users can reach it via
#' `?cohortStudy` without having to know the slot name).
#'
#' The function is a pass-through on `result$catt_df` modulo class: the
#' columns, their values, and their order are unchanged. The returned
#' object carries class `c("cohortStudy", "catt_df", "data.frame")`. The
#' `cohortStudy` class dispatches the broom tidier [tidy.cohortStudy()];
#' the `catt_df` class preserves the helpful-error layer (introduced in
#' fetwfe 1.11.0) that intercepts pre-1.11.0 Title-Case column names
#' (`Cohort`, `Estimated TE`, etc.) with a migration message; the
#' `data.frame` base preserves the standard data-frame methods (`print`,
#' `head`, `nrow`, `dplyr::filter`, etc.).
#'
#' @param result A fitted object from [fetwfe()], [etwfe()], [betwfe()], or
#'   [twfeCovs()] (or their `*WithSimulatedData()` wrapper analogs, which
#'   return the same classes).
#' @return A data frame with class `c("cohortStudy", "catt_df", "data.frame")`
#'   containing one row per treated cohort and columns:
#'   \describe{
#'     \item{cohort}{Character; the cohort label (the calendar time at which
#'       the cohort first received treatment).}
#'     \item{estimate}{Numeric; the per-cohort ATT estimate.}
#'     \item{se}{Numeric; standard error for the per-cohort ATT (`NA` when
#'       the Gram matrix is singular or, for `fetwfe()` / `betwfe()`, the
#'       bridge penalty zeroed out the cohort).}
#'     \item{ci_low, ci_high}{Numeric; the stored lower and upper
#'       confidence-interval bounds, reflecting the fit's `ci_type` --
#'       simultaneous (family-wise) by default, or pointwise `1 - alpha`
#'       Wald bounds when the fit used `ci_type = "pointwise"` (`alpha` is
#'       the value passed at fit time).}
#'     \item{p_value}{Numeric; follows the fit's `ci_type`. Under
#'       `"pointwise"`, the two-sided Wald p-value
#'       (`2 * pnorm(-|estimate / se|)`); under `"simultaneous"` (the
#'       default), the single-step max-T multiplicity-adjusted (family-wise)
#'       p-value matching the simultaneous band (#200). `NA` when `se` is `0`
#'       or `NA`.}
#'     \item{selected}{(`fetwfe()` / `betwfe()` only.) Logical; `TRUE` when
#'       the bridge penalty left the cohort's ATT nonzero. Absent for
#'       `etwfe()` and `twfeCovs()`, which do not perform selection.}
#'   }
#'   Use `tidy(cohortStudy(result))` (with the `broom` package loaded) to
#'   reshape to broom convention (`term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `conf.low`, `conf.high`, optionally
#'   `selected`); see [tidy.cohortStudy()].
#' @seealso [eventStudy()] for the parallel event-time accessor;
#'   [tidy.cohortStudy()] for broom-shape translation.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- fetwfeWithSimulatedData(dat)
#'   cs <- cohortStudy(res)
#'   cs
#'   # Broom-shape translation:
#'   if (requireNamespace("broom", quietly = TRUE)) {
#'     broom::tidy(cs)
#'   }
#' }
#' @export
cohortStudy <- function(result) {
	fitted_classes <- c("fetwfe", "etwfe", "betwfe", "twfeCovs")
	if (!inherits(result, fitted_classes)) {
		stop(
			"cohortStudy() requires a fitted object from fetwfe(), etwfe(), ",
			"betwfe(), or twfeCovs(). Got class: ",
			paste(class(result), collapse = ", "),
			".",
			call. = FALSE
		)
	}
	if (is.null(result$catt_df)) {
		stop(
			"cohortStudy(): the fitted object has no `catt_df` slot. This is ",
			"unexpected; please file an issue at ",
			"https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
	if (!is.data.frame(result$catt_df)) {
		stop(
			"cohortStudy(): the fitted object's `catt_df` slot is not a ",
			"data.frame (got class: ",
			paste(class(result$catt_df), collapse = ", "),
			"). If you have mutated this slot manually, restore the original ",
			"data frame; otherwise please file an issue.",
			call. = FALSE
		)
	}
	out <- result$catt_df
	class(out) <- c("cohortStudy", "catt_df", "data.frame")
	out
}
