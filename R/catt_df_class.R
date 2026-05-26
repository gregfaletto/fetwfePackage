#-------------------------------------------------------------------------------
# S3 helpful-error layer for the `catt_df` class (issue #136).
#
# As of fetwfe 1.11.0 the columns of `result$catt_df` were renamed to a
# snake_case convention matching `eventStudy(result)`:
#
#   Cohort       -> cohort
#   Estimated TE -> estimate
#   SE           -> se
#   ConfIntLow   -> ci_low
#   ConfIntHigh  -> ci_high
#   P_value      -> p_value
#
# The `selected` column (fetwfe / betwfe only) is unchanged.
#
# To make migration of pre-1.11.0 user scripts loud and actionable rather
# than silently returning NULL, every `catt_df` carries class
# c("catt_df", "data.frame"). The three S3 methods below intercept access
# by the old column names via `[[`, `$`, and `[` and stop() with a message
# pointing to the new name. Any other column access falls through to the
# data.frame method via NextMethod().
#
# NSE access (dplyr, ggplot, subset) strips the catt_df class as usual
# and surfaces the standard "object not found" error -- the NEWS rename
# table compensates for that case. The interceptors only catch the
# `[[` / `$` / `[` triplet on the raw frame, which is where the bulk of
# migration friction lives.
#-------------------------------------------------------------------------------

#' @title Old -> new column-name mapping for `catt_df` (internal)
#' @description
#' Single source of truth for the 1.10.0 -> 1.11.0 rename. Consulted by
#' the three S3 helpful-error methods below and asserted on by a
#' regression-guard test in `tests/testthat/test-catt_df-helpful-errors.R`.
#' @keywords internal
#' @noRd
.catt_df_rename_map <- c(
	"Cohort" = "cohort",
	"Estimated TE" = "estimate",
	"SE" = "se",
	"ConfIntLow" = "ci_low",
	"ConfIntHigh" = "ci_high",
	"P_value" = "p_value"
)

#' @title Build the standard migration-message stop() for an old name
#' @keywords internal
#' @noRd
.catt_df_renamed_error <- function(name) {
	new_name <- .catt_df_rename_map[[name]]
	stop(
		sprintf(
			"Column `%s` was renamed to `%s` in fetwfe 1.11.0. See NEWS.md.",
			name,
			new_name
		),
		call. = FALSE
	)
}

#' Double-bracket access on a `catt_df` object
#'
#' Intercepts access by the pre-1.11.0 Title-Case column names
#' (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`,
#' `P_value`) and stops with a migration message pointing to the new
#' snake_case name. All other access falls through to the
#' `data.frame` method via `NextMethod()`.
#'
#' @param x A `catt_df` object (data frame with class
#'   `c("catt_df", "data.frame")`).
#' @param i Index; passed through to `[[.data.frame`. Character
#'   indices matching an old column name are intercepted; other
#'   indices fall through.
#' @param ... Further arguments passed through.
#' @return The column or value, as for `[[.data.frame`.
#' @method [[ catt_df
#' @export
`[[.catt_df` <- function(x, i, ...) {
	if (
		is.character(i) && length(i) == 1L && i %in% names(.catt_df_rename_map)
	) {
		.catt_df_renamed_error(i)
	}
	NextMethod()
}

#' Dollar-sign access on a `catt_df` object
#'
#' Intercepts access by the pre-1.11.0 Title-Case column names and
#' stops with a migration message. All other access falls through
#' to the `data.frame` method via `NextMethod()`.
#'
#' @param x A `catt_df` object.
#' @param name Character; the column name being accessed via `$`.
#' @return The column, as for `$.data.frame`.
#' @method $ catt_df
#' @export
`$.catt_df` <- function(x, name) {
	if (
		is.character(name) &&
			length(name) == 1L &&
			name %in% names(.catt_df_rename_map)
	) {
		.catt_df_renamed_error(name)
	}
	NextMethod()
}

#' Single-bracket access on a `catt_df` object
#'
#' Intercepts column selection by the pre-1.11.0 Title-Case names and
#' stops with a migration message. The check fires when an old name
#' appears in the column-selector position: `j` for the `df[i, j]`
#' two-index form, or `i` for the `df[j]` one-index column-selection
#' form. Row-only access (`df[i, ]`) and access by new column names
#' fall through to the `data.frame` method via `NextMethod()`.
#'
#' @param x A `catt_df` object.
#' @param i,j Row / column selectors; see `[.data.frame`.
#' @param ... Further arguments passed through.
#' @return The selected subset, as for `[.data.frame`.
#' @method [ catt_df
#' @export
`[.catt_df` <- function(x, i, j, ...) {
	# Distinguish x[i] (single-index column selection) from x[i, ]
	# (single-index row selection with explicit trailing comma) and
	# x[i, j] (two-index form). The call form is the only reliable
	# signal: nargs() counts positional arguments, where x[i, ] has
	# nargs() == 3 (j is positionally present but empty) and x[i]
	# has nargs() == 2. `drop` and other named args are not counted
	# by nargs() because nargs() ignores named arguments.
	col_arg <- NULL
	if (nargs() <= 2L) {
		# x[i] form: i is the column selector (possibly missing for x[]).
		if (!missing(i)) {
			col_arg <- i
		}
	} else {
		# x[i, j] or x[, j] form: j is the column selector.
		if (!missing(j)) {
			col_arg <- j
		}
	}
	if (
		!is.null(col_arg) &&
			is.character(col_arg) &&
			any(col_arg %in% names(.catt_df_rename_map))
	) {
		hit <- intersect(col_arg, names(.catt_df_rename_map))[1L]
		.catt_df_renamed_error(hit)
	}
	NextMethod()
}
