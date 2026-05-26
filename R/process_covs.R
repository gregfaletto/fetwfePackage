# Formula-interface support for the `covs` argument (#28).
#
# Normalizes `covs` input to a character vector of column names. Accepts
# either the long-standing character-vector form (`covs = c("x1", "x2")`)
# or a one-sided formula (`covs = ~ x1 + x2`), mirroring
# `did::att_gt(xformla = ~ x1 + x2)`. Called at the top of each
# estimator's body (`fetwfe()` / `etwfe()` / `betwfe()` / `twfeCovs()`)
# so the rest of the package can assume a character vector.

#' @title Normalize the `covs` argument to a character vector
#'
#' @description
#' Internal helper used by the four estimator entry points. Accepts a
#' character vector (returned unchanged) or a one-sided formula
#' (e.g., `~ x1 + x2`), extracting the term labels. Anything else
#' triggers a structured stop. Only additive bare names are supported;
#' interactions, transformations, and `I()` wrappers are rejected so
#' the user gets a clear pointer to compute the derived variable in
#' the data frame first.
#'
#' @param covs Either a character vector of column names (the current
#'   form, including the empty `c()` default), or a one-sided formula.
#'   `NULL` is treated as no covariates (returns `character(0)`).
#' @return A character vector of column names. Length zero when `covs`
#'   is `NULL`, `c()`, or `~ 1`.
#' @keywords internal
#' @noRd
.process_covs_input <- function(covs) {
	if (is.null(covs)) {
		return(character(0))
	}
	if (is.character(covs)) {
		return(covs)
	}
	if (inherits(covs, "formula")) {
		# One-sided formula required: `~ x1 + x2` -- `length(f) == 2`,
		# with `f[[1]]` the `~` operator and `f[[2]]` the RHS.
		if (length(covs) != 2L) {
			stop(
				"`covs` formula must be one-sided (RHS only), e.g. `~ x1 + x2`. ",
				"Got: ",
				deparse(covs),
				".",
				call. = FALSE
			)
		}
		tt <- stats::terms(covs)
		term_labels <- attr(tt, "term.labels")
		# Length-zero (e.g., `~ 1` or `~ 0`) is a valid "no covariates"
		# spelling -- treat like the empty character vector.
		if (length(term_labels) == 0L) {
			return(character(0))
		}
		# Only bare variable names are supported (no interactions,
		# I()-transforms, etc.). The user can compute derived variables
		# in the data frame and pass them by name instead.
		invalid <- grepl("[^A-Za-z0-9._]", term_labels)
		if (any(invalid)) {
			stop(
				"`covs` formula only supports additive bare variable names ",
				"(e.g. `~ x1 + x2`). Got non-bare term(s): ",
				paste(term_labels[invalid], collapse = ", "),
				". For derived variables, compute them in the data frame ",
				"first and pass via the character-vector form ",
				"(`covs = c(\"x1\", \"derived_x2\")`).",
				call. = FALSE
			)
		}
		return(term_labels)
	}
	stop(
		"`covs` must be a character vector of column names or a one-sided ",
		"formula (e.g. `~ x1 + x2`). Got object of class: ",
		paste(class(covs), collapse = ", "),
		".",
		call. = FALSE
	)
}
