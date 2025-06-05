#' Convert data formatted for `att_gt()` to a dataframe suitable for `fetwfe()` / `etwfe()`
#'
#' `attgtToFetwfeDf()` reshapes and renames a panel dataset that is already
#' formatted for `did::att_gt()` (Callaway and Sant'Anna 2021) so that it can be
#' passed directly to fetwfe()` or `etwfe()` from the `fetwfe` package. In
#' particular, it
#'   * creates an *absorbing‑state* treatment dummy that equals 1 from the
#'     first treated period onward* and 0 otherwise,
#'   * (optionally) drops units that are already treated in the very first
#'     period of the sample (because `fetwfe()` removes them internally), and
#'   * returns a tidy dataframe whose column names match the arguments that
#'     `fetwfe()`/`etwfe()` expect.
#'
#' @param data A `data.frame` in **long** format containing at least the four
#'   columns used by `did::att_gt()`: outcome `yname`, time `tname`, unit id
#'   `idname`, and the first‑treatment period `gname` (which is 0 for the
#'   never‑treated group).
#' @param yname  Character scalar. Name of the outcome column.
#' @param tname  Character scalar. Name of the time variable (numeric or
#'   integer). This becomes `time` in the returned dataframe.
#' @param idname Character scalar. Name of the unit identifier. Converted to
#'   character and returned as `unit_var`.
#' @param gname  Character scalar. Name of the *group* variable holding the
#'   first period of treatment. Values must be 0 for never‑treated, or a
#'   positive integer representing the first treated period.
#' @param covars Character vector of additional covariate column names to carry
#'   through (default `character(0)`). These columns are left untouched and
#'   appear *after* the required columns in the returned dataframe.
#' @param drop_first_period_treated Logical. If `TRUE` (default), units that
#'   are already treated in the first sample period are removed *before*
#'   creating the treatment dummy. `fetwfe()` would do this internally, but
#'   dropping them here keeps the returned dataframe cleaner.
#' @param out_names  A named list giving the column names to use in the
#'   resulting dataframe. Defaults are `list(time = "time", unit = "unit",
#'   treatment = "treatment", response = "y")`. Override if you prefer
#'   different names (for instance, to keep the original `yname`). The vector
#'   *must* contain exactly these four names.
#'
#' @return A `data.frame` with columns `time`, `unit`, `treatment`, `y`, and any
#'   covariates requested in `covars`, ready to be fed to
#'   `fetwfe()`/`etwfe()`. All required columns are of the correct type:
#'   `time` is integer, `unit` is character, `treatment` is integer 0/1, and
#'   `y` is numeric.
#' @references Callaway, Brantly and Pedro H.C. Sant'Anna. "Difference-in-
#' Differences with Multiple Time Periods." Journal of Econometrics, Vol. 225,
#' No. 2, pp. 200-230, 2021.
#' \url{https://doi.org/10.1016/j.jeconom.2020.12.001},
#' \url{https://arxiv.org/abs/1803.09015}.
#' @examples
#' ## toy example ---------------------------------------------------------------
#' \dontrun{
#' library(did)  # provides the mpdta example dataframe
#' data(mpdta)
#'
#' head(mpdta)
#'
#' tidy_df <- attgtToFetwfeDf(
#'   data  = mpdta,
#'   yname = "lemp",
#'   tname = "year",
#'   idname = "countyreal",
#'   gname = "first.treat",
#'   covars = c("lpop"))
#'
#' head(tidy_df)
#' }
#'
#' ## Now you can call fetwfe()  ------------------------------------------------
#' # res <- fetwfe(
#' #   pdata      = tidy_df,
#' #   time_var   = "time_var",
#' #   unit_var   = "unit_var",
#' #   treatment  = "treatment",
#' #   response   = "response",
#' #   covs       = c("lpop"))
#'
#' @export
attgtToFetwfeDf <- function(
	data,
	yname,
	tname,
	idname,
	gname,
	covars = character(0),
	drop_first_period_treated = TRUE,
	out_names = list(
		time = "time_var",
		unit = "unit_var",
		treatment = "treatment",
		response = "response"
	)
) {
	.fetwfe_df_core(
		data = data,
		vars = list(y = yname, t = tname, id = idname, g = gname),
		covars = covars,
		drop_first_period_treated = drop_first_period_treated,
		out_names = out_names
	)
}

#' Convert data prepared for `etwfe::etwfe()` to the format required by
#' `fetwfe()` and `fetwfe::etwfe()`
#'
#' `etwfeToFetwfeDf()` reshapes and renames a panel dataset that is already
#' formatted for `etwfe::etwfe()` (McDermott 2024) so that it can be
#' passed directly to fetwfe()` or `etwfe()` from the `fetwfe` package. In
#' particular, it
#'   * creates an *absorbing‑state* treatment dummy that equals 1 from the
#'     first treated period onward* and 0 otherwise,
#'   * (optionally) drops units that are already treated in the very first
#'     period of the sample (because `fetwfe()` removes them internally), and
#'   * returns a tidy dataframe whose column names match the arguments that
#'     `fetwfe()`/`etwfe()` expect.
#'
#' @param data A long-format data.frame that you could already feed to `etwfe()`.
#' @param yvar Character. Column name of the outcome (left-hand side in your `fml`).
#' @param tvar Character. Column name of the time variable that you pass to `etwfe()` as `tvar`.
#' @param idvar Character. Column name of the unit identifier (the variable you would
#'   cluster on, or pass to `etwfe(..., ivar = idvar)` if you were using unit FEs).
#' @param gvar Character. Column name of the “first treated” cohort variable passed to `etwfe()` as `gvar`.
#'   Must be `0` for never-treated units, or the (strictly positive) first treated period.
#' @param covars Character vector of *additional* covariate columns to keep (default `character(0)`).
#' @param drop_first_period_treated Logical. Should units already treated in the very first
#'   sample period be removed?  (`fetwfe()` will drop them internally anyway, but doing it
#'   here keeps the returned dataframe clean.)  Default `TRUE`.
#' @param out_names Named list giving the column names that the returned dataframe should have.
#'   The default (`time`, `unit`, `treatment`, `y`) matches the arguments usually supplied to
#'   `fetwfe()`. **Do not change the *names* of this list** – only the *values* – and keep all four.
#'
#' @return A tidy `data.frame` with (in this order)
#'   * `time`       integer,
#'   * `unit`       character,
#'   * `treatment`  integer 0/1 absorbing-state dummy,
#'   * `response`   numeric outcome,
#'   * any covariates requested in `covars`.
#'   Ready to pass straight to `fetwfe()` or `fetwfe::etwfe()`.
#'
#' @references McDermott G (2024). _etwfe: Extended Two-Way Fixed Effects_.
#' doi:10.32614/CRAN.package.etwfe
#' \url{https://doi.org/10.32614/CRAN.package.etwfe}, R package
#' version 0.5.0, \url{https://CRAN.R-project.org/package=etwfe}.
#' @examples
#' ## toy example ---------------------------------------------------------------
#' \dontrun{
#' library(did)  # provides the mpdta example dataframe
#' data(mpdta)
#'
#' head(mpdta)
#'
#' tidy_df <- etwfeToFetwfeDf(
#'   data  = mpdta,
#'   yvar = "lemp",
#'   tvar = "year",
#'   idvar = "countyreal",
#'   gvar = "first.treat",
#'   covars = c("lpop"))
#'
#' head(tidy_df)
#'
#' }
#' ## Now you can call fetwfe()  ------------------------------------------------
#' # res <- fetwfe(
#' #   pdata      = tidy_df,
#' #   time_var   = "time_var",
#' #   unit_var   = "unit_var",
#' #   treatment  = "treatment",
#' #   response   = "response",
#' #   covs       = c("lpop"))
#'
#' @export
etwfeToFetwfeDf <- function(
	data,
	yvar,
	tvar,
	idvar,
	gvar,
	covars = character(0),
	drop_first_period_treated = TRUE,
	out_names = list(
		time = "time_var",
		unit = "unit_var",
		treatment = "treatment",
		response = "response"
	)
) {
	.fetwfe_df_core(
		data = data,
		vars = list(y = yvar, t = tvar, id = idvar, g = gvar),
		covars = covars,
		drop_first_period_treated = drop_first_period_treated,
		out_names = out_names
	)
}


#' Core converter for "*_ToFetwfeDf" helpers
#'
#' `.fetwfe_df_core()` holds the shared back-end logic used by
#' [etwfeToFetwfeDf()] and [attgtToFetwfeDf()].  It takes a long-format panel
#' dataset, performs a suite of validity checks and coercions, constructs an
#' *absorbing-state* treatment dummy \eqn{ 1\{ t >= g & g > 0 \}}, optionally drops
#' units treated in the very first sample period, and returns a tidy data frame
#' ready for [fetwfe::fetwfe()] / [fetwfe::etwfe()].
#'
#' @section Typical usage:
#' End-users should not call this function directly.  Instead, use one of
#' the thin wrappers that merely translate argument names:
#' * [etwfeToFetwfeDf()] - for data already labelled with `etwfe` conventions.
#' * [attgtToFetwfeDf()] - for data prepared for `did::att_gt()`.
#'
#' @param data A `data.frame` in long panel format.
#' @param vars Named list (or vector) with exactly four elements,
#'   giving the column names of the core variables in `data`:
#'   \describe{
#'     \item{`y`}{Outcome / response variable.}
#'     \item{`t`}{Time variable (numeric or integer).}
#'     \item{`id`}{Unit identifier.}
#'     \item{`g`}{"First treated" period; 0 for never-treated.}
#'   }
#' @param covars Character vector of additional covariate names to carry through
#'   unchanged (default `character(0)`).
#' @param drop_first_period_treated Logical.  If `TRUE` (default), observations
#'   that are already treated in the *earliest* sample period are removed before
#'   creating the treatment dummy.  This mirrors the internal behaviour of
#'   **fetwfe** but keeps the returned data frame cleaner.
#' @param out_names Named list of length four giving the desired column names
#'   in the returned data frame.  Defaults are
#'   `list(time = "time_var", unit = "unit_var", treatment = "treatment",
#'	 response = "response")`.
#'
#' @return A tidy `data.frame` with columns, in this order:
#' \itemize{
#'   \item **`time_var`**   - integer,
#'   \item **`unit_var`**   - character,
#'   \item **`treatment`** - integer 0 / 1 absorbing-state dummy,
#'   \item **`response`**      - numeric outcome,
#'   \item any extra covariates requested via `covars`.
#' }
#' The rows are sorted by `unit` then `time`; row names are dropped.
#'
#' @details
#' The function enforces several data-quality constraints:
#' \enumerate{
#'   \item Each `(unit, time)` pair must appear at most once.
#'   \item The first-treatment period `g` must be constant within a unit
#'     (irreversible treatment).
#'   \item `time` and `g` are coerced to `integer`; `unit` to `character`.
#' }
#' If any unit ever switches from treated back to untreated
#' (`diff(treatment) < 0`), a warning is thrown because those units will be
#' silently dropped by [fetwfe::fetwfe()].
#'
#' @keywords internal
#' @noRd
.fetwfe_df_core <- function(
	data,
	vars,
	covars = character(0),
	drop_first_period_treated = TRUE,
	out_names = list(
		time = "time_var",
		unit = "unit_var",
		treatment = "treatment",
		response = "response"
	)
) {
	stopifnot(is.data.frame(data))
	required <- unlist(vars, use.names = FALSE)
	missing <- setdiff(c(required, covars), names(data))
	if (length(missing))
		stop(
			"Column(s) not found in `data`: ",
			paste(missing, collapse = ", "),
			call. = FALSE
		)

	## ---------- enforce types --------------------------------------------------
	df <- data[, unique(c(required, covars))]
	df[[vars$t]] <- as.integer(df[[vars$t]])
	df[[vars$g]] <- as.integer(df[[vars$g]])
	df[[vars$id]] <- as.character(df[[vars$id]])

	if (anyNA(df[[vars$t]]))
		stop("Missing values in `", vars$t, "`.", call. = FALSE)
	if (anyNA(df[[vars$g]]))
		stop("Missing values in `", vars$g, "`.", call. = FALSE)

	## ---------- uniqueness & irreversibility ----------------------------------
	dup_rows <- duplicated(df[, c(vars$id, vars$t)])
	if (any(dup_rows))
		stop(
			"Each (id, time) pair must appear at most once; found duplicates.",
			call. = FALSE
		)

	g_const <- by(df, df[[vars$id]], function(x) length(unique(x[[vars$g]])))
	if (any(g_const > 1))
		stop(
			"`",
			vars$g,
			"` must be constant within each unit (irreversible treatment).",
			call. = FALSE
		)

	## ---------- drop first-period treated --------------------------------------
	first_period <- min(df[[vars$t]])
	if (drop_first_period_treated) {
		treated_first <- df[[vars$g]] != 0 & df[[vars$g]] <= first_period
		if (any(treated_first)) {
			df <- df[!treated_first, ]
			message(
				"Dropped ",
				sum(treated_first),
				" unit-period(s) treated in the first period."
			)
		}
	}

	## ---------- absorb-state dummy --------------------------------------------
	treat_dummy <- ifelse(
		df[[vars$g]] > 0 & df[[vars$t]] >= df[[vars$g]],
		1L,
		0L
	)
	bad_units <- with(
		df,
		tapply(treat_dummy, df[[vars$id]], function(z) any(diff(z) < 0))
	)
	if (any(bad_units))
		warning(
			sum(bad_units),
			" unit(s) switch from treated back to untreated. ",
			"They will be dropped by `fetwfe()`."
		)

	## ---------- assemble tidy dataframe ---------------------------------------
	res <- data.frame(
		df[[vars$t]], # time
		df[[vars$id]], # unit
		treat_dummy, # treatment dummy
		df[[vars$y]], # outcome
		stringsAsFactors = FALSE
	)
	names(res)[1:4] <- unlist(
		out_names[c("time", "unit", "treatment", "response")],
		use.names = FALSE
	)

	if (length(covars)) res[covars] <- df[covars]

	# guard against factors that slipped in
	res[[out_names$time]] <- as.integer(res[[out_names$time]])
	res[[out_names$unit]] <- as.character(res[[out_names$unit]])
	res[[out_names$treatment]] <- as.integer(res[[out_names$treatment]])
	res[[out_names$response]] <- as.numeric(res[[out_names$response]])

	res <- res[order(res[[out_names$unit]], res[[out_names$time]]), ]
	rownames(res) <- NULL
	res
}
