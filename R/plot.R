# Plot methods for fetwfe / etwfe / betwfe results (#29).
#
# Three S3 `plot.<class>` methods, each a thin wrapper around the shared
# `.plot_estimator()` helper. Returns a ggplot object so users can pipe in
# additional layers / themes. ggplot2 is in Suggests (not Imports) -- the
# helper checks for it at call time and stops with an actionable message
# if missing. Mirrors the visualization style of `did::ggdid()` from the
# Callaway-Sant'Anna `did` package, which is one of the most-cited UX
# wins of that package.
#
# twfeCovs is NOT covered here -- it lacks the full class-method family
# (see #58 for the broader treatment).

# Silence R CMD check's "no visible binding for global variable" notes
# for the bare column names used inside `ggplot2::aes()`. The names are
# stable per the post-#136 snake_case convention.
utils::globalVariables(c(
	"cohort",
	"estimate",
	"se",
	"ci_low",
	"ci_high",
	"p_value",
	"selected",
	"event_time",
	"n_cohorts"
))

#' Plot CATT or event-study estimates from a fitted FETWFE / ETWFE / BETWFE
#'
#' @description
#' Returns a ggplot object showing either event-study coefficients
#' (default; `type = "event_study"`) or per-cohort average treatment
#' effects (`type = "catt"`) from a fitted estimator. Mirrors the
#' visualization style of `did::ggdid()` from the Callaway-Sant'Anna
#' `did` package, providing a single-call route from a fitted object to
#' a publication-ready visualization.
#'
#' For `fetwfe` / `betwfe` (the bridge-penalty estimators) in the CATT
#' view, points are shape- and color-coded by whether the bridge
#' penalty left that cohort's ATT nonzero (`selected = TRUE`) or
#' zeroed it out (`selected = FALSE`). For `etwfe` (no selection),
#' all points are uniformly styled.
#'
#' @param x A fitted object from [fetwfe()], [etwfe()], or [betwfe()].
#'   (`twfeCovs` is not currently supported -- see GitHub issue #58 for
#'   the broader treatment of `twfeCovs` class methods.)
#' @param type Character; either `"event_study"` (default; event-time
#'   pooled coefficients from [eventStudy()]) or `"catt"` (per-cohort
#'   ATTs from `result$catt_df`).
#' @param conf_int Logical; if `TRUE` (default), include confidence
#'   interval error bars.
#' @param alpha Numeric; overrides the fit's alpha for CI computation.
#'   `NULL` (default) uses the fit's alpha (so 95% CIs when the fit's
#'   alpha is 0.05). With the default `NULL`, both views show the fit's
#'   `ci_type` band (simultaneous by default): the `"catt"` view reads the
#'   bounds stored in `catt_df`, and the `"event_study"` view re-derives them
#'   via [eventStudy()] (which inherits the fit's `ci_type`).
#'   **Asymmetry under an explicit `alpha`:** for the `"event_study"` view, the
#'   explicit alpha is forwarded to [eventStudy()], which re-runs the joint
#'   machinery and so still produces a *simultaneous* band at that alpha (when
#'   the fit's `ci_type` is `"simultaneous"`); for the `"catt"` view, an
#'   explicit alpha recomputes `ci_low` / `ci_high` from
#'   `estimate +/- qnorm(1 - alpha/2) * se`, i.e. a *pointwise* band at that
#'   alpha (the stored simultaneous bounds are at the fit's alpha and are not
#'   re-derived at a new alpha for the catt view). To plot simultaneous catt
#'   bands at a different alpha, refit at that alpha or call
#'   [simultaneousCIs()] directly.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A `ggplot` object. Users can customize further via standard
#'   ggplot layer-addition syntax (e.g.,
#'   `plot(res) + ggplot2::theme_classic()`).
#'
#' @seealso [cohortStudy()] for the per-cohort accessor;
#'   [eventStudy()] for the event-time accessor;
#'   [plot.etwfe()], [plot.betwfe()] for the parallel methods.
#'
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- fetwfeWithSimulatedData(dat)
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot(res)                          # default: event-study coefficients
#'     plot(res, type = "catt")           # per-cohort ATTs
#'     plot(res, conf_int = FALSE)        # point estimates only
#'     plot(res, alpha = 0.1)             # 90% CIs
#'   }
#' }
#'
#' @export
plot.fetwfe <- function(
	x,
	type = c("event_study", "catt"),
	conf_int = TRUE,
	alpha = NULL,
	...
) {
	.plot_estimator(
		x,
		type = match.arg(type),
		conf_int = conf_int,
		alpha = alpha
	)
}

#' Plot CATT or event-study estimates from a fitted ETWFE
#'
#' @description
#' Parallel to [plot.fetwfe()]. ETWFE does not perform selection, so all
#' points are uniformly styled (no `selected = TRUE / FALSE` encoding).
#'
#' @inheritParams plot.fetwfe
#' @return A `ggplot` object.
#' @seealso [plot.fetwfe()] for the full documentation; [eventStudy()];
#'   [cohortStudy()].
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- etwfeWithSimulatedData(dat)
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot(res)
#'   }
#' }
#' @export
plot.etwfe <- function(
	x,
	type = c("event_study", "catt"),
	conf_int = TRUE,
	alpha = NULL,
	...
) {
	.plot_estimator(
		x,
		type = match.arg(type),
		conf_int = conf_int,
		alpha = alpha
	)
}

#' Plot CATT or event-study estimates from a fitted BETWFE
#'
#' @description
#' Parallel to [plot.fetwfe()]. BETWFE uses the same bridge-penalty
#' selection mechanism, so `selected` is encoded the same way (TRUE
#' = nonzero, FALSE = zeroed by the bridge penalty).
#'
#' @inheritParams plot.fetwfe
#' @return A `ggplot` object.
#' @seealso [plot.fetwfe()] for the full documentation; [eventStudy()];
#'   [cohortStudy()].
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   res <- betwfeWithSimulatedData(dat)
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot(res)
#'   }
#' }
#' @export
plot.betwfe <- function(
	x,
	type = c("event_study", "catt"),
	conf_int = TRUE,
	alpha = NULL,
	...
) {
	.plot_estimator(
		x,
		type = match.arg(type),
		conf_int = conf_int,
		alpha = alpha
	)
}

# ----------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------

# Dispatch to per-view builder; check ggplot2 availability and run the
# shared estimator-object precondition guard (`.check_for_plot()` in
# R/class_helpers.R asserts a fitted-class object, mirroring the
# precondition the predecessor event-study-only plot helpers had).
.plot_estimator <- function(x, type, conf_int, alpha) {
	.check_for_plot(x)
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop(
			"plot() for fetwfe / etwfe / betwfe results requires the ",
			"`ggplot2` package. Install it with ",
			"`install.packages(\"ggplot2\")`.",
			call. = FALSE
		)
	}
	if (type == "catt") {
		return(.plot_catt(x, conf_int = conf_int, alpha = alpha))
	}
	# type == "event_study" (validated by match.arg upstream)
	.plot_event_study(x, conf_int = conf_int, alpha = alpha)
}

# Per-cohort ATT plot. One point per cohort, with optional CI error bars.
# Selected/unselected cohorts (for fetwfe / betwfe) are shape- and
# color-coded.
.plot_catt <- function(x, conf_int, alpha) {
	catt <- x$catt_df
	# The catt_df S3 class adds helpful-error interceptors on $/`[`/`[[`
	# for the pre-1.11.0 Title-Case names. Aesthetics access modern
	# snake_case columns, so no interception fires.
	#
	# If the user passed an explicit alpha, recompute CIs from
	# estimate +/- qnorm(1 - alpha/2) * se. Otherwise use the CI columns
	# already in catt_df (computed at fit time using the fit's alpha).
	if (!is.null(alpha)) {
		z <- stats::qnorm(1 - alpha / 2)
		# Use unclass() so .subset2() / direct $-access can mutate the
		# data.frame without triggering catt_df's helpful-error layer on
		# new column writes (which doesn't intercept these but is safer
		# to avoid).
		catt$ci_low <- catt$estimate - z * catt$se
		catt$ci_high <- catt$estimate + z * catt$se
	}
	p <- ggplot2::ggplot(
		catt,
		ggplot2::aes(x = cohort, y = estimate)
	)
	if (isTRUE(conf_int)) {
		p <- p +
			ggplot2::geom_errorbar(
				ggplot2::aes(ymin = ci_low, ymax = ci_high),
				width = 0.2,
				color = "gray40"
			)
	}
	# Selection encoding: only for fetwfe / betwfe (catt_df has
	# `selected` column then). For etwfe / twfeCovs (no selection),
	# uniform styling.
	if ("selected" %in% names(catt)) {
		p <- p +
			ggplot2::geom_point(
				ggplot2::aes(shape = selected, color = selected),
				size = 3
			) +
			ggplot2::scale_shape_manual(
				values = c("FALSE" = 1, "TRUE" = 16),
				name = "Selected"
			) +
			ggplot2::scale_color_manual(
				values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
				name = "Selected"
			)
	} else {
		p <- p +
			ggplot2::geom_point(size = 3, color = "steelblue")
	}
	p +
		ggplot2::geom_hline(
			yintercept = 0,
			linetype = "dashed",
			color = "gray50"
		) +
		ggplot2::labs(
			x = "Cohort (first-treated period)",
			y = "Treatment effect estimate",
			title = "Cohort Average Treatment Effects"
		) +
		ggplot2::theme_minimal()
}

# Event-study plot. Points connected by a line, with optional CI bars.
# Forwards `alpha` to eventStudy() so the CIs reflect the requested
# confidence level.
.plot_event_study <- function(x, conf_int, alpha) {
	es <- eventStudy(x, alpha = alpha)
	p <- ggplot2::ggplot(
		es,
		ggplot2::aes(x = event_time, y = estimate)
	) +
		ggplot2::geom_line(color = "steelblue") +
		ggplot2::geom_point(size = 3, color = "steelblue")
	if (isTRUE(conf_int)) {
		p <- p +
			ggplot2::geom_errorbar(
				ggplot2::aes(ymin = ci_low, ymax = ci_high),
				width = 0.2,
				color = "gray40"
			)
	}
	p +
		ggplot2::geom_hline(
			yintercept = 0,
			linetype = "dashed",
			color = "gray50"
		) +
		ggplot2::labs(
			x = "Event time (periods since treatment)",
			y = "Treatment effect estimate",
			title = "Event-Study Average Treatment Effects"
		) +
		ggplot2::theme_minimal()
}
