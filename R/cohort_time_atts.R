# Per-(cohort, time) ATT accessor: the fully disaggregated sibling of
# `cohortStudy()` (per-cohort) and `eventStudy()` (per-event-time). Returns one
# row per treatment-effect *cell* (cohort g first treated, calendar time t),
# i.e. the `num_treats` underlying parameters with no averaging over either
# axis. See issue #280 (umbrella gregfaletto/fetwfe#48) and
# `.plans/feat-cohort-time-atts-280/PLAN.md`.
#
# Supports fetwfe / etwfe / betwfe (the estimators with a per-cohort-time
# treatment-effect block), mirroring eventStudy(). twfeCovs is excluded: the
# two-way-fixed-effects-with-covariates estimator carries a single treatment
# effect per cohort (num_treats == G, no time disaggregation), so its finest
# granularity is already cohortStudy().
#
# Each cell is a single (g, e) with aggregation weight 1, so this is the
# eventStudy per-event-time machinery specialized to a singleton validity set
# `V_e = {g}`: the second variance term `var_2` (cohort-probability sampling)
# vanishes identically, and the per-cell SE is the pure regression `var_1`
# `sqrt(sig_eps_sq * psi' gram_inv psi / (N * T))`. The SE is driven entirely by
# `psi` (NOT gated on `estimate == 0`): a cell fused/zeroed away has an all-zero
# `psi`, hence `se = 0` and `(0, 0)` bounds, which is the correct degenerate
# answer rather than a missing one.

#' Per-(cohort, time) average treatment effects
#'
#' @description
#' Extracts the fully disaggregated treatment-effect estimates from a fitted
#' FETWFE / ETWFE / BETWFE object: one row for every `(cohort, time)` cell, with
#' no averaging over cohorts or over event time. This is the finest-grained view
#' of the estimated effects --- the `num_treats` underlying parameters
#' themselves --- complementing [cohortStudy()] (which averages each cohort's
#' cells over time) and [eventStudy()] (which averages over cohorts at each
#' event time).
#'
#' Like [eventStudy()], this accessor is not available for [twfeCovs()] objects:
#' that estimator has a single treatment-effect parameter per cohort (no
#' per-time disaggregation), so its finest granularity is already
#' [cohortStudy()].
#'
#' Standard errors are the per-cell regression standard errors
#' \eqn{\sqrt{\sigma_\varepsilon^2 \, \psi' G^{-1} \psi / (NT)}}, recomputed at
#' call time from the fit's stored design (the same Gram-matrix machinery
#' [eventStudy()] uses; nothing is added to the fitted object). Because each
#' cell is a single cohort-time parameter, the cohort-probability sampling
#' variance that contributes to the aggregated [cohortStudy()] /
#' [eventStudy()] standard errors is identically zero here, so a cell's SE is a
#' single coefficient's regression SE.
#'
#' Confidence intervals and p-values are \emph{pointwise} \eqn{1 - \alpha} Wald
#' quantities (`estimate +/- z * se`). For simultaneous (family-wise) bands
#' over the cell family, use
#' `simultaneousCIs(result, family = "all_post_treatment")`.
#'
#' @param result A fitted object from [fetwfe()], [etwfe()], or [betwfe()] (or
#'   their `*WithSimulatedData()` wrapper analogs, which return the same
#'   classes). [twfeCovs()] objects are not supported (see Description).
#' @param alpha Numeric in `(0, 1)`; the pointwise confidence level is
#'   `1 - alpha`. Defaults to the `alpha` stored on the fit.
#'
#' @return A data frame with class `c("cohortTimeATTs", "data.frame")`
#'   containing one row per `(cohort, time)` treatment-effect cell, sorted by
#'   cohort then time, with columns:
#'   \describe{
#'     \item{cohort}{Character; the cohort label (the calendar time at which the
#'       cohort first received treatment), matching [cohortStudy()].}
#'     \item{time}{Numeric; the calendar time of the cell, equal to the cohort's
#'       adoption time plus the event time (`0, 1, ...`). Real panels carry their
#'       actual calendar times. For synthetic `genCoefs()` / `simulateData()`
#'       fixtures (whose panel runs `1, ..., T`, so the stored first year is `1`)
#'       this coincides with the 1-based panel-time index. (Only a hand-built or
#'       legacy fit with no stored first year falls back to that panel-time index
#'       directly.)}
#'     \item{estimate}{Numeric; the cell's ATT estimate.}
#'     \item{se}{Numeric; the pointwise standard error. `0` for a cell zeroed
#'       out by the fusion/bridge penalty while other cells survive
#'       (`fetwfe()` / `betwfe()`); `NA` when standard errors are unavailable ---
#'       the fit was computed with `q >= 1`, the Gram matrix on the selected
#'       support is singular, or the penalty zeroed the *entire* treatment block
#'       (no cells selected, so there is no support to recompute the Gram from;
#'       this matches [eventStudy()]).}
#'     \item{ci_low, ci_high}{Numeric; the pointwise `1 - alpha` Wald bounds
#'       `estimate -/+ qnorm(1 - alpha/2) * se`. `(0, 0)` for a fused-away cell;
#'       `NA` when `se` is `NA`.}
#'     \item{p_value}{Numeric; the two-sided pointwise Wald p-value
#'       `2 * pnorm(-|estimate / se|)`. `NA` when `se` is `0` or `NA`.}
#'     \item{selected}{(`fetwfe()` / `betwfe()` only.) Logical; `TRUE` when the
#'       bridge penalty left the cell's estimate nonzero. Absent for `etwfe()`,
#'       which does not perform selection.}
#'   }
#'   Use `tidy(cohortTimeATTs(result))` (with the `broom` package loaded) to
#'   reshape to broom convention; see [tidy.cohortTimeATTs()].
#'
#' @details
#' The cell standard error is computed from `psi`, the cell's row of the
#' (selected) treatment-effect design --- for `fetwfe()` the relevant row of the
#' inverse fusion transform \eqn{D^{-1}} in the transformed (theta) coordinate
#' space, for `betwfe()` / `etwfe()` a unit selector in the original (beta)
#' coordinate space restricted to the selected support. It is
#' never gated on the point estimate: a cell whose estimate is exactly zero
#' because the penalty fused it away has an all-zero `psi` and therefore
#' `se = 0` (the correct degenerate value), while a cell whose estimate happens
#' to be near zero for other reasons still receives its proper nonzero SE.
#'
#' @seealso [cohortStudy()] for the per-cohort (time-averaged) accessor;
#'   [eventStudy()] for the per-event-time (cohort-averaged) accessor;
#'   [simultaneousCIs()] for simultaneous (family-wise) bands over the cell
#'   family (`family = "all_post_treatment"`); [tidy.cohortTimeATTs()] for
#'   broom-shape translation.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
#'   res <- fetwfeWithSimulatedData(dat)
#'   cta <- cohortTimeATTs(res)
#'   cta
#'   # Broom-shape translation:
#'   if (requireNamespace("broom", quietly = TRUE)) {
#'     broom::tidy(cta)
#'   }
#' }
#' @export
cohortTimeATTs <- function(result, alpha = NULL) {
	# Restrict to the estimators with a per-cohort-time treatment-effect block
	# (mirrors eventStudy()'s dispatcher). twfeCovs has one effect per cohort
	# (num_treats == G), so there is nothing to disaggregate over time.
	if (!inherits(result, c("fetwfe", "etwfe", "betwfe"))) {
		if (inherits(result, "twfeCovs")) {
			stop(
				"cohortTimeATTs() is not available for twfeCovs objects: the ",
				"two-way-fixed-effects-with-covariates estimator has a single ",
				"treatment-effect parameter per cohort (no per-time ",
				"disaggregation). Use cohortStudy() for its per-cohort effects.",
				call. = FALSE
			)
		}
		stop(
			"cohortTimeATTs() requires a fitted object from fetwfe(), etwfe(), ",
			"or betwfe(). Got class: ",
			paste(class(result), collapse = ", "),
			".",
			call. = FALSE
		)
	}

	# Validate the object and derive the SE-validity gate (the #73 contract:
	# SEs are NA when the fit was computed with q >= 1). Reuses the same
	# precondition eventStudy() / simultaneousCIs() use.
	contract <- .check_for_event_study(result)
	x <- result

	if (is.null(alpha)) {
		alpha <- x$alpha
	}
	if (
		!is.numeric(alpha) ||
			length(alpha) != 1L ||
			is.na(alpha) ||
			alpha <= 0 ||
			alpha >= 1
	) {
		stop(
			"cohortTimeATTs(): `alpha` must be a single number in (0, 1).",
			call. = FALSE
		)
	}
	z <- stats::qnorm(1 - alpha / 2)

	# ---- Shared quantities ------------------------------------------------
	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	N <- x$N
	T <- x$T
	G <- x$G
	sig_eps_sq <- x$sig_eps_sq
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	tes <- beta_hat[treat_inds]

	# Canonical cohort labels (adoption time), in g-order, straight from the
	# fit's per-cohort table so the `cohort` column matches cohortStudy().
	if (is.null(x$catt_df) || is.null(x$catt_df$cohort)) {
		stop(
			"cohortTimeATTs(): the fitted object has no `catt_df$cohort` slot. ",
			"Re-fit with a current version of fetwfe()/etwfe()/betwfe()/twfeCovs().",
			call. = FALSE
		)
	}
	c_names <- as.character(x$catt_df$cohort)
	stopifnot(length(c_names) == G)

	# Cohort offsets + first-treatment indices (same resolution eventStudy()
	# uses: real panels map integer-coercible cohort labels to panel-time
	# offsets; synthetic fixtures fall back to the consecutive-cohort
	# assumption).
	offs <- .resolve_event_study_offsets_and_first_inds(x, G = G, T = T)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds
	first_year <- x$internal$first_year

	# ---- Class-specific support + design ----------------------------------
	is_fetwfe <- inherits(x, "fetwfe")
	is_betwfe <- inherits(x, "betwfe")
	# Only the bridge estimators perform selection, so only they get a
	# `selected` column (mirrors getCohortATTsFinal()'s include_selected).
	do_selected <- is_fetwfe || is_betwfe

	# For FETWFE the per-cell psi lives in the transformed (theta) coordinate
	# space, mapped through the treatment block of the inverse fusion transform
	# D^{-1}. For the OLS-family and BETWFE it is a unit selector in the
	# original (beta) coordinate space.
	d_inv_treat_sel <- NULL
	if (is_fetwfe) {
		X_final <- x$internal$X_final
		y_final <- x$internal$y_final
		p <- x$p
		theta_hat_slopes <- x$internal$theta_hat[2:(p + 1)]
		sel_feat_inds <- which(theta_hat_slopes != 0)
		sel_treat_inds_shifted <- which(theta_hat_slopes[treat_inds] != 0)
		# Treatment block of D^{-1}, threading any user-supplied fusion matrix
		# (#236) so this accessor uses the SAME transform as the fit.
		d_inv_treat <- .gen_inv_treat_block(
			num_treats = num_treats,
			first_inds = first_inds,
			G = G,
			fusion_structure = x$fusion_structure,
			d_inv_treat = x$internal$d_inv_treat
		)
		if (length(sel_treat_inds_shifted) > 0) {
			d_inv_treat_sel <- d_inv_treat[,
				sel_treat_inds_shifted,
				drop = FALSE
			]
		}
	} else if (is_betwfe) {
		X_final <- x$X_final
		y_final <- x$y_final
		sel_feat_inds <- which(beta_hat != 0)
		sel_treat_inds_shifted <- which(beta_hat[treat_inds] != 0)
	} else {
		# etwfe: unpenalized OLS, no selection.
		X_final <- x$X_final
		y_final <- x$y_final
		sel_feat_inds <- NA
		sel_treat_inds_shifted <- seq_len(num_treats)
	}

	# ---- Recompute the Gram inverse on the selected support ---------------
	calc_ses <- contract$has_valid_ses
	gram_inv <- NULL
	if (calc_ses && length(sel_treat_inds_shifted) > 0) {
		res_gram <- getGramInv(
			N = N,
			T = T,
			X_final = X_final,
			sel_feat_inds = sel_feat_inds,
			treat_inds = treat_inds,
			num_treats = num_treats,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			calc_ses = TRUE
		)
		gram_inv <- res_gram$gram_inv
		calc_ses <- res_gram$calc_ses
	} else {
		calc_ses <- FALSE
	}

	# ---- Optional cluster-robust sandwich ---------------------------------
	sandwich_full <- NULL
	treat_block_mask <- NULL
	if (identical(se_type, "cluster") && calc_ses) {
		sel_arg <- if (any(!is.na(sel_feat_inds))) sel_feat_inds else NULL
		res <- .assemble_cluster_robust_sandwich(
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			treat_inds = treat_inds,
			sel_feat_inds = sel_arg
		)
		sandwich_full <- res$sandwich_full
		treat_block_mask <- res$treat_block_mask
	}

	# ---- Per-cell enumeration ---------------------------------------------
	# Pre-size to num_treats: every treatment-effect parameter is exactly one
	# (cohort, time) cell.
	cohorts <- character(num_treats)
	times <- numeric(num_treats)
	estimates <- numeric(num_treats)
	ses <- numeric(num_treats)
	row <- 0L

	for (g in seq_len(G)) {
		cell_inds <- .cohort_block_inds(g, G, first_inds, num_treats)
		for (m in seq_along(cell_inds)) {
			idx <- cell_inds[m]
			e <- m - 1L
			row <- row + 1L

			cohorts[row] <- c_names[g]
			# Calendar time of the cell = adoption time + event time. Every fit
			# produced through the public API stamps `first_year` (= the first
			# panel time; see design_matrix.R), so the calendar branch runs:
			# synthetic genCoefs()/simulateData() fixtures carry first_year = 1
			# (with consecutive offsets 2, ..., G + 1), so it reduces to the
			# 1-based panel-time index. The is.null() branch is a defensive
			# fall-back for hand-built/legacy objects lacking the slot.
			times[row] <- if (is.null(first_year)) {
				cohort_offsets_int[g] + e
			} else {
				first_year + cohort_offsets_int[g] + e - 1
			}
			estimates[row] <- tes[idx]

			if (!calc_ses) {
				ses[row] <- NA_real_
				next
			}

			# psi for this single cell (aggregation weight 1).
			if (is_fetwfe) {
				psi_cell <- d_inv_treat_sel[idx, ]
			} else {
				psi_full_treat <- numeric(num_treats)
				psi_full_treat[idx] <- 1
				psi_cell <- psi_full_treat[sel_treat_inds_shifted]
			}

			# var_1 only: var_2 (cohort-probability sampling) is identically
			# zero for a singleton cell, so the SE is psi' (.) psi.
			if (identical(se_type, "cluster")) {
				psi_full <- numeric(ncol(sandwich_full))
				psi_full[treat_block_mask] <- psi_cell
				# Floor the cluster-sandwich quadratic form at zero before
				# sqrt(); see the matching guard in .event_study_fetwfe().
				var_1 <- .floor_cluster_quad(
					as.numeric(t(psi_full) %*% sandwich_full %*% psi_full),
					"cohort_time_atts/var_1"
				)
			} else {
				var_1 <- sig_eps_sq *
					as.numeric(t(psi_cell) %*% gram_inv %*% psi_cell) /
					(N * T)
			}

			ses[row] <- sqrt(.combine_att_variance(
				var_1,
				0,
				is_indep,
				se_type
			))
		}
	}
	stopifnot(row == num_treats)

	# ---- Assemble -------------------------------------------------------
	out <- data.frame(
		cohort = cohorts,
		time = times,
		estimate = estimates,
		se = ses,
		ci_low = estimates - z * ses,
		ci_high = estimates + z * ses,
		p_value = .compute_p_values(estimates, ses),
		stringsAsFactors = FALSE
	)
	if (do_selected) {
		out$selected <- estimates != 0
	}

	class(out) <- c("cohortTimeATTs", "data.frame")
	out
}
