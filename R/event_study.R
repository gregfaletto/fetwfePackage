# Event-study aggregation and plotting for fetwfe / etwfe / betwfe.
# See vignettes/etwfe_betwfe_vignette.Rmd for usage; see
# `.plans/feat-event-study-plot/plan.md` for the design rationale.

# Bare column names referenced inside ggplot2::aes() in .plot_event_study()
# are intentional (NSE); declare them here so R CMD check doesn't emit a
# "no visible binding for global variable" NOTE.
utils::globalVariables(c("event_time", "estimate", "ci_low", "ci_high"))

#' Compute pooled event-time treatment-effect estimates
#'
#' @description
#' For a fitted object from `fetwfe()`, `etwfe()`, or `betwfe()`, computes the
#' pooled event-time treatment-effect estimates `tau_E(e)`, defined as
#' cohort-weighted averages of the cell-level treatment-effect estimates at
#' each post-treatment event time `e = t - g` (where `t` is calendar time and
#' `g` is the cohort's first-treated calendar time). Weights are
#' sample-cohort-size weights (matching `did::aggte(type = "dynamic")`
#' convention).
#'
#' Standard errors combine two terms, mirroring the package's existing
#' overall-ATT SE machinery: `var_1(e)` from regression-coefficient noise
#' (computed via the same `gram_inv` machinery the package uses for cohort
#' SEs, or the cluster-robust sandwich under `se_type = "cluster"`), and
#' `var_2(e)` from cohort-probability noise (analog of the existing
#' `getSecondVarTermOLS` / `getSecondVarTermDataApp` machinery, with the
#' multinomial Jacobian restricted to cohorts valid at event time `e`).
#' Combined as `sqrt(var_1 + var_2)` by default (asymptotically exact
#' under paper Theorem `te.asym.norm.thm`(\eqn{c'}) / Assumption (\eqn{\Psi}-IF), which
#' the package's default cohort-sample-proportions estimator satisfies);
#' the conservative Cauchy-Schwarz bound `sqrt(var_1 + var_2 +
#' 2 sqrt(var_1 * var_2))` is available via `se_type = "conservative"`
#' (for users with non-(\eqn{\Psi}-IF) propensity-score estimators). When
#' `indep_counts` was supplied at fit time, the tight formula applies
#' regardless of `se_type` (two-sample regime, Theorem (b)).
#'
#' @param x A fitted object of class `"fetwfe"`, `"etwfe"`, or `"betwfe"`.
#' @param alpha (Optional) Significance level for confidence intervals.
#'   Defaults to `x$alpha` (the alpha used at fit time).
#' @param ci_type (Optional) Character; one of `"simultaneous"` or
#'   `"pointwise"`, or `NULL` (default). Controls whether the returned
#'   `ci_low` / `ci_high` are the event-study-family simultaneous (family-wise,
#'   uniform) band or the per-event-time pointwise Wald intervals. `NULL`
#'   inherits the fit's stored `ci_type` (so a fit made with the default
#'   `ci_type = "simultaneous"` yields simultaneous event-study bands, which
#'   `print` / `summary` / `plot` then display); for objects fitted before
#'   version 1.16.0 (no `ci_type` slot) `NULL` falls back to `"pointwise"`. The
#'   `se` column is identical under both settings; the interval bounds, their
#'   critical-value multiplier, and the `p_value` differ (under
#'   `"simultaneous"` the `p_value` is the single-step max-T multiplicity-
#'   adjusted dual of the band, under `"pointwise"` the per-effect Wald
#'   p-value; #200). Degenerate event times (no valid contributing cohorts)
#'   carry `NA` bounds and `NA` `p_value` under both settings.
#'
#'   **High-dimensional (`p >= NT`):** the `"simultaneous"` band here is the
#'   analytic post-selection band on the selected support (or `NA` bounds on a
#'   `gls = FALSE` fit), *not* the uniformly-valid desparsified band. For the
#'   valid high-dimensional event-study band call [simultaneousCIs()] with
#'   `family = "event_study"`, `method = "bootstrap"` (Theorem
#'   `debiased.highdim.joint.thm`, validated near-nominally in Faletto 2025).
#' @return A data frame with class `c("eventStudy", "data.frame")` and
#'   columns:
#'   \describe{
#'     \item{event_time}{Integer; event time `e = t - g`, ranging from 0
#'       to `T - 2`.}
#'     \item{n_cohorts}{Integer; number of cohorts contributing to the
#'       pooled estimate at event time `e`.}
#'     \item{estimate}{Numeric; the pooled event-time ATT estimate.}
#'     \item{se}{Numeric; combined standard error.}
#'     \item{ci_low}{Numeric; lower bound of the (1 - alpha) Wald CI.}
#'     \item{ci_high}{Numeric; upper bound of the (1 - alpha) Wald CI.}
#'     \item{p_value}{Numeric; follows the fit's `ci_type`. Under
#'       `"pointwise"`, the two-sided Wald p-value
#'       (`2 * pnorm(-|estimate / se|)`); under `"simultaneous"` (the
#'       default), the single-step max-T multiplicity-adjusted (family-wise)
#'       p-value matching the simultaneous band (#200). `NA` when `se` is `0`
#'       or `NA`.}
#'   }
#'   Only post-treatment event times (`e >= 0`) are included; pre-treatment
#'   placebo periods would require an extended regression specification and
#'   are out of scope for this initial release.
#' @seealso [cohortStudy()] for the parallel per-cohort accessor;
#'   [cohortTimeATTs()] for the fully disaggregated per-(cohort, time) accessor.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
#'   res <- fetwfeWithSimulatedData(dat)
#'   eventStudy(res)
#' }
#' @export
eventStudy <- function(x, alpha = NULL, ci_type = NULL) {
	if (inherits(x, "fetwfe")) {
		return(.event_study_fetwfe(x, alpha, ci_type))
	}
	if (inherits(x, "etwfe") || inherits(x, "betwfe")) {
		return(.event_study_etwfe_betwfe(x, alpha, ci_type))
	}
	stop(
		"eventStudy() requires an object of class 'fetwfe', 'etwfe', or ",
		"'betwfe'. Got class: ",
		paste(class(x), collapse = ", "),
		".",
		call. = FALSE
	)
}

#' Resolve the effective `ci_type` for an `eventStudy()` call
#'
#' @description `NULL` inherits the fit's stored `ci_type` slot (or
#'   `"pointwise"` for pre-1.16.0 objects with no slot); a non-`NULL` value is
#'   `match.arg`-validated and overrides for that call. Factored so both
#'   `.event_study_*` dispatchers resolve it identically (#197).
#' @keywords internal
#' @noRd
.resolve_event_study_ci_type <- function(x, ci_type) {
	if (is.null(ci_type)) {
		if (is.null(x$ci_type)) {
			"pointwise"
		} else {
			x$ci_type
		}
	} else {
		match.arg(ci_type, c("simultaneous", "pointwise"))
	}
}

#' Resolve cohort offsets and first-treatment-effect indices
#'
#' The package-wide cohort-offset resolver: **not** event-study-specific (despite
#' living in this file), it is used by `event_study()`, `cohortTimeATTs()`,
#' `simultaneousCIs()`, and `debiasedATT()` (#318). Wraps
#' `.derive_cohort_offsets_from_fit(x)` and the conditional dispatch between
#' `getFirstInds(G, T)` (the consecutive-cohort assumption, used as a fall-back
#' when `cohort_probs` carry no integer-coercible names --- true of synthetic
#' genCoefs-based fixtures) and `getFirstIndsFromOffsets(...)` (the
#' scattered-cohort path used on real panels like `bacondecomp::divorce`).
#'
#' Extracted in v1.13.3 (#174) to eliminate a byte-identical 10-line block
#' that had appeared in both `.event_study_etwfe_betwfe()` and
#' `.event_study_fetwfe()`. Future changes to the fall-back contract
#' (e.g., adding a "cohort_probs names are factor levels" branch) land in
#' this one helper.
#'
#' @param x A fitted estimator object passing `.check_for_event_study(x)`.
#' @param G Integer; the number of treated cohorts.
#' @param T Integer; the number of time periods.
#' @return A list with two elements: `cohort_offsets_int` (integer vector of
#'   length `G`) and `first_inds` (integer vector of length `G`).
#' @keywords internal
#' @noRd
.resolve_cohort_offsets_and_first_inds <- function(x, G, T) {
	cohort_offsets_int <- .derive_cohort_offsets_from_fit(x)
	if (is.null(cohort_offsets_int)) {
		list(
			cohort_offsets_int = seq.int(2L, G + 1L),
			first_inds = getFirstInds(G = G, T = T)
		)
	} else {
		list(
			cohort_offsets_int = cohort_offsets_int,
			first_inds = getFirstIndsFromOffsets(
				cohort_offsets_int = cohort_offsets_int,
				T = T
			)
		)
	}
}

#' True per-event-time treatment effects from cohort-time cell effects
#'
#' @description The point-estimate aggregation of `.event_study_fetwfe()` (the
#'   `estimates[k]` loop) factored out for reuse by [getTes()]: given the
#'   per-`(g, t)` treatment-effect cells and the cohort structure, returns the
#'   event-time effects `tau_E(e) = sum_{g in V_e} w_g * cell(g, e)`, where `V_e`
#'   is the set of cohorts still observed at event time `e`, `w_g` is cohort
#'   `g`'s probability normalized within `V_e`, and
#'   `cell(g, e) = cell_effects[first_inds[g] + e]`. Feeding a fit's estimated
#'   cells reproduces `eventStudy(fit)$estimate` (verified to machine precision);
#'   feeding a DGP's true cells (`coefs$beta[treat_inds]`) gives the true
#'   event-time effects. An event time with no contributing cohort returns `NA`
#'   (the truth is undefined there; the estimator reports `0`).
#' @param cell_effects Numeric length `num_treats`; per-`(g, t)` effects in
#'   `treat_inds` order.
#' @param first_inds Integer length `G`; index within `cell_effects` of each
#'   cohort's event-time-0 cell.
#' @param cohort_offsets_int Integer length `G`; each cohort's adoption offset
#'   (calendar adoption time), used to decide which cohorts reach event time `e`.
#' @param cohort_probs Numeric length `G`; cohort probabilities (any positive
#'   scaling; normalized within `V_e`).
#' @param T Integer; number of time periods. Event times run `0:(T - 2)`.
#' @return Named numeric length `T - 1`; the event-time effects with names
#'   `as.character(0:(T - 2))`, `NA` where no cohort reaches that event time.
#' @keywords internal
#' @noRd
.true_event_time_effects <- function(
	cell_effects,
	first_inds,
	cohort_offsets_int,
	cohort_probs,
	T
) {
	event_times <- 0:(T - 2L)
	out <- rep(NA_real_, length(event_times))
	for (k in seq_along(event_times)) {
		e <- event_times[k]
		V_e <- which(cohort_offsets_int <= T - e)
		if (length(V_e) == 0L) {
			next
		}
		probs_Ve <- cohort_probs[V_e]
		S_V <- sum(probs_Ve)
		if (!is.finite(S_V) || S_V <= 0) {
			next
		}
		weights_Ve <- probs_Ve / S_V
		out[k] <- sum(weights_Ve * cell_effects[first_inds[V_e] + e])
	}
	names(out) <- as.character(event_times)
	out
}

#' Event-study aggregation for ETWFE / BETWFE
#' @keywords internal
#' @noRd
.event_study_etwfe_betwfe <- function(x, alpha = NULL, ci_type = NULL) {
	# Method-entry precondition (#86). Validates the fitted object and
	# derives `has_valid_ses` — the contract gate that fixes #73
	# (eventStudy reporting finite SEs when the fit's att_se is NA
	# for q >= 1).
	contract <- .check_for_event_study(x)

	if (is.null(alpha)) {
		alpha <- x$alpha
	}
	z <- stats::qnorm(1 - alpha / 2)
	# Resolve the band type (#197): NULL inherits the fit's stored ci_type
	# (pointwise for pre-1.16.0 objects); an explicit value overrides.
	ci_type <- .resolve_event_study_ci_type(x, ci_type)

	if (is.null(x$cohort_probs_overall)) {
		stop(
			"eventStudy(): cohort_probs_overall missing from object. Re-fit with etwfe() / betwfe() at version 1.7.0 or later."
		)
	}

	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	N <- x$N
	T <- x$T
	G <- x$G
	sig_eps_sq <- x$sig_eps_sq
	cohort_probs_overall <- x$cohort_probs_overall
	X_final <- x$X_final
	y_final <- x$y_final
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	# v1.13.3 (#174): if cohort labels on the fit are parseable as panel-
	# time values (and the fit exposes its `first_year`), use the actual
	# offsets so scattered-cohort panels (`bacondecomp::divorce` is the
	# canonical reproducer) work; otherwise fall back to the consecutive-
	# cohort assumption to preserve byte-identical results on synthetic
	# fixtures whose `cohort_probs` carry no integer-coercible names.
	# `cohort_offsets_int` is also used for the per-event-time validity
	# set `V_e <- which(cohort_offsets_int <= T - e)`; for consecutive
	# cohorts this reduces to the pre-#174 `seq_len(G) <= T - 1L - e`.
	offs <- .resolve_cohort_offsets_and_first_inds(x, G = G, T = T)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds
	tes <- beta_hat[treat_inds]

	# Determine the selected support
	if (inherits(x, "etwfe")) {
		sel_feat_inds <- NA
		sel_treat_inds_shifted <- seq_len(num_treats)
	} else {
		# betwfe: selection is in beta-space
		sel_feat_inds <- which(beta_hat != 0)
		sel_treat_inds_shifted <- which(beta_hat[treat_inds] != 0)
	}

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

	# Cluster-robust sandwich (recomputed from existing slots)
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

	max_event <- T - 2L
	event_times <- 0:max_event
	estimates <- numeric(length(event_times))
	ses <- numeric(length(event_times))
	n_cohorts <- integer(length(event_times))

	for (k in seq_along(event_times)) {
		e <- event_times[k]
		# v1.13.3 (#174): use the actual cohort offsets to determine
		# which cohorts contribute at event time `e` (cohort `g` has a
		# treatment-effect cell at `(g, e)` iff `cohort_offsets[g] + e
		# <= T`, equivalently `cohort_offsets[g] <= T - e`). For
		# consecutive offsets this reduces to the pre-#174
		# `seq_len(G) <= T - 1L - e`.
		V_e <- which(cohort_offsets_int <= T - e)
		n_cohorts[k] <- length(V_e)
		if (length(V_e) == 0L) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}

		# Cohort-size weights, renormalized over V_e (sample-cohort-size convention)
		probs_Ve <- cohort_probs_overall[V_e]
		S_V <- sum(probs_Ve)
		if (S_V <= 0) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}
		weights_Ve <- probs_Ve / S_V

		# psi_e_tes (length num_treats): weight at idx(g, e) for g in V_e, 0 elsewhere
		psi_e_tes <- numeric(num_treats)
		for (j in seq_along(V_e)) {
			psi_e_tes[first_inds[V_e[j]] + e] <- weights_Ve[j]
		}

		# Point estimate
		estimates[k] <- as.numeric(crossprod(psi_e_tes, tes))

		if (!calc_ses) {
			ses[k] <- NA_real_
			next
		}

		# var_1(e) on the selected support
		if (identical(se_type, "cluster")) {
			psi_full <- numeric(ncol(sandwich_full))
			psi_e_sel <- psi_e_tes[sel_treat_inds_shifted]
			psi_full[treat_block_mask] <- psi_e_sel
			# Issue #84 item 9: floor the cluster-sandwich quadratic form
			# at zero. The Liang-Zeger sandwich is PSD in exact arithmetic,
			# so `t(psi) %*% sandwich %*% psi >= 0` always — but the N-loop
			# (now a `rowsum` aggregation, see item 18) involves enough
			# floating-point summation that a near-zero quadratic form
			# could come out negative by a few ulps. That would NaN the
			# downstream sqrt() and surface as a baffling "NA SE" report.
			# Floor at zero so any near-zero quadratic form rounds to
			# exactly zero (the mathematically correct answer for a
			# psi orthogonal to the sandwich's image). Issue #139 layers
			# a two-tier (warning / error) diagnostic on top of the floor
			# via `.floor_cluster_quad()`; on well-conditioned data the
			# behavior is unchanged.
			var_1_e <- .floor_cluster_quad(
				as.numeric(t(psi_full) %*% sandwich_full %*% psi_full),
				"event_study_etwfe_betwfe/var_1_e"
			)
		} else {
			psi_e_sel <- psi_e_tes[sel_treat_inds_shifted]
			var_1_e <- sig_eps_sq *
				as.numeric(t(psi_e_sel) %*% gram_inv %*% psi_e_sel) /
				(N * T)
		}

		# var_2(e) via existing getSecondVarTermOLS with single-cell psi_mat and
		# masked cohort_probs_overall (zeros outside V_e).
		psi_e_mat <- matrix(0, nrow = num_treats, ncol = G)
		for (g in V_e) {
			psi_e_mat[first_inds[g] + e, g] <- 1
		}
		masked_probs <- numeric(G)
		masked_probs[V_e] <- cohort_probs_overall[V_e]
		if (length(V_e) == 1L) {
			# Single-cohort: var_2 vanishes by construction (weight = 1)
			var_2_e <- 0
		} else {
			var_2_e <- getSecondVarTermOLS(
				psi_mat = psi_e_mat,
				tes = tes,
				cohort_probs_overall = masked_probs,
				num_treats = num_treats,
				N = N,
				T = T,
				G = G
			)
		}

		# Combine the two variance pieces. Mirrors the overall-ATT logic in
		# `getTeResultsOLS()` / `getTeResults2()` (R/variance_machinery.R):
		# the tight Gaussian formula is the default for the same-data path
		# under (Psi-IF) (Theorem (c$'$), paper line 1233 onwards); the
		# Cauchy-Schwarz upper bound is opt-in via `se_type = "conservative"`
		# for the rare non-(Psi-IF) propensity-score-estimator case.
		ses[k] <- sqrt(.combine_att_variance(
			var_1_e,
			var_2_e,
			is_indep,
			se_type
		))
	}

	# Simultaneous band (#197): when ci_type == "simultaneous" and SEs are
	# available, override the pointwise bounds with the event-study-family
	# simultaneous band (degrading to pointwise on any error, and preserving
	# NA on degenerate rows). `calc_ses` here is the LOCAL value (FALSE if the
	# Gram on the selected support was singular), which is the right gate.
	sb <- if (identical(ci_type, "simultaneous") && calc_ses) {
		.event_study_simultaneous_bounds(x, alpha, estimates, ses)
	} else {
		NULL
	}
	.assemble_event_study_df(
		event_times,
		n_cohorts,
		estimates,
		ses,
		z,
		ci_low = sb$ci_low,
		ci_high = sb$ci_high,
		p_value = sb$adjusted_p_values
	)
}

#' Event-study aggregation for FETWFE
#' @keywords internal
#' @noRd
.event_study_fetwfe <- function(x, alpha = NULL, ci_type = NULL) {
	# Method-entry precondition (#86). Validates the fitted object and
	# derives `has_valid_ses` — the contract gate that fixes #73
	# (eventStudy reporting finite SEs when the fit's att_se is NA
	# for q >= 1).
	contract <- .check_for_event_study(x)

	if (is.null(alpha)) {
		alpha <- x$alpha
	}
	z <- stats::qnorm(1 - alpha / 2)
	# Resolve the band type (#197): NULL inherits the fit's stored ci_type
	# (pointwise for pre-1.16.0 objects); an explicit value overrides.
	ci_type <- .resolve_event_study_ci_type(x, ci_type)

	if (is.null(x$internal) || is.null(x$internal$theta_hat)) {
		stop(
			"eventStudy(): theta_hat missing from x$internal. Re-fit with fetwfe() at version 1.7.0 or later."
		)
	}

	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	N <- x$N
	T <- x$T
	G <- x$G
	d <- x$d
	p <- x$p
	sig_eps_sq <- x$sig_eps_sq
	cohort_probs_overall <- x$cohort_probs_overall
	X_final <- x$internal$X_final
	y_final <- x$internal$y_final
	theta_hat_full <- x$internal$theta_hat
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	# v1.13.3 (#174): offset resolution lives in the shared helper
	# `.resolve_cohort_offsets_and_first_inds()` at the top of this
	# file. The fall-back keeps synthetic-fixture results byte-identical
	# to the pre-#174 path.
	offs <- .resolve_cohort_offsets_and_first_inds(x, G = G, T = T)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds
	tes <- beta_hat[treat_inds]

	# Selected support in theta-space (slopes only; drop intercept)
	theta_hat_slopes <- theta_hat_full[2:(p + 1)]
	sel_feat_inds <- which(theta_hat_slopes != 0)
	sel_treat_inds_shifted <- which(theta_hat_slopes[treat_inds] != 0)
	theta_hat_treat_sel <- theta_hat_slopes[treat_inds][sel_treat_inds_shifted]

	# d_inv_treat: the treatment-block of D^{-1}, then restrict columns to selected.
	# A user-supplied custom fusion matrix (#236), if any, is stored inverted on
	# `x$internal$d_inv_treat`; pass it through so this accessor uses the SAME
	# transform as the fit (otherwise it would fall back to the built-in dispatch
	# and silently disagree with the fitted effects).
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
	} else {
		d_inv_treat_sel <- NULL
	}

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

	# Cluster-robust sandwich (recomputed in theta-space on the selected support)
	sandwich_full <- NULL
	treat_block_mask <- NULL
	if (identical(se_type, "cluster") && calc_ses) {
		res <- .assemble_cluster_robust_sandwich(
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T,
			treat_inds = treat_inds,
			sel_feat_inds = sel_feat_inds
		)
		sandwich_full <- res$sandwich_full
		treat_block_mask <- res$treat_block_mask
	}

	max_event <- T - 2L
	event_times <- 0:max_event
	estimates <- numeric(length(event_times))
	ses <- numeric(length(event_times))
	n_cohorts <- integer(length(event_times))

	for (k in seq_along(event_times)) {
		e <- event_times[k]
		# v1.13.3 (#174): see matching block in `.event_study_etwfe_betwfe()`
		# above for rationale.
		V_e <- which(cohort_offsets_int <= T - e)
		n_cohorts[k] <- length(V_e)
		if (length(V_e) == 0L) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}

		probs_Ve <- cohort_probs_overall[V_e]
		S_V <- sum(probs_Ve)
		if (S_V <= 0) {
			estimates[k] <- 0
			ses[k] <- NA_real_
			next
		}
		weights_Ve <- probs_Ve / S_V

		# Point estimate: weighted sum of per-cell tes (which already
		# incorporates d_inv_treat times theta_hat under the FETWFE fit).
		psi_e_tes <- numeric(num_treats)
		for (j in seq_along(V_e)) {
			psi_e_tes[first_inds[V_e[j]] + e] <- weights_Ve[j]
		}
		estimates[k] <- as.numeric(crossprod(psi_e_tes, tes))

		if (!calc_ses) {
			ses[k] <- NA_real_
			next
		}

		# psi_e_theta on the selected theta-space treatment support:
		# psi_e_theta[col] = sum_{g in V_e} weights_Ve[g] * d_inv_treat_sel[idx(g,e), col]
		psi_e_theta <- numeric(length(sel_treat_inds_shifted))
		for (j in seq_along(V_e)) {
			g <- V_e[j]
			psi_e_theta <- psi_e_theta +
				weights_Ve[j] * d_inv_treat_sel[first_inds[g] + e, ]
		}

		# var_1(e)
		if (identical(se_type, "cluster")) {
			psi_full <- numeric(ncol(sandwich_full))
			psi_full[treat_block_mask] <- psi_e_theta
			# Issue #84 item 9: floor the cluster-sandwich quadratic form
			# at zero. See the matching guard in `.event_study_etwfe_betwfe`
			# (just above) for rationale. Issue #139 layers a two-tier
			# (warning / error) diagnostic on top of the floor via
			# `.floor_cluster_quad()`; on well-conditioned data the
			# behavior is unchanged.
			var_1_e <- .floor_cluster_quad(
				as.numeric(t(psi_full) %*% sandwich_full %*% psi_full),
				"event_study_fetwfe/var_1_e"
			)
		} else {
			var_1_e <- sig_eps_sq *
				as.numeric(
					t(psi_e_theta) %*% gram_inv %*% psi_e_theta
				) /
				(N * T)
		}

		# var_2(e) via the FETWFE-specific helper
		var_2_e <- .event_study_var2_fetwfe(
			e = e,
			V_e = V_e,
			theta_hat_treat_sel = theta_hat_treat_sel,
			d_inv_treat_sel = d_inv_treat_sel,
			cohort_probs_overall = cohort_probs_overall,
			first_inds = first_inds,
			num_treats = num_treats,
			N = N,
			T = T,
			G = G
		)

		# Combine the two variance pieces (see matching comment in
		# `.event_study_etwfe_betwfe` above): tight Gaussian by default
		# under (Psi-IF), conservative bound only via `se_type =
		# "conservative"`.
		ses[k] <- sqrt(.combine_att_variance(
			var_1_e,
			var_2_e,
			is_indep,
			se_type
		))
	}

	# Simultaneous band (#197): see the matching block in
	# `.event_study_etwfe_betwfe()`. `calc_ses` is the LOCAL value (FALSE if
	# the selected-support Gram was singular).
	sb <- if (identical(ci_type, "simultaneous") && calc_ses) {
		.event_study_simultaneous_bounds(x, alpha, estimates, ses)
	} else {
		NULL
	}
	.assemble_event_study_df(
		event_times,
		n_cohorts,
		estimates,
		ses,
		z,
		ci_low = sb$ci_low,
		ci_high = sb$ci_high,
		p_value = sb$adjusted_p_values
	)
}

#' Per-event-time second-variance term for FETWFE event-study
#'
#' Structurally parallel to `getSecondVarTermDataApp` but with the
#' cohort-block time-averaging Jacobian replaced by a single-row
#' selection at `idx(g, e)`, and `cohort_probs_overall` masked to zero
#' outside `V_e`. When `|V_e| = 1` the Jacobian rows vanish exactly and
#' `var_2(e) = 0`.
#'
#' Design choice: this helper is structurally parallel to
#' `getSecondVarTermDataApp` in `R/variance_machinery.R` (the cohort-block time-
#' averaging is replaced by a single-row selection at `idx(g, e)`, and
#' `cohort_probs_overall` is masked to zero outside `V_e`). Both helpers
#' implement paper Theorem 6.3's Jacobian formula (`paper_arxiv.tex:2577-
#' 2592`) where the off-diagonal coefficient `J_{rs} = -pi_s / S^2` uses the
#' column-index marginal cohort probability. (Prior to v1.8.0 the off-
#' diagonal coefficient was indexed by the outer-loop row; see issue #46.)
#' @keywords internal
#' @noRd
.event_study_var2_fetwfe <- function(
	e,
	V_e,
	theta_hat_treat_sel,
	d_inv_treat_sel,
	cohort_probs_overall,
	first_inds,
	num_treats,
	N,
	T,
	G
) {
	# Masked cohort_probs_overall: zero outside V_e
	masked <- numeric(G)
	masked[V_e] <- cohort_probs_overall[V_e]
	S_V <- sum(masked)
	if (S_V <= 0 || length(V_e) <= 1L) {
		# |V_e| = 1: single-cohort weight is 1; Jacobian rows vanish.
		return(0)
	}

	# Multinomial Sigma_pi_hat on the masked cohort probabilities.
	Sigma_pi_hat <- .multinomial_cov(masked)

	# Per-event-time Jacobian: G rows, length(sel_treat_inds_shifted) cols.
	# Rows for g not in V_e are zero (and are zero-killed by Sigma_pi_hat).
	# Diagonal: (S_V - pi_g)/S_V^2 * d_inv_treat_sel[idx(g, e), ]
	# Off-diagonal: subtract pi_{g'}/S_V^2 * d_inv_treat_sel[idx(g', e), ]
	# for g' in V_e \ {g}. Off-diagonal coefficient uses the COLUMN-index
	# pi_{g'}, matching paper Theorem 6.3. Consolidated into `.build_jacobian()`
	# (#192, WORKFLOW_LESSONS §14 Class A); the per-effect-masked path returns
	# this construction byte-identically.
	jacobian_e <- .build_jacobian(
		cohort_probs_overall = cohort_probs_overall,
		G = G,
		d_inv_treat_sel = d_inv_treat_sel,
		mode = "per_effect_masked",
		V_e = V_e,
		first_inds = first_inds,
		e = e
	)

	T *
		as.numeric(
			t(theta_hat_treat_sel) %*%
				t(jacobian_e) %*%
				Sigma_pi_hat %*%
				jacobian_e %*%
				theta_hat_treat_sel
		) /
		(N * T)
}

#' Assemble the event-study output data frame
#'
#' @description Builds the `eventStudy`-classed data frame. By default the
#'   confidence-interval bounds are the pointwise Wald bounds
#'   `estimate +/- z * se`. When `ci_low` / `ci_high` are supplied (the
#'   `ci_type = "simultaneous"` path, #197), those override the pointwise
#'   computation --- the single-function-with-optional-args form keeps the two
#'   bound-construction paths in ONE place (WORKFLOW_LESSONS section 14, avoiding a
#'   sibling `_with_bounds()` copy). `se` is always the same pointwise quantity;
#'   `p_value` defaults to the pointwise Wald p but is overridden when the
#'   `p_value` argument is supplied --- the `ci_type = "simultaneous"`
#'   single-step max-T adjusted dual (#200), masked to `NA` on the same
#'   degenerate rows as the bounds.
#' @param event_times,n_cohorts,estimates,ses,z See the per-event-time loops.
#' @param ci_low,ci_high Optional numeric vectors (length `length(estimates)`)
#'   to use as the CI bounds instead of `estimate +/- z * se`. `NULL` (default)
#'   computes the pointwise bounds.
#' @param p_value Optional numeric vector (length `length(estimates)`) of
#'   adjusted p-values to use instead of the pointwise Wald p (the
#'   `ci_type = "simultaneous"` path). `NULL` (default) computes the pointwise
#'   Wald p-value.
#' @keywords internal
#' @noRd
.assemble_event_study_df <- function(
	event_times,
	n_cohorts,
	estimates,
	ses,
	z,
	ci_low = NULL,
	ci_high = NULL,
	p_value = NULL
) {
	if (is.null(ci_low)) {
		ci_low <- estimates - z * ses
	}
	if (is.null(ci_high)) {
		ci_high <- estimates + z * ses
	}
	if (is.null(p_value)) {
		p_value <- .compute_p_values(estimates, ses)
	}
	out <- data.frame(
		event_time = as.integer(event_times),
		n_cohorts = as.integer(n_cohorts),
		estimate = estimates,
		se = ses,
		ci_low = ci_low,
		ci_high = ci_high,
		p_value = p_value
	)
	class(out) <- c("eventStudy", "data.frame")
	out
}

#' Compute event-study-family simultaneous CI bounds at fit time
#'
#' @description Internal helper for the `ci_type = "simultaneous"` default
#'   (#197). Calls the #192 event-study-family worker
#'   `.simultaneous_cis_impl(family = "event_study")` and returns its
#'   `simultaneous_ci_low` / `simultaneous_ci_high`, aligned to event-time
#'   order `0:(T-2)` (the same order the `eventStudy()` loop produces).
#'   Returns `NULL` to fall back to pointwise on any error or shape mismatch.
#'   `suppressMessages()` silences the `se_type = "conservative"`
#'   Bonferroni-substitution `message()` (same verbose-gating rationale as the
#'   cohort helper).
#'
#'   BLOCKER FIX (#197 round 1): preserves `NA` on every degenerate event-time
#'   row (where `eventStudy()`'s own loop set `se = NA`). The worker routes a
#'   zero-variance / empty-valid-cohort-set effect through its zero-variance
#'   handling and returns `0/0`; `eventStudy()` represents it as `NA/NA`.
#'   Blindly mapping the worker's `0/0` would silently turn `NA/NA` into `0/0`
#'   --- a representation change contradicting the documented invariant
#'   "bounds are `NA` under both settings when SEs are unavailable". So
#'   wherever `ses` is not finite, the simultaneous bounds are forced back to
#'   `NA`.
#' @param x A fully-classed `fetwfe`/`etwfe`/`betwfe` object.
#' @param alpha Numeric significance level.
#' @param estimates,ses The per-event-time point estimates and SEs the
#'   `eventStudy()` loop already computed (used for the `NA`-preservation mask).
#' @return A list with `ci_low` / `ci_high` / `adjusted_p_values` (length
#'   `length(estimates)`), or `NULL` to fall back to pointwise. Under
#'   `ci_type = "simultaneous"` the adjusted p-values are the single-step max-T
#'   duals of the band (#200); degenerate rows are `NA` to match the bounds.
#' @keywords internal
#' @noRd
.event_study_simultaneous_bounds <- function(x, alpha, estimates, ses) {
	sci <- tryCatch(
		suppressMessages(
			.simultaneous_cis_impl(
				x = x,
				family = "event_study",
				alpha = alpha,
				contrasts = NULL,
				has_valid_ses = TRUE,
				# Silent band-attach helper (like .apply_simultaneous_catt_band):
				# keep the high-dim degenerate notice a message() rather than a
				# warning(). The user gets the warning by calling simultaneousCIs()
				# directly. (In practice an all-zeroed fit has all-NA event-time SEs
				# so this path is skipped; FALSE keeps it silent defensively.) #304.
				warn_degenerate_highdim = FALSE
			)
		),
		error = function(e) NULL
	)
	if (is.null(sci)) {
		return(NULL)
	}
	# Positional alignment: `.build_psi_tes_for_family(family = "event_study")`
	# iterates over event_times 0:(T-2), the SAME order the eventStudy() loop
	# produces. Assert the row count defensively before mapping.
	if (nrow(sci$ci) != length(estimates)) {
		return(NULL)
	}
	ci_low <- sci$ci$simultaneous_ci_low
	ci_high <- sci$ci$simultaneous_ci_high
	# BLOCKER FIX (#197 round 1): preserve NA on degenerate event-time rows.
	# `eventStudy()`'s loop sets se = NA -> NA bounds for an empty/degenerate
	# valid-cohort-set event time; the worker returns 0/0 there. Mapping the
	# worker's 0/0 would silently turn NA/NA into 0/0.
	na_rows <- !is.finite(ses)
	ci_low[na_rows] <- NA_real_
	ci_high[na_rows] <- NA_real_
	# Single-step max-T adjusted p-values (#200) -- the dual of the band.
	# Force NA on the same degenerate rows the bounds are NA'd on, mirroring
	# `.compute_p_values()` (which returns NA wherever se is NA).
	adjusted_p_values <- sci$adjusted_p_values
	adjusted_p_values[na_rows] <- NA_real_
	list(
		ci_low = ci_low,
		ci_high = ci_high,
		adjusted_p_values = adjusted_p_values
	)
}

# `plot.<class>` methods and the rendering helpers were consolidated
# into `R/plot.R` (#29 PR). The event-study-only versions that used to
# live here have been superseded by a richer dispatch that supports
# both `type = "catt"` and `type = "event_study"`.
