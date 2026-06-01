# Parametric simultaneous confidence intervals over families of treatment
# effects for fetwfe / etwfe / betwfe / twfeCovs (#192).
#
# The simultaneous (family-wise) critical value c_{1 - alpha} is the
# (1 - alpha) quantile of max_k |Z_k| under Z ~ N(0, cov2cor(Sigma)), computed
# deterministically via mvtnorm::qmvnorm(). The joint covariance Sigma is
# reconstructed at call time from the fit's stored slots (mirroring how
# eventStudy() reconstructs per-event-time variances) using the K-effect
# generalizations of the package's variance machinery
# (.assemble_joint_cov_var1 / .assemble_joint_cov_var2 in
# R/variance_machinery.R). See .plans/feat-simultaneous-cis-192/PLAN.md.

# Bare column names referenced inside ggplot2::aes() in
# plot.simultaneous_cis() are intentional (NSE); declare them here so R CMD
# check doesn't emit a "no visible binding for global variable" NOTE (matches
# the pattern in R/plot.R and R/event_study.R).

#' @importFrom mvtnorm qmvnorm GenzBretz
NULL

utils::globalVariables(c(
	"effect_idx",
	"estimate",
	"effect",
	"simultaneous_ci_low",
	"simultaneous_ci_high",
	"pointwise_ci_low",
	"pointwise_ci_high"
))

#' Parametric simultaneous (1 - alpha) confidence intervals over a family of
#' treatment effects
#'
#' @description
#' Computes simultaneous (family-wise) confidence intervals for a user-
#' specified family of treatment effects from a fitted FETWFE / ETWFE /
#' BETWFE / twfeCovs object. The simultaneous critical value `c_{1 - alpha}`
#' is the `(1 - alpha)` quantile of `max_k |Z_k|` where `Z` follows a
#' multivariate normal with correlation matrix `cov2cor(Sigma)`; it is
#' computed deterministically via `mvtnorm::qmvnorm()`. Under Faletto (2025)
#' Theorem (c') tight Gaussianity and Assumption (Psi-IF), the family of
#' psi-linear effects is asymptotically multivariate normal with a covariance
#' that is estimable from the package's existing variance machinery; under the
#' paper's fixed-dim framing no high-dimensional correction is needed.
#'
#' The pointwise critical value `qnorm(1 - alpha/2)` (per-effect coverage) and
#' the Bonferroni-conservative critical value `qnorm(1 - alpha/(2K))` (family-
#' wise coverage with no correlation assumption) are returned for side-by-side
#' comparison; the simultaneous critical value is always between them when the
#' effects are positively correlated (as is typical in difference-in-
#' differences, where effects share the regression-coefficient variance piece).
#'
#' @param result A fitted object of class `"fetwfe"`, `"etwfe"`, `"betwfe"`,
#'   or `"twfeCovs"`.
#' @param family Character; one of `"event_study"`, `"cohort"`,
#'   `"all_post_treatment"`, or `"custom"`. See Details for each family's
#'   resolution.
#' @param alpha Numeric in `(0, 1)`; significance level. Default `0.05`.
#' @param contrasts For `family = "custom"`, a `K x num_treats` matrix whose
#'   rows give the `K` linear combinations of the underlying per-`(r, t)`
#'   treatment-effect vector (the `multcomp::glht()` convention). Ignored for
#'   the other families.
#' @return An object of S3 class `"simultaneous_cis"`: a list with
#'   \describe{
#'     \item{ci}{A data frame with columns `effect`, `estimate`,
#'       `simultaneous_ci_low`, `simultaneous_ci_high`, `pointwise_ci_low`,
#'       `pointwise_ci_high` (one row per effect in the family).}
#'     \item{critical_value}{The simultaneous critical value `c_{1 - alpha}`
#'       (or, under `se_type = "conservative"`, the Bonferroni critical value
#'       `qnorm(1 - alpha/(2K))` -- see Details).}
#'     \item{pointwise_critical_value}{`qnorm(1 - alpha/2)`, for reference.}
#'     \item{bonferroni_critical_value}{`qnorm(1 - alpha/(2K))`, for reference.}
#'     \item{family}{The requested family (character).}
#'     \item{alpha}{The significance level used.}
#'     \item{K}{The number of effects in the family (integer).}
#'   }
#' @details
#' **Family resolution and `K`.** `"event_study"` resolves to one effect per
#' post-treatment event time `e = 0, ..., T - 2` (`K = T - 1`); `"cohort"` to
#' one effect per treated cohort (`K = R`); `"all_post_treatment"` to one
#' effect per `(r, t)` cell (`K = num_treats`); `"custom"` to the
#' `K = nrow(contrasts)` user-supplied contrasts.
#'
#' **Joint covariance.** The `K x K` covariance `Sigma = Sigma_1 + Sigma_2` is
#' reconstructed at call time from the fit's stored slots (design matrix,
#' selected support, `theta_hat` / `beta_hat`, `cohort_probs_overall`,
#' `sig_eps_sq`). `Sigma_1` is the regression-coefficient piece and `Sigma_2`
#' the cohort-probability piece, generalizing the package's per-point variance
#' machinery (the same machinery `eventStudy()` uses). By construction
#' `sqrt(diag(Sigma))` equals the package's existing per-point standard errors
#' for the corresponding effects. The `Sigma` blocks are not persisted on the
#' fit; re-derivation is sub-second.
#'
#' **Paper grounding.** Theorem (c') tight Gaussianity (Faletto 2025,
#' `paper_arxiv.tex:1233`) guarantees the joint asymptotic normality; Assumption
#' (Psi-IF) (assumption equation `paper_arxiv.tex:2013`; in-prose discussion at
#' paper line 1268) is the influence-function condition the package's default
#' cohort-sample-proportions estimator satisfies; the fixed-dim framing follows
#' the paper's AE point 1(d).
#'
#' **Conservative fallback.** Under `se_type = "conservative"` the function
#' falls back to Bonferroni-corrected pointwise CIs (the Cauchy-Schwarz upper
#' bound used for the conservative scalar SE does not generalize to a `K x K`
#' covariance matrix) and emits a brief `message()`. The `$critical_value`
#' field is set to the Bonferroni value in this branch.
#'
#' **Numerical integration.** The critical value is computed via
#' `mvtnorm::qmvnorm(..., algorithm = mvtnorm::GenzBretz())` (mvtnorm's default
#' quasi-Monte Carlo integrator; sub-second through `K` up to about 100, so no
#' `K` cap is of practical concern for FETWFE families). `mvtnorm` is an
#' `Imports` dependency (as of version 1.16.0, when simultaneous bands became
#' the default reported confidence interval; see the `ci_type` argument of
#' [fetwfe()]). The function uses it only when `K > 1` and
#' `se_type != "conservative"` (the `K = 1` and conservative paths bypass the
#' dependency), and retains a defensive `stop()` with an actionable message if
#' it is somehow unavailable (e.g., a corrupted install).
#'
#' **Determinism contract.** The function is deterministic in its inputs: the
#' same fit plus the same `family`, `alpha`, and `contrasts` always produces
#' the same critical value across calls. This is achieved by wrapping the
#' internal `mvtnorm::qmvnorm()` call with a save/restore of the caller's
#' `.Random.seed` and a fixed internal `set.seed(1L)` immediately before the
#' call. The function does NOT mutate the caller's `.Random.seed` (the
#' save/restore via `on.exit()` leaves the caller's RNG state identical pre- and
#' post-call), matching the convention adopted by
#' `R/fetwfe_core.R::getBetaCV()` in PR #181 / v1.13.5. Users do not need to
#' call `set.seed()` before `simultaneousCIs()` to get reproducible results,
#' and downstream RNG-using code observes no perturbation.
#' @references
#' Faletto, G. (2025). Fused Extended Two-Way Fixed Effects for Difference-in-
#'   Differences with Staggered Adoptions. arXiv:2312.05985.
#'
#' Hothorn, T., Bretz, F., & Westfall, P. (2008). Simultaneous Inference in
#'   General Parametric Models. \emph{Biometrical Journal} 50(3), 346-363.
#' @seealso [eventStudy()] for the per-point event-study estimates and Wald
#'   intervals that `family = "event_study"` provides simultaneous bands over.
#' @examples
#' \donttest{
#'   coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
#'   fit <- fetwfeWithSimulatedData(sim)
#'   sci <- simultaneousCIs(fit, family = "event_study", alpha = 0.05)
#'   print(sci)
#' }
#' @export
simultaneousCIs <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL
) {
	UseMethod("simultaneousCIs")
}

#' @export
simultaneousCIs.fetwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses
	)
}

#' @export
simultaneousCIs.etwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses
	)
}

#' @export
simultaneousCIs.betwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses
	)
}

#' @note `twfeCovs` is documented as biased in its own roxygen header. The
#'   simultaneous CIs returned here are still mathematically well-defined on
#'   the `twfeCovs` CATT vector under the package's variance machinery, but the
#'   underlying point estimates are not unbiased. Prefer
#'   `simultaneousCIs.fetwfe` / `.etwfe` / `.betwfe` for inference. `twfeCovs`
#'   estimates a single pooled effect per cohort (no per-`(r, t)` cell or
#'   event-time structure), so only `family = "cohort"` and
#'   `family = "custom"` are defined for it; `"event_study"` and
#'   `"all_post_treatment"` raise an error.
#' @export
simultaneousCIs.twfeCovs <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses
	)
}

#' @title Shared worker for simultaneousCIs()
#' @description Reconstructs the K x K joint covariance from the fit's stored
#'   slots, computes the simultaneous critical value via `mvtnorm::qmvnorm()`,
#'   and assembles the `simultaneous_cis` result. Class-aware shims read the
#'   FETWFE-vs-OLS-family slot shapes. See the file-level comment and
#'   `.plans/feat-simultaneous-cis-192/PLAN.md`.
#' @param x A validated fitted estimator object.
#' @param family Character (already `match.arg`-resolved by the caller).
#' @param alpha Numeric significance level.
#' @param contrasts For `family = "custom"`, the `K x num_treats` contrast
#'   matrix; `NULL` otherwise.
#' @param has_valid_ses Logical; from `.check_for_simultaneous_cis()`.
#' @return An object of S3 class `"simultaneous_cis"`.
#' @keywords internal
#' @noRd
.simultaneous_cis_impl <- function(
	x,
	family,
	alpha,
	contrasts,
	has_valid_ses
) {
	# --- 1. Argument validation (mvtnorm guard deferred to step 11). ---
	stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)

	# --- 2. Read fit slots via class-aware shims. ---
	is_fetwfe <- inherits(x, "fetwfe")
	X_final <- if (is_fetwfe) x$internal$X_final else x$X_final
	y_final <- if (is_fetwfe) x$internal$y_final else x$y_final
	theta_hat_full <- if (is_fetwfe) x$internal$theta_hat else NULL
	beta_hat <- x$beta_hat
	treat_inds <- x$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs_overall <- x$cohort_probs_overall
	sig_eps_sq <- x$sig_eps_sq
	N <- x$N
	T_ <- x$T
	R <- x$R
	p <- x$p
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	tes <- beta_hat[treat_inds]

	# --- 3. Resolve cohort-time offsets + first_inds (same helper eventStudy
	#        uses; preserves scattered-cohort support). ---
	offs <- .resolve_event_study_offsets_and_first_inds(x, R = R, T = T_)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds

	# twfeCovs estimates a single pooled treatment effect PER COHORT (so
	# `treat_inds` has length R = num_treats, and each effect is already a
	# cohort ATT), not one effect per (r, t) cell. There is therefore no
	# event-time or per-cell structure to expand: only the `cohort` and
	# `custom` families are well-defined (this is why eventStudy() excludes
	# twfeCovs, R/event_study.R:73-75). Override first_inds to the per-cohort
	# singletons and reject the per-cell families with a clear message.
	if (inherits(x, "twfeCovs")) {
		if (!family %in% c("cohort", "custom")) {
			stop(
				"simultaneousCIs(): family = '",
				family,
				"' is not defined for a twfeCovs object, which estimates a ",
				"single pooled effect per cohort (no per-(r, t) cell or ",
				"event-time structure). Use family = 'cohort' or ",
				"family = 'custom'.",
				call. = FALSE
			)
		}
		first_inds <- seq_len(R)
	}

	# --- 4. Build the K x num_treats psi_tes matrix (row k = the contrast
	#        that picks effect k out of the (r, t) treatment-effect vector). ---
	psi_tes_mat <- .build_psi_tes_for_family(
		family = family,
		contrasts = contrasts,
		R = R,
		T = T_,
		num_treats = num_treats,
		cohort_offsets_int = cohort_offsets_int,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall
	)
	K <- nrow(psi_tes_mat)
	estimates <- as.numeric(psi_tes_mat %*% tes)

	pointwise_crit <- stats::qnorm(1 - alpha / 2)
	bonferroni_crit <- stats::qnorm(1 - alpha / (2 * K))
	effect_labels <- .effect_labels_for_family(
		family = family,
		K = K,
		R = R,
		T = T_,
		cohort_offsets_int = cohort_offsets_int
	)

	# --- 5. Edge case 1: no valid SEs -> stop(). ---
	if (!isTRUE(has_valid_ses)) {
		stop(
			"simultaneousCIs(): the fitted object was constructed with ",
			"calc_ses = FALSE; standard errors are not available. Re-fit ",
			"with q < 1 (FETWFE/BETWFE) and a satisfied rank condition.",
			call. = FALSE
		)
	}

	# --- 6. Determine the selected support (theta-space for FETWFE, beta-space
	#        for the OLS family + BETWFE). Mirrors event_study.R. ---
	if (is_fetwfe) {
		theta_hat_slopes <- theta_hat_full[2:(p + 1)]
		sel_feat_inds <- which(theta_hat_slopes != 0)
		sel_treat_inds_shifted <- which(theta_hat_slopes[treat_inds] != 0)
		theta_sel <- theta_hat_slopes[treat_inds][sel_treat_inds_shifted]
		d_inv_treat <- genInvTwoWayFusionTransformMat(
			n_vars = num_treats,
			first_inds = first_inds,
			R = R
		)
		d_inv_treat_sel <- if (length(sel_treat_inds_shifted) > 0) {
			d_inv_treat[, sel_treat_inds_shifted, drop = FALSE]
		} else {
			NULL
		}
	} else if (inherits(x, "betwfe")) {
		sel_feat_inds <- which(beta_hat != 0)
		sel_treat_inds_shifted <- which(beta_hat[treat_inds] != 0)
		# OLS-family unification: tes lives in beta-space, so the fusion-
		# inverse map is the identity restricted to the selected cells.
		theta_sel <- tes[sel_treat_inds_shifted]
		d_inv_treat_sel <- if (length(sel_treat_inds_shifted) > 0) {
			diag(num_treats)[, sel_treat_inds_shifted, drop = FALSE]
		} else {
			NULL
		}
	} else {
		# etwfe / twfeCovs: pure OLS, all cells selected.
		sel_feat_inds <- NA
		sel_treat_inds_shifted <- seq_len(num_treats)
		theta_sel <- tes
		d_inv_treat_sel <- diag(num_treats)
	}

	# --- 7. Degenerate support: all effects zeroed out by selection. Return
	#        zero estimates / zero SEs / pointwise critical value. ---
	if (length(sel_treat_inds_shifted) == 0L) {
		message(
			"simultaneousCIs(): the bridge penalty zeroed out every treatment ",
			"effect; estimates and standard errors are all zero."
		)
		ses <- rep(0, K)
		ci <- data.frame(
			effect = effect_labels,
			estimate = estimates,
			simultaneous_ci_low = estimates,
			simultaneous_ci_high = estimates,
			pointwise_ci_low = estimates,
			pointwise_ci_high = estimates,
			stringsAsFactors = FALSE
		)
		out <- list(
			ci = ci,
			critical_value = pointwise_crit,
			pointwise_critical_value = pointwise_crit,
			bonferroni_critical_value = bonferroni_crit,
			family = family,
			alpha = alpha,
			K = K
		)
		class(out) <- "simultaneous_cis"
		return(out)
	}

	# --- 8. Re-run getGramInv() and (if cluster) the cluster-robust sandwich,
	#        reusing the eventStudy machinery. ---
	gram_sel_feat <- if (any(!is.na(sel_feat_inds))) sel_feat_inds else NA
	gram_sel_treat <- if (any(!is.na(sel_feat_inds))) {
		sel_treat_inds_shifted
	} else {
		NA
	}
	res_gram <- getGramInv(
		N = N,
		T = T_,
		X_final = X_final,
		sel_feat_inds = gram_sel_feat,
		treat_inds = treat_inds,
		num_treats = num_treats,
		sel_treat_inds_shifted = gram_sel_treat,
		calc_ses = TRUE
	)
	gram_inv <- res_gram$gram_inv
	if (!isTRUE(res_gram$calc_ses)) {
		stop(
			"simultaneousCIs(): the Gram matrix on the selected support is not ",
			"invertible; the assumptions needed for inference are not ",
			"satisfied, so simultaneous confidence intervals cannot be ",
			"computed.",
			call. = FALSE
		)
	}

	sandwich_full <- NULL
	treat_block_mask <- NULL
	if (identical(se_type, "cluster")) {
		sel_arg <- if (any(!is.na(sel_feat_inds))) sel_feat_inds else NULL
		res_cl <- .assemble_cluster_robust_sandwich(
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T_,
			treat_inds = treat_inds,
			sel_feat_inds = sel_arg
		)
		sandwich_full <- res_cl$sandwich_full
		treat_block_mask <- res_cl$treat_block_mask
	}

	# --- 9. Shape Psi (p_sel x K) for Sigma_1, and the per-effect J_list for
	#        Sigma_2. `psi_tes_mat %*% d_inv_treat_sel` maps each effect's
	#        per-cell contrast into the selected support (theta-space for
	#        FETWFE, beta-space for the OLS family). ---
	Psi <- t(psi_tes_mat %*% d_inv_treat_sel)

	Sigma_1 <- .assemble_joint_cov_var1(
		Psi = Psi,
		gram_inv = gram_inv,
		sig_eps_sq = sig_eps_sq,
		N = N,
		T = T_,
		se_type = se_type,
		sandwich_full = sandwich_full,
		treat_block_mask = treat_block_mask
	)

	# --- 10. Sigma_2: per-effect Jacobian list. The cohort-probability piece
	#         is non-zero only when an effect pools across cohorts (its valid
	#         set V_k has more than one cohort). ---
	J_list <- .build_j_list_for_family(
		family = family,
		K = K,
		R = R,
		T = T_,
		num_treats = num_treats,
		cohort_offsets_int = cohort_offsets_int,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall,
		d_inv_treat_sel = d_inv_treat_sel
	)
	Sigma_pi_hat <- .multinomial_cov(cohort_probs_overall[1:R])
	Sigma_2 <- .assemble_joint_cov_var2(
		J_list = J_list,
		theta_sel = theta_sel,
		Sigma_pi_hat = Sigma_pi_hat,
		N = N,
		T = T_
	)

	# --- 11. Combine + critical value. ---
	if (is_indep || !identical(se_type, "conservative")) {
		# Tight Gaussian (default / cluster / indep): Sigma = Sigma_1 + Sigma_2.
		Sigma <- Sigma_1 + Sigma_2
		ses <- sqrt(pmax(diag(Sigma), 0))

		# Degenerate (zero-variance) effects: when the bridge penalty zeroes
		# out an effect's entire contribution to the selected support its SE is
		# exactly 0, its CI collapses to a point at the estimate, and it
		# contributes no family-wise risk. Such effects cannot enter the
		# correlation matrix (cov2cor() needs positive diagonal entries), so
		# the simultaneous critical value is computed over the non-degenerate
		# sub-family. A relative tolerance mirrors the existing rank-deficiency
		# threshold style in the variance machinery.
		var_tol <- .Machine$double.eps^0.5 * max(diag(Sigma), 1)
		nondeg <- diag(Sigma) > var_tol
		if (sum(nondeg) <= 1L) {
			# Zero or one non-degenerate effect: no joint correlation to
			# integrate; the per-effect critical value is exact.
			crit <- pointwise_crit
		} else {
			Sigma_nd <- Sigma[nondeg, nondeg, drop = FALSE]
			# Numerical-instability guard around cov2cor().
			rho <- tryCatch(
				stats::cov2cor(Sigma_nd),
				error = function(e) {
					stop(
						"simultaneousCIs(): could not convert the joint ",
						"covariance to a correlation matrix (possibly rank-",
						"deficient). Original error: ",
						conditionMessage(e),
						call. = FALSE
					)
				}
			)
			if (any(!is.finite(rho)) || any(abs(rho) > 1 + 1e-8)) {
				stop(
					"simultaneousCIs(): the joint covariance produced an ",
					"invalid correlation matrix; this is likely a rank-",
					"deficiency or numerical-stability issue.",
					call. = FALSE
				)
			}
			if (!requireNamespace("mvtnorm", quietly = TRUE)) {
				stop(
					"simultaneousCIs(): the `mvtnorm` package is required for ",
					"K > 1. Install it with `install.packages(\"mvtnorm\")`.",
					call. = FALSE
				)
			}
			# Paper grounding (WORKFLOW_LESSONS §12): Theorem (c') tight
			# Gaussianity (Faletto 2025, paper_arxiv.tex:1233) gives the joint
			# Z = sqrt(N) * (theta_hat - theta_0) -> N(0, V); the simultaneous
			# (1 - alpha) critical value is the (1 - alpha) quantile of
			# max_k |Z_k / sqrt(V_kk)|, computed by mvtnorm::qmvnorm() over the
			# correlation matrix cov2cor(V) (Hothorn, Bretz, Westfall 2008).
			#
			# mvtnorm::GenzBretz() (mvtnorm's default; quasi-Monte Carlo) is
			# sub-second through K ~ 100. Byte-determinism across calls is
			# preserved by save/restoring the caller's .Random.seed + a fixed
			# internal set.seed(1L) immediately before the qmvnorm() call. This
			# matches the in-package precedent at PR #181 / v1.13.5
			# (R/fetwfe_core.R::getBetaCV() lines 503-524), which uses the same
			# pattern for internal RNG use without mutating caller-side RNG
			# state.
			old_rng <- if (
				exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
			) {
				.GlobalEnv$.Random.seed
			} else {
				NULL
			}
			on.exit(
				{
					if (is.null(old_rng)) {
						if (
							exists(
								".Random.seed",
								envir = .GlobalEnv,
								inherits = FALSE
							)
						) {
							rm(".Random.seed", envir = .GlobalEnv)
						}
					} else {
						assign(".Random.seed", old_rng, envir = .GlobalEnv)
					}
				},
				add = TRUE
			)
			set.seed(1L)
			qmv <- mvtnorm::qmvnorm(
				p = 1 - alpha,
				corr = rho,
				tail = "both.tails",
				algorithm = mvtnorm::GenzBretz()
			)
			crit <- qmv$quantile
		}
	} else {
		# Conservative same-data: the Cauchy-Schwarz upper bound does not
		# generalize to a K x K matrix. Fall back to Bonferroni-corrected
		# pointwise CIs on the conservative scalar SEs.
		v1 <- pmax(diag(Sigma_1), 0)
		v2 <- pmax(diag(Sigma_2), 0)
		ses <- sqrt(v1 + v2 + 2 * sqrt(v1 * v2))
		crit <- if (K == 1L) pointwise_crit else bonferroni_crit
		message(
			"simultaneousCIs(): under se_type = 'conservative' the Cauchy-",
			"Schwarz upper bound does not extend to a joint K x K covariance ",
			"matrix; falling back to Bonferroni-corrected pointwise CIs. To ",
			"use the tight simultaneous critical value, re-fit with the ",
			"default se_type (which is valid under (Psi-IF))."
		)
	}

	# --- 12. Assemble the K-row $ci data frame + result. ---
	ci <- data.frame(
		effect = effect_labels,
		estimate = estimates,
		simultaneous_ci_low = estimates - crit * ses,
		simultaneous_ci_high = estimates + crit * ses,
		pointwise_ci_low = estimates - pointwise_crit * ses,
		pointwise_ci_high = estimates + pointwise_crit * ses,
		stringsAsFactors = FALSE
	)

	out <- list(
		ci = ci,
		critical_value = crit,
		pointwise_critical_value = pointwise_crit,
		bonferroni_critical_value = bonferroni_crit,
		family = family,
		alpha = alpha,
		K = K
	)
	class(out) <- "simultaneous_cis"
	out
}

#' @title Recompute a fit's cohort-family CIs as simultaneous bands (fit-time)
#' @description Internal helper for the `ci_type = "simultaneous"` default
#'   (#197). Calls the #192 worker `.simultaneous_cis_impl()` for the cohort
#'   family and returns the simultaneous lower/upper bounds aligned to the
#'   fit's `catt_df` cohort row order. Degrades gracefully: when standard
#'   errors are unavailable (`calc_ses = FALSE`), the family is degenerate, or
#'   `K = 1`, the worker returns bounds equal to the existing pointwise bounds
#'   (no widening) and any error short-circuits to `NULL` (leave `catt_df`
#'   unchanged).
#' @param x A fully-classed `fetwfe`/`etwfe`/`betwfe`/`twfeCovs` object.
#' @param alpha Numeric; the alpha the fit's `catt_df` was built at (`x$alpha`
#'   for fetwfe/etwfe/betwfe, `0.05` for twfeCovs, which has no `alpha` slot).
#'   Passed explicitly so all four callers agree.
#' @param has_valid_ses Logical; `res$calc_ses` from the core (read off the
#'   classed object by the caller).
#' @return A list with elements `ci_low` and `ci_high` (numeric vectors of
#'   length `R`, aligned to `x$catt_df` cohort order), or `NULL` to signal
#'   "leave `catt_df` unchanged".
#' @keywords internal
#' @noRd
.apply_simultaneous_catt_band <- function(x, alpha, has_valid_ses) {
	# If no valid SEs, nothing to widen -> leave the pointwise bounds (which
	# are NA / 0 in this case) untouched.
	if (!isTRUE(has_valid_ses)) {
		return(NULL)
	}
	# Cohort family. Wrap in tryCatch: a rank-deficient or degenerate family
	# should NOT abort the fit -> fall back to pointwise (NULL).
	# suppressMessages(): under se_type = "conservative" the worker emits a
	# Bonferroni-substitution message(); we do NOT want unprompted chatter on
	# every conservative fit (verbose-gating convention). The user retains the
	# note by calling simultaneousCIs() directly. See Decision Log "Q1 / R2".
	sci <- tryCatch(
		suppressMessages(
			.simultaneous_cis_impl(
				x = x,
				family = "cohort",
				alpha = alpha,
				contrasts = NULL,
				has_valid_ses = TRUE
			)
		),
		error = function(e) NULL
	)
	if (is.null(sci)) {
		return(NULL)
	}
	# Positional alignment: `.build_psi_tes_for_family(family = "cohort")`
	# iterates r = 1:R in cohort-block order, the SAME order
	# `getCohortATTsFinal()` builds `catt_df`. So `sci$ci` is already
	# row-aligned to `catt_df`. The `effect` labels differ
	# ("Cohort <offset>" vs `catt_df$cohort` = `c_names`), so rely on
	# position, not labels; assert the row count defensively.
	if (nrow(sci$ci) != nrow(x$catt_df)) {
		return(NULL)
	}
	list(
		ci_low = sci$ci$simultaneous_ci_low,
		ci_high = sci$ci$simultaneous_ci_high
	)
}

#' @title Apply the fit's `ci_type` to its cohort-family bounds (fit-time)
#' @description Internal finalizer for the `ci_type` default (#197). Given a
#'   fully-classed estimator object, if `ci_type == "simultaneous"` it
#'   overwrites `catt_df`'s `ci_low`/`ci_high` with the cohort-family
#'   simultaneous band (via `.apply_simultaneous_catt_band()`), re-validates,
#'   and returns the updated object. For `ci_type == "pointwise"` (or when the
#'   band degrades to `NULL`) it is a no-op pass-through. Hoisted out of the
#'   four entry-point tails to avoid a 4-site copy (WORKFLOW_LESSONS §14).
#' @param out A fully-classed `fetwfe`/`etwfe`/`betwfe`/`twfeCovs` object.
#' @param alpha Numeric; the alpha the fit used (`0.05` for twfeCovs).
#' @return `out` with `catt_df` bounds overwritten when simultaneous, else
#'   `out` unchanged.
#' @keywords internal
#' @noRd
.finalize_ci_type <- function(out, alpha) {
	if (!identical(out$ci_type, "simultaneous")) {
		return(out)
	}
	band <- .apply_simultaneous_catt_band(
		out,
		alpha = alpha,
		has_valid_ses = out$calc_ses
	)
	if (!is.null(band)) {
		# Overwrite ONLY the two interval-bound columns; se / p_value /
		# selected / estimate are untouched. (The catt_df S3 layer only
		# intercepts the OLD Title-Case names; ci_low/ci_high fall through.)
		cd <- out$catt_df
		cd$ci_low <- band$ci_low
		cd$ci_high <- band$ci_high
		out$catt_df <- cd
	}
	# Re-validate the final object (defense-in-depth; the ci_type slot + the
	# widened-band contract C9/C10 are now in place). `.assert_estimator_object()`
	# is class-dispatched off `out`, so no per-class name is hardcoded here.
	.assert_estimator_object(out)
	out
}

#' @title Build the K x num_treats psi_tes contrast matrix for a family
#' @description Row `k` of the returned matrix is the contrast on the per-
#'   `(r, t)` treatment-effect vector (`beta_hat[treat_inds]`) that defines
#'   effect `k`. The point estimate for effect `k` is
#'   `psi_tes_mat[k, ] %*% tes`.
#' @param family,contrasts,R,T,num_treats,cohort_offsets_int,first_inds,cohort_probs_overall
#'   See `.simultaneous_cis_impl()`.
#' @return A numeric `K x num_treats` matrix.
#' @keywords internal
#' @noRd
.build_psi_tes_for_family <- function(
	family,
	contrasts,
	R,
	T,
	num_treats,
	cohort_offsets_int,
	first_inds,
	cohort_probs_overall
) {
	if (family == "custom") {
		if (is.null(contrasts)) {
			stop(
				"simultaneousCIs(): family = 'custom' requires a `contrasts` ",
				"matrix (K x num_treats).",
				call. = FALSE
			)
		}
		if (
			!is.matrix(contrasts) ||
				!is.numeric(contrasts) ||
				nrow(contrasts) < 1L ||
				ncol(contrasts) != num_treats ||
				any(!is.finite(contrasts))
		) {
			stop(
				"simultaneousCIs(): `contrasts` must be a numeric matrix with ",
				"at least one row, exactly num_treats = ",
				num_treats,
				" columns, and all finite entries.",
				call. = FALSE
			)
		}
		return(unname(contrasts))
	}

	if (family == "cohort") {
		# Effect r = cohort r's ATT = mean(tes over cohort r's block).
		psi_tes <- matrix(0, nrow = R, ncol = num_treats)
		for (r in 1:R) {
			block <- .cohort_block_inds(r, R, first_inds, num_treats)
			psi_tes[r, block] <- 1 / length(block)
		}
		return(psi_tes)
	}

	if (family == "all_post_treatment") {
		# One effect per (r, t) cell: the identity contrast.
		return(diag(num_treats))
	}

	# family == "event_study": one effect per post-treatment event time
	# e = 0, ..., T - 2. Effect e pools the cells (r, e) over the valid cohort
	# set V_e, with sample-cohort-size weights.
	event_times <- 0:(T - 2L)
	K <- length(event_times)
	psi_tes <- matrix(0, nrow = K, ncol = num_treats)
	for (kk in seq_len(K)) {
		e <- event_times[kk]
		V_e <- which(cohort_offsets_int <= T - e)
		if (length(V_e) == 0L) {
			next
		}
		probs_Ve <- cohort_probs_overall[V_e]
		S_V <- sum(probs_Ve)
		if (S_V <= 0) {
			next
		}
		weights_Ve <- probs_Ve / S_V
		for (j in seq_along(V_e)) {
			psi_tes[kk, first_inds[V_e[j]] + e] <- weights_Ve[j]
		}
	}
	psi_tes
}

#' @title Build the per-effect Jacobian list for a family's Sigma_2 block
#' @description Returns a length-K list of per-effect cohort-weight Jacobians
#'   (via `.build_jacobian()`), matching each family's valid cohort set. Only
#'   `family = "event_study"` produces non-zero Jacobians; the cohort and
#'   all-post-treatment families have no cohort-probability weighting in a
#'   single effect, so their Jacobians (and hence `Sigma_2`) are zero. For
#'   `family = "custom"`, `Sigma_2` is conservatively set to zero (the user's
#'   contrasts operate directly on the fixed per-cell effects).
#' @param family,K,R,T,num_treats,cohort_offsets_int,first_inds,cohort_probs_overall,d_inv_treat_sel
#'   See `.simultaneous_cis_impl()`.
#' @return A length-K list of numeric matrices (each `R x ncol(d_inv_treat_sel)`).
#' @keywords internal
#' @noRd
.build_j_list_for_family <- function(
	family,
	K,
	R,
	T,
	num_treats,
	cohort_offsets_int,
	first_inds,
	cohort_probs_overall,
	d_inv_treat_sel
) {
	p_sel <- ncol(d_inv_treat_sel)
	zero_jac <- matrix(0, nrow = R, ncol = p_sel)

	if (family == "event_study") {
		event_times <- 0:(T - 2L)
		return(lapply(seq_len(K), function(kk) {
			e <- event_times[kk]
			V_e <- which(cohort_offsets_int <= T - e)
			.build_jacobian(
				cohort_probs_overall = cohort_probs_overall,
				R = R,
				d_inv_treat_sel = d_inv_treat_sel,
				mode = "per_effect_masked",
				V_e = V_e,
				first_inds = first_inds,
				e = e
			)
		}))
	}

	# cohort / all_post_treatment / custom: no per-effect cohort-probability
	# weighting -> Sigma_2 = 0. A single cohort's own ATT (cohort family) and a
	# single (r, t) cell (all_post_treatment) do not depend on the estimated
	# cohort probabilities; custom contrasts operate on the fixed per-cell
	# effects. (Mirrors getCohortATTsFinal(), whose per-cohort SEs carry no
	# var_2 term.)
	rep(list(zero_jac), K)
}

#' @title Build the `effect` column labels for a family
#' @param family,K,R,T,cohort_offsets_int See `.simultaneous_cis_impl()`.
#' @return A character vector of length `K`.
#' @keywords internal
#' @noRd
.effect_labels_for_family <- function(
	family,
	K,
	R,
	T,
	cohort_offsets_int
) {
	if (family == "event_study") {
		return(paste0("e", 0:(T - 2L)))
	}
	if (family == "cohort") {
		return(paste0("Cohort ", cohort_offsets_int[seq_len(R)]))
	}
	if (family == "all_post_treatment") {
		return(paste0("cell_", seq_len(K)))
	}
	# custom
	paste0("contrast_", seq_len(K))
}

#' @export
print.simultaneous_cis <- function(x, ...) {
	cat("Parametric simultaneous (1 - alpha) confidence intervals\n")
	cat(sprintf("Family: %s, K = %d, alpha = %g\n", x$family, x$K, x$alpha))
	cat(sprintf(
		"Simultaneous critical value: %.4f (pointwise %.4f, Bonferroni %.4f)\n",
		x$critical_value,
		x$pointwise_critical_value,
		x$bonferroni_critical_value
	))
	cat("\n")
	print(x$ci, row.names = FALSE, right = TRUE)
	invisible(x)
}

#' @importFrom generics tidy
#' @export
tidy.simultaneous_cis <- function(x, ...) {
	data.frame(
		term = x$ci$effect,
		estimate = x$ci$estimate,
		conf.low = x$ci$simultaneous_ci_low,
		conf.high = x$ci$simultaneous_ci_high,
		pointwise.conf.low = x$ci$pointwise_ci_low,
		pointwise.conf.high = x$ci$pointwise_ci_high,
		stringsAsFactors = FALSE
	)
}

#' @note Requires the `ggplot2` package (in Suggests). Install via
#'   `install.packages("ggplot2")` if it is not already installed. Matches
#'   `R/plot.R`'s `.plot_estimator` precedent for ggplot2-dependent plot
#'   methods.
#' @export
plot.simultaneous_cis <- function(x, ...) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop(
			"plot.simultaneous_cis requires ggplot2. Install via ",
			"install.packages('ggplot2').",
			call. = FALSE
		)
	}
	df <- x$ci
	df$effect_idx <- seq_len(nrow(df))
	# Two interval widths drawn together: the wider simultaneous band
	# (error bars) and the narrower pointwise band (thick line range), so the
	# family-wise vs per-effect coverage trade-off is visible. Bare column
	# names inside aes() are declared via utils::globalVariables() at the top
	# of this file (matches R/plot.R).
	ggplot2::ggplot(
		df,
		ggplot2::aes(x = effect_idx, y = estimate)
	) +
		ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
		ggplot2::geom_errorbar(
			ggplot2::aes(
				ymin = simultaneous_ci_low,
				ymax = simultaneous_ci_high
			),
			width = 0.2
		) +
		ggplot2::geom_linerange(
			ggplot2::aes(
				ymin = pointwise_ci_low,
				ymax = pointwise_ci_high
			),
			linewidth = 1
		) +
		ggplot2::geom_point(size = 2) +
		ggplot2::scale_x_continuous(
			breaks = df$effect_idx,
			labels = df$effect
		) +
		ggplot2::labs(
			x = "Effect",
			y = "Estimate",
			title = sprintf(
				"Simultaneous (thin) vs pointwise (thick) %g%% CIs",
				100 * (1 - x$alpha)
			)
		)
}
