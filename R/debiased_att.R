# Debiased overall-ATT accessor (#291): a uniformly-valid complement to the
# fused plug-in interval (fit$att_hat / fit$att_se). The fused plug-in SE is the
# within-selection variance and under-covers the *aggregated* overall ATT; the
# debiased estimator restores nominal coverage at ETWFE efficiency. Hand-off from
# the FETWFE methodology paper, Theorem `debiased.att.thm` (eqs
# `debiased.att.def` / `debiased.ols.identity` / `debiased.att.se`). The exact
# algorithm mirrors the validated simulation reference
# `simulations/method_functions.R::debiased_fetwfe` in the paper repo.

#' Plug-in cohort-weight variance V2 (the `att_var_2` channel, without Omega)
#'
#' @description The debiased SE's cohort-weight channel
#'   `V2 = (1/N_tau) sum_g pi_hat_g (catt_g - att)^2` (paper eq `debiased.att.se`)
#'   --- the variance from estimating the cohort sample proportions. It needs only
#'   the per-cohort effects `catt_g`, the within-treated weights `pi_hat_g`, the
#'   overall ATT, and the treated-unit count `N_tau` --- NOT the noise variances
#'   or Omega. This equals the value the oracle SE machinery stores as
#'   `att_var_2` (verified identical, single-sample, including unequal cohorts);
#'   `debiasedATT()` recomputes it here when the fit carries no `att_var_2`
#'   (a `gls = FALSE` high-dimensional fit, #307).
#' @param fit A fitted `fetwfe` object.
#' @return Numeric scalar; the plug-in V2.
#' @keywords internal
#' @noRd
.plugin_v2 <- function(fit) {
	att <- fit$att_hat
	pi_g <- fit$cohort_probs
	catt <- fit$catt_df$estimate
	# N_tau = number of treated units = sum_g N_g, recovered from the overall
	# cohort proportions (`cohort_probs_overall` sums to the treated fraction).
	N_tau <- round(sum(fit$cohort_probs_overall) * fit$N)
	if (
		!is.numeric(att) ||
			length(att) != 1L ||
			is.na(att) ||
			!is.numeric(pi_g) ||
			!is.numeric(catt) ||
			length(pi_g) != length(catt) ||
			length(pi_g) == 0L ||
			anyNA(pi_g) ||
			anyNA(catt) ||
			!is.numeric(N_tau) ||
			length(N_tau) != 1L ||
			is.na(N_tau) ||
			N_tau <= 0
	) {
		stop(
			"debiasedATT(): could not compute the plug-in cohort-weight variance ",
			"(V2) from the fit (missing `catt_df` / `cohort_probs` / `att_hat`). ",
			"Re-fit with a current version of fetwfe().",
			call. = FALSE
		)
	}
	sum(pi_g * (catt - att)^2) / N_tau
}

#' Debiased overall-ATT with a uniformly-valid standard error
#'
#' @description
#' Computes a *debiased* estimate of the overall average treatment effect on the
#' treated (ATT) from a fitted [fetwfe()] object, together with a
#' uniformly-valid standard error. This is a complement to the fused plug-in
#' interval reported on the fit (`fit$att_hat` / `fit$att_se`): the fused SE is
#' the within-selection variance and tends to *under-cover* the aggregated
#' overall ATT (it omits the between-selection term), whereas `debiasedATT()`
#' returns an interval with asymptotically nominal coverage at roughly the
#' efficiency of unrestricted ETWFE.
#'
#' The debiased point estimate **differs** from the fused `fit$att_hat`: when
#' `p < NT`, by the OLS identity (paper eq. `debiased.ols.identity`) it equals the
#' unrestricted ETWFE/OLS estimate in the ATT direction. You therefore get a dual
#' offering --- `fit$att_hat` / `fit$att_se` (fused, efficient, pointwise) and
#' `debiasedATT(fit)` (debiased, ETWFE-efficient, uniformly valid). It is exposed
#' as a separate accessor (rather than a `se_type` option) precisely because the
#' debiased SE accompanies a different point estimate than the fused one.
#'
#' When `p >= NT` the unrestricted ETWFE is infeasible, so the OLS identity no
#' longer applies; the accessor instead builds the debiasing direction by a
#' nodewise (desparsified-lasso) relaxed inverse (paper Theorem
#' `debiased.highdim.thm`). The point estimate, SE formula, and return value are
#' otherwise identical --- "one estimator, two regimes." This high-dimensional
#' branch is **experimental** (see Assumptions).
#'
#' @details
#' The standard error has two channels that **add** under marginal cohort
#' assignment / Assumption (Psi-IF): a regression (outcome) channel,
#' per-unit-clustered, and a cohort-weight channel for the sampling noise in the
#' estimated cohort weights. See *Assumptions* below.
#'
#' The cohort-weight channel is the fit's own plug-in propensity variance
#' (`att_var_2`), the same quantity `fit$att_se` uses --- so it automatically
#' carries the correct single-sample formula and, for `indep_counts` fits, the
#' two-sample formula (the weight channel concerns the estimated cohort weights,
#' not the selection, so it is shared by the fused and debiased estimators). The
#' fit's `se_type` does **not** affect the regression channel ---
#' `debiasedATT()` always computes its own per-unit-clustered version. Fits made
#' with `add_ridge = TRUE` are **not** supported (the ridge-augmented design and
#' post-fit un-shrink are outside the validated construction); the accessor errors
#' on them.
#'
#' @section Assumptions:
#' The uniformly-valid SE relies on the hypotheses of paper Theorem
#' `debiased.att.thm`:
#' \itemize{
#'   \item **Marginal cohort assignment / (Psi-IF).** Cohort membership is
#'     independent of the outcome noise, so the regression and cohort-weight
#'     variance channels are uncorrelated and add. The package's default cohort
#'     weight estimator (sample proportions) satisfies this. Under *confounded*
#'     assignment (cohort membership correlated with treatment-effect
#'     heterogeneity) the two channels correlate and the additive SE can
#'     mis-state uncertainty; a combined per-unit score is then required (not yet
#'     implemented).
#'   \item **Correctly specified conditional mean** of the outcome.
#'     Neyman-orthogonality removes first-order sensitivity to the selected
#'     nuisance coefficients, but not to mean-misspecification. The regression
#'     channel is the unit-clustered sandwich, valid as the number of units
#'     `N -> infinity` with **no** model for within-unit dependence (so a
#'     `gls = FALSE` fit, which skips GLS whitening, yields a valid cluster-robust
#'     SE; #307, paper Decision D1). GLS whitening (`gls = TRUE`) buys
#'     *efficiency* (an asymptotically smaller `var_reg`), not validity.
#'   \item **(Low-dimensional `p < NT` only) Asymptotically negligible ridge,**
#'     `lambda = o((NT)^(-1/2))` (the theoretical condition; the leading case is
#'     the exact inverse `lambda = 0`). The implementation does **not** use a
#'     vanishing schedule --- it adds a fixed numerical stabilizer
#'     `1e-6 * mean(diag(Sigma))` to the Gram before solving, which is negligible
#'     and cancels in the OLS-identity case.
#'   \item **Growing number of clusters,** `N -> infinity`. The CLT is over the
#'     `N` independent units, not the `NT` rows. With few treated units the
#'     cluster approximation is poor; prefer a wild-cluster bootstrap when `N` is
#'     small.
#'   \item **Regularity / two regimes.** When `p < NT` (Theorem
#'     `debiased.att.thm`) the full design Gram is nonsingular and the debiasing
#'     direction is the exact inverse; the accessor reduces to debiased ETWFE.
#'     When `p >= NT` (Theorem `debiased.highdim.thm`) the Gram is singular and the
#'     direction is the **nodewise (desparsified-lasso) relaxed inverse** of
#'     equation `debiased.highdim.v` (the same estimate and SE, only `v` differs
#'     --- "one estimator, two regimes"); it requires sparsity
#'     (`s_N log(p_N) / sqrt(N) -> 0`), a restricted-eigenvalue condition over
#'     sparse cones, `||v*||_1 = O(1)`, and a bounded limiting variance
#'     `a' Sigma_theta^(-1) a`. The high-dimensional branch realizes Theorem
#'     `debiased.highdim.thm` with the theory-scaled penalty, but its
#'     finite-sample **coverage is not yet validated by simulation** (the penalty
#'     constant `lambda_c` is not yet tuned); treat the `p >= NT` path as
#'     **experimental** and inspect the returned `feasibility` / `converged`
#'     diagnostics.
#' }
#'
#' The two regimes use different nuisances. The **fixed-`p`** branch uses the
#' fit's bridge (`q < 1`) `theta_hat`: Theorem `debiased.att.thm` needs only
#' *consistency* of the nuisance, and with the exact inverse this is the OLS
#' identity. The **high-dimensional** branch instead fits an internal `q = 1`
#' fused lasso (`grpreg::cv.grpreg(penalty = "gBridge", gamma = 1)`, CV-selected
#' `lambda`) as the nuisance, because Theorem `debiased.highdim.thm` controls the
#' orthogonalization remainder through the nuisance's `l1` *rate*
#' (`||theta_hat - theta*||_1 = O_p(s_N lambda_theta)`) --- the rate the `q = 1`
#' fused lasso enjoys under a restricted-eigenvalue condition (the standard
#' desparsified-lasso rate; cf. Negahban et al. 2012), the `p >= NT` extension
#' the paper points to; the `q < 1` bridge is super-efficient / non-uniform ---
#' the wrong nuisance for a uniformly-valid CI (#303). The
#' *input* fit must still be `q < 1` (it supplies the cohort-weight variance
#' channel `att_var_2`); only the high-dimensional debiasing nuisance is `q = 1`.
#' That nuisance's cross-validation uses a fixed, data-derived seed, so it is
#' deterministic across calls (`debiasedATT()` and the high-dimensional
#' [simultaneousCIs()] band center agree exactly) --- but it is not tunable via
#' the fit's `cv_seed`.
#'
#' @param fit A fitted object from [fetwfe()]. In the fixed-p (`p < NT`) regime
#'   the SE needs either bridge selection (`q < 1`, `gls = TRUE`) or an un-whitened
#'   fit (`gls = FALSE`, any `q`, the Omega-free cluster-robust path; #312); a
#'   GLS-whitened `q >= 1` fit is rejected. The high-dimensional (`p >= NT`) path
#'   accepts any fit. [etwfe()] / [betwfe()] / [twfeCovs()] are not supported (the
#'   construction is specific to the FETWFE transformed-coefficient space).
#' @param alpha Numeric in `(0, 1)`; the confidence level is `1 - alpha`.
#'   Defaults to the `alpha` stored on the fit.
#' @param lambda_c The leading constant of the high-dimensional (`p >= NT`)
#'   nodewise penalty `lambda_node = lambda_c * max(|a|) * sqrt(log(p) / N)` (`N` =
#'   number of units). Either a single positive number (a fixed constant; default
#'   `1.0` = the theory scale) or the string `"cv"`, which selects the constant
#'   **per fit** by cross-validating the unit-level Riesz loss
#'   `0.5 v'Sigma v - a'v` over the KKT-feasible region (the desparsified-lasso /
#'   auto-DML standard; van de Geer; Chernozhukov-Newey-Singh), falling back to
#'   the theory scale `1.0` when no grid penalty is feasible (#295). Larger values
#'   shrink the debiasing direction more. **Ignored when `p < NT`.** The default
#'   stays the fixed `1.0` (the `"cv"` mechanism is opt-in pending simulation
#'   validation of its coverage); treat the `p >= NT` path as experimental.
#' @param riesz_max_iter,riesz_tol Integer / numeric; coordinate-descent controls
#'   for the high-dimensional nodewise solver. **Ignored when `p < NT`.**
#'
#' @return An object of S3 class `"debiased_att"` (a named list, with a `print`
#'   and a `tidy` method) with elements:
#'   \describe{
#'     \item{att}{Numeric scalar; the debiased overall-ATT point estimate.}
#'     \item{se}{Numeric scalar; the uniformly-valid standard error,
#'       `sqrt(var_reg + var_weight)`.}
#'     \item{ci_low, ci_high}{Numeric; the `1 - alpha` Wald interval
#'       `att +/- qnorm(1 - alpha/2) * se`.}
#'     \item{var_reg}{Numeric; the regression (outcome) channel's contribution to
#'       `se^2`, per-unit-clustered.}
#'     \item{var_weight}{Numeric; the cohort-weight channel's contribution to
#'       `se^2` --- the fit's plug-in propensity variance (`att_var_2`), which
#'       carries the correct single- or two-sample (`indep_counts`) formula.}
#'   }
#'   `var_reg` and `var_weight` are the two **`se^2`-scale** (unit-scale)
#'   components, so `se = sqrt(var_reg + var_weight)`. They are *not* the paper's
#'   CLT-scale variances `V_1^full` / `V_2` from Theorem `debiased.att.thm`, which
#'   are larger by a factor on the order of `NT`.
#'
#'   For high-dimensional (`p >= NT`) fits the list additionally contains
#'   `feasibility` (`= ||Sigma v - a||_inf`, the nodewise KKT certificate, which
#'   should be `<= lambda_node`), `converged` (the coordinate-descent flag),
#'   `lambda_node` (the penalty used), `lambda_c` (the leading constant used),
#'   and `lambda_c_selection` (`"fixed"` or `"cv"`). When `lambda_c = "cv"` it
#'   also contains `lambda_cv`, the cross-validation diagnostics (the grid, the
#'   per-grid feasibility flags and CV losses, and whether the theory-scale
#'   fallback fired). These are absent for `p < NT` fits.
#'
#' @seealso [fetwfe()] for the fused fit; [cohortTimeATTs()] / [eventStudy()] /
#'   [cohortStudy()] for the disaggregated fused effects.
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 1)
#'   dat <- simulateData(coefs, N = 200, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
#'   fit <- fetwfeWithSimulatedData(dat, q = 0.5)
#'   fit$att_hat # fused (pointwise)
#'   debiasedATT(fit) # debiased (uniformly valid)
#' }
#' @export
debiasedATT <- function(
	fit,
	alpha = NULL,
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	if (!inherits(fit, "fetwfe")) {
		stop(
			"debiasedATT() requires a fitted object from fetwfe(). The debiased ",
			"ATT construction is specific to FETWFE's transformed-coefficient ",
			"space and is not defined for etwfe()/betwfe()/twfeCovs(). Got class: ",
			paste(class(fit), collapse = ", "),
			".",
			call. = FALSE
		)
	}

	if (is.null(alpha)) {
		alpha <- fit$alpha
	}
	if (
		!is.numeric(alpha) ||
			length(alpha) != 1L ||
			is.na(alpha) ||
			alpha <= 0 ||
			alpha >= 1
	) {
		stop(
			"debiasedATT(): `alpha` must be a single number in (0, 1).",
			call. = FALSE
		)
	}

	# High-dimensional nodewise-solver controls. These only affect the p >= NT
	# path, but validate them unconditionally so a malformed value fails loudly
	# rather than silently producing a degenerate debiasing direction.
	# `lambda_c` is either the string "cv" (cross-validate the penalty constant per
	# fit, #295) or a single positive number (the fixed constant; 1.0 = theory
	# scale, the default and the CV fallback).
	if (
		!(identical(lambda_c, "cv") ||
			(is.numeric(lambda_c) &&
				length(lambda_c) == 1L &&
				!is.na(lambda_c) &&
				lambda_c > 0))
	) {
		stop(
			"debiasedATT(): `lambda_c` must be \"cv\" (cross-validate) or a single ",
			"positive number.",
			call. = FALSE
		)
	}
	if (
		!is.numeric(riesz_max_iter) ||
			length(riesz_max_iter) != 1L ||
			is.na(riesz_max_iter) ||
			riesz_max_iter < 1
	) {
		stop(
			"debiasedATT(): `riesz_max_iter` must be a single positive integer.",
			call. = FALSE
		)
	}
	if (
		!is.numeric(riesz_tol) ||
			length(riesz_tol) != 1L ||
			is.na(riesz_tol) ||
			riesz_tol <= 0
	) {
		stop(
			"debiasedATT(): `riesz_tol` must be a single positive number.",
			call. = FALSE
		)
	}

	X <- fit$internal$X_final
	y <- as.numeric(fit$internal$y_final)
	theta_hat <- fit$internal$theta_hat
	n <- nrow(X)
	p <- ncol(X)
	# Regime: p < NT (fixed-p, exact/ridged inverse) vs p >= NT (high-dimensional
	# nodewise desparsified direction, paper Theorem `debiased.highdim.thm`).
	highdim <- p >= n

	# Standard-error gate, regime-aware.
	#  - Fixed-p (p < NT): two valid SE paths. (i) A GLS-whitened bridge fit
	#    (`gls = TRUE`, q < 1, `calc_ses = TRUE`) uses the OLS-identity oracle SE.
	#    (ii) An un-whitened fit (`gls = FALSE`, `sig_eps_sq = NA`) uses the exact-
	#    inverse OLS direction + unit-clustered residual sandwich -- Omega-free and
	#    valid for ANY q, the fixed-p analog of the high-dim path (#312, paper
	#    Decision D1). Only a WHITENED q >= 1 fit (`calc_ses = FALSE` but
	#    `!is.na(sig_eps_sq)`) has no validated debiased SE and is rejected below.
	#  - High-dim (p >= NT): the SE is the desparsified cluster-robust sandwich
	#    (V1 unit-clustered + V2 plug-in), which needs neither the oracle SE
	#    machinery nor Omega -- so a `gls = FALSE` fit (`calc_ses = FALSE`, no GLS
	#    whitening; #307) is the EXPECTED input and is accepted here. The high-dim
	#    branch re-fits its OWN q=1 nuisance internally (#303), so it depends on
	#    neither the input fit's `calc_ses` NOR its `q`: ANY p >= NT fit is accepted,
	#    and the debiased estimate is IDENTICAL across the input fit's q (verified
	#    |diff| = 0, q = 1 vs q < 1) -- the high-dim center uses the re-fit q=1
	#    nuisance, not the input `theta_hat` (which seeds only the `a_theta` identity
	#    check below). A high-dim q >= 1 fit is newly accepted (vs prior versions).
	# A WHITENED (`gls = TRUE`) fixed-p fit with q >= 1: the exact inverse gives the
	# GLS direction, whose validated SE is the oracle `calc_ses = TRUE` path that
	# q >= 1 lacks -- so no validated debiased SE exists, and it is rejected. A
	# `gls = FALSE` / un-whitened fit (`is.na(sig_eps_sq)`) is NOT rejected: it falls
	# through to the exact-inverse OLS + unit-clustered sandwich, valid for any q
	# (#312). `is.na(sig_eps_sq)` is an exact proxy for un-whitened: `gls = FALSE`
	# leaves it NA, while REML estimates it under `gls = TRUE` even when q >= 1.
	if (!highdim && !isTRUE(fit$internal$calc_ses) && !is.na(fit$sig_eps_sq)) {
		stop(
			"debiasedATT() requires a fit computed with q < 1 (bridge selection) in ",
			"the fixed-p (p < NT) regime when the design is GLS-whitened ",
			"(`gls = TRUE`), so a valid debiased standard error exists. This fit has ",
			"calc_ses = FALSE (q >= 1). Re-fit with q < 1 (e.g. q = 0.5), or with ",
			"`gls = FALSE` for the cluster-robust SE.",
			call. = FALSE
		)
	}

	# add_ridge = TRUE fits store a ridge-row-augmented `y_final` (length NT + p)
	# against an unaugmented `X_final` (NT rows), and un-shrink the coefficients
	# post-fit. The debiased construction is not validated for that design, so
	# reject it cleanly rather than letting the dimension mismatch / un-shrink
	# trip the internal identity guard with a cryptic message.
	if (length(y) != n) {
		stop(
			"debiasedATT() does not support fits made with add_ridge = TRUE ",
			"(the ridge-augmented design and post-fit un-shrink are outside the ",
			"validated debiased construction). Re-fit with add_ridge = FALSE.",
			call. = FALSE
		)
	}

	G <- fit$G
	T <- fit$T
	d <- fit$d
	treat_inds <- fit$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs <- fit$cohort_probs

	# The cohort-weight SE channel V2 = `(1/N_tau) sum_g pi_g (catt_g - att)^2` --
	# the variance from estimating the cohort weights pi_hat_g. Prefer the fit's
	# stored plug-in `att_var_2` when present (the same value `fit$att_se` uses,
	# byte-identical on single-sample fits and carrying the correct two-sample
	# formula for `indep_counts` fits; single source of truth, #291 review). It is
	# ABSENT on a `gls = FALSE` fit (`calc_ses = FALSE`, no oracle SE machinery;
	# #307) -- there, recompute the identical plug-in directly via `.plugin_v2()`,
	# which needs only `catt_g` / `att_hat` / the cohort weights, NOT Omega.
	#
	# CAVEAT (#303 review, second-order, deferred -> #295): V2 plugs in the BRIDGE
	# `catt_g` / `att_hat`, while the high-dim center is the q=1 debiased estimate
	# -- the same center(q=1)/variance-channel(bridge) mismatch as the event_study
	# `F_pi` (#309), here for the overall-ATT V2. Second-order (bridge and q=1
	# `catt_g` are both consistent), so the additive SE stays asymptotically
	# correct; the finite-sample effect is part of the #295 coverage validation.
	att_var_2 <- fit$internal$variance_components$att_var_2
	var_weight <- if (
		is.null(att_var_2) ||
			!is.numeric(att_var_2) ||
			length(att_var_2) != 1L ||
			is.na(att_var_2)
	) {
		.plugin_v2(fit)
	} else {
		att_var_2
	}

	# Balanced, unit-major panel is required for the per-unit cluster sums
	# (fetwfe() enforces a balanced panel and sorts rows unit-major).
	if (n %% T != 0L) {
		stop(
			"debiasedATT(): internal design has ",
			n,
			" rows, not a multiple of T = ",
			T,
			" (expected a balanced panel). Please file an issue.",
			call. = FALSE
		)
	}

	# ---- theta-space ATT weight vector a_theta = c(0, A' a_beta) ----
	# Resolve the ACTUAL cohort offsets (#174 / #318): `getFirstInds(G, T)` and
	# `sizes <- (T-1):(T-G)` assume CONSECUTIVE adoption (offsets 2..G+1), which
	# mis-assigns the treatment cells under non-consecutive (scattered-cohort)
	# adoption -- e.g. bacondecomp::divorce -- so the identity guard below would
	# (correctly) fire on a wrong weight vector. The #174 dispatcher reads the real
	# offsets from the fit's cohort labels (and falls back to `getFirstInds()` when
	# the labels are not integer-coercible). Either way it returns exactly
	# `getFirstInds()` for any CONSECUTIVE layout -- whether by deriving the
	# consecutive labels `2..G+1` (e.g. genCoefs/simulateData fixtures) or by the
	# no-label fallback -- so consecutive fits stay byte-identical.
	# `genFullInvFusionTransformMat()` already takes `first_inds`, so the scattered
	# offsets propagate into the transform automatically.
	first_inds <- .resolve_cohort_offsets_and_first_inds(
		fit,
		G = G,
		T = T
	)$first_inds
	# Cohort g's treatment-cell block runs from `first_inds[g]` to the next cohort's
	# first index; a_beta loads treat_inds with weight cohort_probs[g] / (#cells in
	# cohort g) so a_beta' beta = sum_g pi_g * catt_g. For consecutive cohorts the
	# block sizes reduce to (T-1):(T-G).
	sizes <- diff(c(first_inds, num_treats + 1L))
	cohort_of_treat <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- treat_inds[cohort_of_treat == g]
		a_beta[idx] <- cohort_probs[g] / length(idx)
	}
	# Use the SAME fusion transform the fit used (#236 custom matrices included),
	# NOT a hard-coded "cohort" transform -- otherwise event-study / custom fits
	# violate the ATT-direction identity below.
	A <- genFullInvFusionTransformMat(
		first_inds = first_inds,
		T = T,
		G = G,
		d = d,
		num_treats = num_treats,
		fusion_structure = fit$fusion_structure,
		d_inv_treat = fit$internal$d_inv_treat
	)
	a_theta <- c(0, as.numeric(crossprod(A, a_beta)))

	# Invariant: a_theta picks out the same overall ATT the fit reports (the
	# plug-in functional). A failure means the weight vector / transform is
	# inconsistent with the fit.
	if (abs(sum(a_theta * theta_hat) - fit$att_hat) >= 1e-6) {
		stop(
			"debiasedATT(): the ATT-direction identity failed --- the reconstructed ",
			"ATT weight does not match the fit's att_hat. The fit's fusion structure ",
			"may be unsupported. Please file an issue.",
			call. = FALSE
		)
	}

	# ---- debiased estimate: plug-in functional + orthogonal correction ----
	# Two regimes, ONE estimator: only the debiasing direction `v` differs (paper
	# Theorem `debiased.highdim.thm`, "one estimator, two regimes"); everything
	# downstream (the correction, both SE channels) is identical.
	Sig <- crossprod(X) / n
	riesz_diag <- NULL
	if (!highdim) {
		# Fixed-p (p < NT): exact / tiny-ridge inverse. Byte-identical to the
		# validated fixed-p reference; `lambda_c` / `riesz_*` are ignored here.
		v <- solve(Sig + (1e-6 * mean(diag(Sig))) * diag(p), a_theta[-1])
		# Fixed-p nuisance: the bridge theta_hat. Theorem `debiased.att.thm` needs
		# only nuisance *consistency*, and with the exact inverse this reduces to
		# the OLS identity. Byte-unchanged.
		theta_nuis <- theta_hat
	} else {
		# High-dimensional (p >= NT): nodewise (desparsified-lasso) relaxed inverse.
		# `Sig` is singular, so the exact inverse is invalid; `v` is the l1-penalized
		# Riesz representer with theory penalty
		# lambda_node = lambda_c * max(|a|) * sqrt(log p / N), N = clusters = n / T
		# (NOT NT). max(|a|) rescales to the units of the KKT constraint.
		a_th <- a_theta[-1]
		N_units <- n / T
		# Resolve the penalty constant: CV-select per fit (#295, Decision D2) when
		# `lambda_c = "cv"`, else use the supplied fixed constant. The band reuses
		# this same constant (one `lambda_c` for point + band) -- see
		# `.simultaneous_cis_bootstrap()`.
		lambda_cv <- NULL
		lambda_c_used <- lambda_c
		if (identical(lambda_c, "cv")) {
			lambda_cv <- .cv_lambda_node(
				Sig,
				a_th,
				X,
				N_units,
				T,
				riesz_max_iter = riesz_max_iter,
				riesz_tol = riesz_tol
			)
			lambda_c_used <- lambda_cv$lambda_c
		}
		lambda_node <- lambda_node_default(
			p = p,
			N = N_units,
			const = lambda_c_used,
			scale = max(abs(a_th))
		)
		v <- riesz_lasso(
			Sig,
			a_th,
			lambda_node,
			max_iter = riesz_max_iter,
			tol = riesz_tol
		)
		feasibility <- attr(v, "feasibility")
		converged <- attr(v, "converged")
		# The KKT certificate ||Sig v - a||_inf <= lambda_node is exactly the
		# relaxed-inverse feasibility the Theorem 6.6 remainder bound needs.
		if (!.riesz_feasible(feasibility, lambda_node)) {
			warning(
				"debiasedATT(): the high-dimensional nodewise direction did not meet ",
				"its feasibility constraint (||Sig v - a||_inf = ",
				signif(feasibility, 3),
				" > lambda_node = ",
				signif(lambda_node, 3),
				"); the standard error may be unreliable. Raise `riesz_max_iter` ",
				"and/or `lambda_c`.",
				call. = FALSE
			)
		} else if (!converged) {
			warning(
				"debiasedATT(): the high-dimensional nodewise solver hit ",
				"`riesz_max_iter` (the feasibility constraint is met, but more ",
				"iterations would tighten the direction).",
				call. = FALSE
			)
		}
		riesz_diag <- list(
			feasibility = feasibility,
			converged = converged,
			lambda_node = lambda_node,
			lambda_c = lambda_c_used,
			lambda_c_selection = if (identical(lambda_c, "cv")) {
				"cv"
			} else {
				"fixed"
			},
			lambda_cv = if (!is.null(lambda_cv)) {
				lambda_cv[c("cv_loss", "feasible", "mult_grid", "fallback")]
			} else {
				NULL
			}
		)
		# High-dim nuisance (#303): an internal q=1 fused lasso, NOT the reused
		# q<1 bridge theta_hat. Its l1 rate is what Theorem `debiased.highdim.thm`
		# controls the orthogonalization remainder through (the q<1 bridge is
		# super-efficient / non-uniform -- the wrong nuisance for a uniform CI).
		# Deterministic (fixed seed) so the simultaneousCIs() high-dim band center
		# matches debiasedATT()$att exactly.
		theta_nuis <- .fit_q1_nuisance(X, y, N_units, T)
	}
	# Nuisance: the q=1 fused lasso in high-dim (theta_nuis set above), the bridge
	# theta_hat in fixed-p. The plug-in functional and the residuals use it; the
	# debiasing direction `v` and the `a_theta` identity check (above) do not.
	resid <- as.numeric(y - theta_nuis[1] - X %*% theta_nuis[-1])
	score <- as.numeric((X %*% v) * resid)
	att_db <- sum(a_theta * theta_nuis) + mean(score)

	# ---- channel 1: regression (outcome) variance, per-unit clustered ----
	unit <- rep(seq_len(n / T), each = T)
	unit_scores <- tapply(score, unit, sum)
	var_reg <- sum(unit_scores^2) / n^2

	# ---- channel 2: cohort-weight variance ----
	# `var_weight` is the fit's stored `att_var_2` or its `.plugin_v2()`
	# recomputation (read above): the same cohort-weight variance, with the correct
	# single- / two-sample formula.

	se_db <- sqrt(var_reg + var_weight)
	z <- stats::qnorm(1 - alpha / 2)
	out <- list(
		att = att_db,
		se = se_db,
		ci_low = att_db - z * se_db,
		ci_high = att_db + z * se_db,
		var_reg = var_reg,
		var_weight = var_weight
	)
	# High-dimensional fits append the nodewise-direction diagnostics (the fixed-p
	# return is unchanged).
	if (highdim) {
		out$feasibility <- riesz_diag$feasibility
		out$converged <- riesz_diag$converged
		out$lambda_node <- riesz_diag$lambda_node
		out$lambda_c <- riesz_diag$lambda_c
		out$lambda_c_selection <- riesz_diag$lambda_c_selection
		# CV diagnostics (per-draw logging for the #88 coverage gate); NULL for the
		# fixed-`lambda_c` path.
		if (!is.null(riesz_diag$lambda_cv)) {
			out$lambda_cv <- riesz_diag$lambda_cv
		}
	}
	# Lightweight S3 class for a `print` / `tidy` method (#326). The list contents
	# are unchanged -- the regime is inferred from the presence of `lambda_node` --
	# so existing `$`-accessor and `names()` code is unaffected.
	class(out) <- "debiased_att"
	out
}

#' @title Print a debiased overall-ATT estimate
#' @description Compact display of [debiasedATT()]'s estimate, standard error, and
#'   Wald confidence interval. In the high-dimensional (`p >= NT`) regime it also
#'   surfaces the experimental caveat and the nodewise-direction diagnostics
#'   (`lambda_c` / `lambda_node`, KKT feasibility, convergence) the documentation
#'   tells users to inspect; the fixed-`p` print shows only the estimate block.
#' @param x A `"debiased_att"` object from [debiasedATT()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.debiased_att <- function(x, ...) {
	highdim <- !is.null(x$lambda_node)
	cat(sprintf(
		"Debiased overall ATT (%s regime)\n",
		if (highdim) "high-dimensional" else "fixed-p"
	))
	level <- if (is.finite(x$se) && x$se > 0) {
		round(100 * (2 * stats::pnorm((x$ci_high - x$att) / x$se) - 1))
	} else {
		NA_real_
	}
	ci_label <- if (is.na(level)) "Wald CI" else sprintf("%g%% CI", level)
	cat(sprintf("  Estimate:   %.4f\n", x$att))
	cat(sprintf("  Std. Error: %.4f\n", x$se))
	cat(sprintf("  %s: [%.4f, %.4f]\n", ci_label, x$ci_low, x$ci_high))
	if (highdim) {
		.print_highdim_diagnostics(x)
	}
	invisible(x)
}

#' @title Tidy a debiased overall-ATT estimate
#' @description One-row `broom::tidy()` data frame for a [debiasedATT()] result,
#'   following the package's `broom` column convention (`term`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`) -- the same columns `tidy.fetwfe()` /
#'   `tidy.eventStudy()` return.
#' @param x A `"debiased_att"` object from [debiasedATT()].
#' @param ... Ignored.
#' @return A one-row data frame with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`.
#' @importFrom generics tidy
#' @export
tidy.debiased_att <- function(x, ...) {
	data.frame(
		term = "overall ATT",
		estimate = x$att,
		std.error = x$se,
		conf.low = x$ci_low,
		conf.high = x$ci_high,
		stringsAsFactors = FALSE
	)
}

#' Run debiasedATT() on simulated data
#'
#' @description
#' Convenience wrapper mirroring the `*WithSimulatedData()` pattern: fits
#' [fetwfe()] on a `"FETWFE_simulated"` object (via [fetwfeWithSimulatedData()])
#' and returns [debiasedATT()] on the fit.
#'
#' @param simulated_obj An object of class `"FETWFE_simulated"` (from
#'   [simulateData()]).
#' @param q Numeric; the `L_q` bridge penalty exponent. Must be `< 1` (the
#'   debiased SE requires bridge selection). Defaults to `0.5`.
#' @param alpha Numeric in `(0, 1)` or `NULL`; confidence level `1 - alpha`.
#'   Default `NULL`, inheriting the `alpha` stored on the fit.
#' @param ... Further arguments forwarded to [fetwfeWithSimulatedData()] (e.g.
#'   `fusion_structure`, `lambda_selection`, `cv_seed`).
#' @return An object of class `"debiased_att"` (see [debiasedATT()]).
#' @seealso [debiasedATT()], [fetwfeWithSimulatedData()].
#' @examples
#' \dontrun{
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 1)
#'   dat <- simulateData(coefs, N = 200, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
#'   debiasedATTWithSimulatedData(dat)
#' }
#' @export
debiasedATTWithSimulatedData <- function(
	simulated_obj,
	q = 0.5,
	alpha = NULL,
	...
) {
	fit <- fetwfeWithSimulatedData(simulated_obj, q = q, ...)
	debiasedATT(fit, alpha = alpha)
}
