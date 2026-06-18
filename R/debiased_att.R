# Debiased overall-ATT accessor (#291): a uniformly-valid complement to the
# fused plug-in interval (fit$att_hat / fit$att_se). The fused plug-in SE is the
# within-selection variance and under-covers the *aggregated* overall ATT; the
# debiased estimator restores nominal coverage at ETWFE efficiency. Hand-off from
# the FETWFE methodology paper, Theorem `debiased.att.thm` (eqs
# `debiased.att.def` / `debiased.ols.identity` / `debiased.att.se`). The exact
# algorithm mirrors the validated simulation reference
# `simulations/method_functions.R::debiased_fetwfe` in the paper repo.

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
#'   \item **Correctly specified (GLS-whitened) conditional mean** of the
#'     outcome. Neyman-orthogonality removes first-order sensitivity to the
#'     selected nuisance coefficients, but not to mean-misspecification.
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
#' @param fit A fitted object from [fetwfe()] computed with `q < 1` (so the
#'   bridge selection produces a valid standard error). [etwfe()] / [betwfe()] /
#'   [twfeCovs()] are not supported (the construction is specific to the FETWFE
#'   transformed-coefficient space).
#' @param alpha Numeric in `(0, 1)`; the confidence level is `1 - alpha`.
#'   Defaults to the `alpha` stored on the fit.
#' @param lambda_c Numeric `> 0`; the leading constant of the high-dimensional
#'   (`p >= NT`) nodewise penalty `lambda_node = lambda_c * max(|a|) *
#'   sqrt(log(p) / N)` (`N` = number of units). Larger values shrink the
#'   debiasing direction more (more conservative). **Ignored when `p < NT`.**
#'   Default `1.0`. The coverage-optimal value is regime- and data-dependent and
#'   not yet established by simulation (see *Assumptions*); treat this as a
#'   tunable, experimental knob.
#' @param riesz_max_iter,riesz_tol Integer / numeric; coordinate-descent controls
#'   for the high-dimensional nodewise solver. **Ignored when `p < NT`.**
#'
#' @return A named list with elements:
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
#'   should be `<= lambda_node`), `converged` (the coordinate-descent flag), and
#'   `lambda_node` (the penalty used). These are absent for `p < NT` fits.
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
	if (
		!is.numeric(lambda_c) ||
			length(lambda_c) != 1L ||
			is.na(lambda_c) ||
			lambda_c <= 0
	) {
		stop(
			"debiasedATT(): `lambda_c` must be a single positive number.",
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

	# q >= 1 gate: a valid debiased SE requires the bridge-selection regime
	# (q < 1), which is exactly when the fit computed standard errors. `q` is not
	# stored on the object, so detect via the canonical `calc_ses` flag.
	if (!isTRUE(fit$internal$calc_ses)) {
		stop(
			"debiasedATT() requires a fit computed with q < 1 (bridge selection), ",
			"so a valid debiased standard error exists. This fit has ",
			"calc_ses = FALSE (q >= 1). Re-fit with q < 1 (e.g. q = 0.5).",
			call. = FALSE
		)
	}

	X <- fit$internal$X_final
	y <- as.numeric(fit$internal$y_final)
	theta_hat <- fit$internal$theta_hat
	n <- nrow(X)
	p <- ncol(X)

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

	# Regime is detected from `p` vs `NT` below (see the `v` construction): p < NT
	# uses the exact/ridged inverse; p >= NT uses the nodewise desparsified-lasso
	# direction of paper Theorem `debiased.highdim.thm`. Both share everything
	# else (estimate + SE).

	G <- fit$G
	T <- fit$T
	d <- fit$d
	treat_inds <- fit$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs <- fit$cohort_probs

	# The cohort-weight SE channel is the package's plug-in propensity variance
	# `att_var_2` = `(1/N_T) sum_g pi_g (catt_g - att)^2` -- the variance from
	# estimating the cohort weights pi_hat_g. Reuse the fit's stored value rather
	# than recomputing: it is the same `att_var_2` `fit$att_se` uses, byte-identical
	# on single-sample fits and already carrying the correct two-sample formula for
	# `indep_counts` fits (single source of truth, #291 review).
	#
	# CAVEAT (#303 review, second-order, deferred -> #295): this channel plugs in
	# the BRIDGE `catt_g` / `att_hat`, while the high-dim center is now the q=1
	# debiased estimate -- the same center(q=1)/variance-channel(bridge) mismatch as
	# the event_study propensity channel `F_pi` (#309), here for the overall-ATT V2.
	# It is second-order (the bridge and q=1 `catt_g` are both consistent for the
	# true cohort effects), so the additive SE stays asymptotically correct; the
	# finite-sample effect is part of the #295 high-dim coverage validation.
	var_weight <- fit$internal$variance_components$att_var_2
	if (
		is.null(var_weight) ||
			!is.numeric(var_weight) ||
			length(var_weight) != 1L ||
			is.na(var_weight)
	) {
		stop(
			"debiasedATT(): the fit does not carry a cohort-weight variance ",
			"component (`internal$variance_components$att_var_2`). Re-fit with a ",
			"current version of fetwfe() (>= 1.12.0).",
			call. = FALSE
		)
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
	first_inds <- getFirstInds(G = G, T = T)
	# Cohort g has T - g post-treatment cells; a_beta loads treat_inds with weight
	# cohort_probs[g] / (#cells in cohort g) so a_beta' beta = sum_g pi_g * catt_g.
	sizes <- (T - 1):(T - G)
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
	highdim <- p >= n
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
		lambda_node <- lambda_node_default(
			p = p,
			N = N_units,
			c = lambda_c,
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
		if (feasibility > lambda_node * (1 + riesz_tol)) {
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
			lambda_node = lambda_node
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
	# `var_weight` is the fit's plug-in `att_var_2` (read above): the same
	# cohort-weight variance, with the correct single- / two-sample formula.

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
	}
	out
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
#' @param alpha Numeric in `(0, 1)`; confidence level `1 - alpha`. Default
#'   `0.05`.
#' @param ... Further arguments forwarded to [fetwfeWithSimulatedData()] (e.g.
#'   `fusion_structure`, `lambda_selection`, `cv_seed`).
#' @return The list returned by [debiasedATT()].
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
	alpha = 0.05,
	...
) {
	fit <- fetwfeWithSimulatedData(simulated_obj, q = q, ...)
	debiasedATT(fit, alpha = alpha)
}
