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

#' @importFrom mvtnorm qmvnorm pmvnorm GenzBretz
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
#'   rows give the `K` linear combinations of the underlying per-`(g, t)`
#'   treatment-effect vector (the `multcomp::glht()` convention; `num_treats`
#'   is the number of estimated effects, equal to `G` for `twfeCovs`). Ignored
#'   for the other families. Note that the `"custom"` family omits the
#'   cohort-probability variance term (`Sigma_2 = 0`), so a custom contrast
#'   that pools across cohorts in a probability-weighted way is
#'   anti-conservative (its band can under-cover); use `family = "cohort"` for
#'   cohort-pooled effects.
#' @param method Character; how the simultaneous critical value is computed.
#'   `"analytic"` (default) uses the exact multivariate-normal sup-t quantile via
#'   `mvtnorm::qmvnorm()`. `"bootstrap"` uses the multiplier bootstrap of
#'   Chernozhukov, Chetverikov & Kato (2013): it perturbs the per-unit influence
#'   functions and reads the sup-t critical value off the bootstrap
#'   distribution. The two are asymptotically equivalent; the bootstrap scales
#'   better to large effect families and is heteroskedasticity/cluster-robust by
#'   construction (it always reports cluster-robust per-unit standard errors,
#'   regardless of the fit's `se_type`). `method = "bootstrap"` supports all four
#'   families. For `family = "event_study"` the bootstrap additionally perturbs
#'   the per-unit cohort-probability (propensity) influence function with its own
#'   independent multiplier stream -- a two-channel bootstrap matching the
#'   analytic variance `Sigma = Sigma_1 + Sigma_2` (the event-study family is the
#'   one whose `Sigma_2` is non-zero, because each event-time effect pools across
#'   cohorts weighted by the estimated cohort probabilities). When the (full)
#'   design is **high-dimensional (`p >= NT`)** -- where the analytic Gram inverse
#'   need not exist -- the bootstrap uses the full-design **desparsified**
#'   construction of `debiasedATT()` (per-effect nodewise directions) generalized
#'   to the family. This desparsified `p >= NT` path is **experimental**
#'   (`fetwfe()` fits only; coverage is not yet simulation-validated): inspect the
#'   returned `feasibility` / `converged` diagnostics. A non-`fetwfe()` `p >= NT`
#'   fit (e.g. `betwfe()`) instead falls back to the fixed-`p` selected-support
#'   band (valid when its selected support is low-dimensional). The desparsified
#'   path covers all
#'   four families (a high-dimensional `family = "event_study"` fit additionally
#'   carries the propensity channel `F_pi`). In the high-dimensional regime the
#'   band is centered on the **debiased** estimate (the Theorem 6.6 correction,
#'   equal to `debiasedATT()`'s point estimate for the matching contrast), not
#'   the post-selection bridge estimate; fixed-p bands center on the (unbiased)
#'   bridge estimate as before.
#' @param B Integer; number of multiplier-bootstrap replicates
#'   (`method = "bootstrap"` only). Default `1000`.
#' @param seed Optional integer; if supplied, the bootstrap draws are
#'   reproducible (the ambient random-number stream is saved and restored around
#'   them). If `NULL` (default), the draws come from the ambient generator, so
#'   results vary run to run. Ignored when `method = "analytic"`.
#' @param multiplier Character; the multiplier-weight distribution for the
#'   bootstrap: `"rademacher"` (default, `+/-1`) or `"mammen"` (the Mammen 1993
#'   two-point distribution). Ignored when `method = "analytic"`.
#' @param lambda_c,riesz_max_iter,riesz_tol Controls for the high-dimensional
#'   (`p >= NT`) bootstrap, where each effect's debiasing direction is a nodewise
#'   (desparsified-lasso) `riesz_lasso()` solve. `lambda_c` is the leading
#'   constant of the penalty `lambda_node = lambda_c * max(|a|) * sqrt(log p / N)`:
#'   either a single positive number (the default `1.0` = theory scale) or the
#'   string `"cv"`, which selects the constant **per fit** by cross-validation
#'   (the same CV `debiasedATT()` uses, on the overall-ATT direction, so one
#'   constant serves the point estimate and every band effect; #295). **Smaller
#'   fixed values can leave directions infeasible (a warning fires and those bands
#'   are unreliable).** Ignored when `p < NT` or `method = "analytic"`. The default
#'   stays the fixed `1.0` (the `"cv"` mechanism is opt-in pending coverage
#'   validation).
#' @return An object of S3 class `"simultaneous_cis"`: a list with
#'   \describe{
#'     \item{ci}{A data frame with columns `effect`, `estimate`,
#'       `simultaneous_ci_low`, `simultaneous_ci_high`, `pointwise_ci_low`,
#'       `pointwise_ci_high` (one row per effect in the family).}
#'     \item{adjusted_p_values}{Numeric vector of length `K`: the single-step
#'       max-T multiplicity-adjusted (family-wise) p-value for each effect, the
#'       exact dual of the simultaneous band (a coefficient lies outside the
#'       `(1 - alpha)` band iff its adjusted p-value is `< alpha`). Computed via
#'       `mvtnorm::pmvnorm()` over the same correlation matrix the band uses
#'       (or, under `se_type = "conservative"`, the Bonferroni adjustment
#'       `min(1, K * pointwise_p)`). `NA` for degenerate (zero-variance)
#'       effects. (#200)}
#'     \item{critical_value}{The simultaneous critical value `c_{1 - alpha}`
#'       (or, when the fit used `se_type = "conservative"`, the Bonferroni critical value
#'       `qnorm(1 - alpha/(2K))` -- see Details).}
#'     \item{pointwise_critical_value}{`qnorm(1 - alpha/2)`, for reference.}
#'     \item{bonferroni_critical_value}{`qnorm(1 - alpha/(2K))`, for reference.}
#'     \item{family}{The requested family (character).}
#'     \item{alpha}{The significance level used.}
#'     \item{K}{The number of effects in the family (integer).}
#'   }
#'   For `method = "bootstrap"` the list additionally carries `method`, `B`,
#'   `seed`, `multiplier`, and `regime` (`"fixed-p"` or `"high-dimensional"`); the
#'   `critical_value` is the bootstrap sup-t quantile and the standard errors
#'   backing `ci` are the cluster-robust per-unit ones (so for a non-`cluster` fit
#'   the bootstrap band may differ from the analytic homoskedastic band). A
#'   `"high-dimensional"` fit additionally carries per-effect `feasibility`,
#'   `converged`, and `lambda_node` (the nodewise-direction diagnostics), plus
#'   `lambda_c` (the leading constant used) and `lambda_c_selection` (`"fixed"`
#'   or `"cv"`).
#' @details
#' **Family resolution and `K`.** `"event_study"` resolves to one effect per
#' post-treatment event time `e = 0, ..., T - 2` (`K = T - 1`); `"cohort"` to
#' one effect per treated cohort (`K = G`); `"all_post_treatment"` to one
#' effect per `(g, t)` cell (`K = num_treats`); `"custom"` to the
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
#' **Degenerate (zero-variance) effects.** An effect whose entire contribution
#' to the selected support is zeroed by the bridge penalty -- or, in
#' scattered-cohort panels, an event time with an empty valid-cohort set -- has
#' a standard error of exactly 0 by construction, so its simultaneous and
#' pointwise CIs collapse to a point at the estimate and it is excluded from the
#' joint correlation matrix (it adds no family-wise risk; the critical value is
#' computed over the non-degenerate sub-family). This `se = 0` convention is the
#' simultaneous-CI analog of the `NA` standard error `eventStudy()` reports for
#' the same structurally-degenerate event times; both assign the effect an
#' estimate of 0.
#'
#' **Paper grounding.** Theorem (c') tight Gaussianity (Faletto 2025,
#' `paper_arxiv.tex:1233`) guarantees the joint asymptotic normality; Assumption
#' (Psi-IF) (assumption equation `paper_arxiv.tex:2013`; in-prose discussion at
#' paper line 1268) is the influence-function condition the package's default
#' cohort-sample-proportions estimator satisfies; the fixed-dim framing follows
#' the paper's AE point 1(d).
#'
#' **Conservative fallback.** When the fit was made with `se_type = "conservative"`, the function
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
#'   coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
#'   sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
#'   fit <- fetwfeWithSimulatedData(sim)
#'   sci <- simultaneousCIs(fit, family = "event_study", alpha = 0.05)
#'   print(sci)
#' }
#' @export
simultaneousCIs <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL,
	method = c("analytic", "bootstrap"),
	B = 1000L,
	seed = NULL,
	multiplier = c("rademacher", "mammen"),
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	UseMethod("simultaneousCIs")
}

#' @export
simultaneousCIs.fetwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL,
	method = c("analytic", "bootstrap"),
	B = 1000L,
	seed = NULL,
	multiplier = c("rademacher", "mammen"),
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses,
		method = match.arg(method),
		B = B,
		seed = seed,
		multiplier = match.arg(multiplier),
		lambda_c = lambda_c,
		riesz_max_iter = riesz_max_iter,
		riesz_tol = riesz_tol
	)
}

#' @export
simultaneousCIs.etwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL,
	method = c("analytic", "bootstrap"),
	B = 1000L,
	seed = NULL,
	multiplier = c("rademacher", "mammen"),
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses,
		method = match.arg(method),
		B = B,
		seed = seed,
		multiplier = match.arg(multiplier),
		lambda_c = lambda_c,
		riesz_max_iter = riesz_max_iter,
		riesz_tol = riesz_tol
	)
}

#' @export
simultaneousCIs.betwfe <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL,
	method = c("analytic", "bootstrap"),
	B = 1000L,
	seed = NULL,
	multiplier = c("rademacher", "mammen"),
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses,
		method = match.arg(method),
		B = B,
		seed = seed,
		multiplier = match.arg(multiplier),
		lambda_c = lambda_c,
		riesz_max_iter = riesz_max_iter,
		riesz_tol = riesz_tol
	)
}

#' @note `twfeCovs` is documented as biased in its own roxygen header. The
#'   simultaneous CIs returned here are still mathematically well-defined on
#'   the `twfeCovs` CATT vector under the package's variance machinery, but the
#'   underlying point estimates are not unbiased. Prefer
#'   `simultaneousCIs.fetwfe` / `.etwfe` / `.betwfe` for inference. `twfeCovs`
#'   estimates a single pooled effect per cohort (no per-`(g, t)` cell or
#'   event-time structure), so only `family = "cohort"` and
#'   `family = "custom"` are defined for it; `"event_study"` and
#'   `"all_post_treatment"` raise an error.
#' @export
simultaneousCIs.twfeCovs <- function(
	result,
	family = c("event_study", "cohort", "all_post_treatment", "custom"),
	alpha = 0.05,
	contrasts = NULL,
	method = c("analytic", "bootstrap"),
	B = 1000L,
	seed = NULL,
	multiplier = c("rademacher", "mammen"),
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9
) {
	contract <- .check_for_simultaneous_cis(result)
	family <- match.arg(family)
	.simultaneous_cis_impl(
		x = result,
		family = family,
		alpha = alpha,
		contrasts = contrasts,
		has_valid_ses = contract$has_valid_ses,
		method = match.arg(method),
		B = B,
		seed = seed,
		multiplier = match.arg(multiplier),
		lambda_c = lambda_c,
		riesz_max_iter = riesz_max_iter,
		riesz_tol = riesz_tol
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
	has_valid_ses,
	method = "analytic",
	B = 1000L,
	seed = NULL,
	multiplier = "rademacher",
	lambda_c = 1.0,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9,
	# Severity of the high-dim (p >= NT) all-zero degenerate early-exit: a direct
	# user accessor call warns (default); the silent fit-time band precompute
	# (`.apply_simultaneous_catt_band()`, wrapped in suppressMessages()) passes
	# FALSE so it stays a message() and fitting does not chatter. See #304.
	warn_degenerate_highdim = TRUE
) {
	# --- 1. Argument validation (mvtnorm guard deferred to step 11). ---
	stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
	method <- match.arg(method, c("analytic", "bootstrap"))
	# `lambda_c` is validated unconditionally (symmetric with debiasedATT()). It is
	# only USED on the high-dim bootstrap path, but a malformed value should error
	# regardless of `method` rather than be silently ignored under "analytic".
	if (
		!(identical(lambda_c, "cv") ||
			(is.numeric(lambda_c) &&
				length(lambda_c) == 1L &&
				!is.na(lambda_c) &&
				lambda_c > 0))
	) {
		stop(
			"simultaneousCIs(): `lambda_c` must be a single positive number ",
			"or the string \"cv\".",
			call. = FALSE
		)
	}
	if (method == "bootstrap") {
		if (
			!is.numeric(B) ||
				length(B) != 1L ||
				is.na(B) ||
				B < 1 ||
				B != round(B)
		) {
			stop(
				"simultaneousCIs(): `B` must be a single positive integer.",
				call. = FALSE
			)
		}
		B <- as.integer(B)
		if (
			!is.null(seed) &&
				(!is.numeric(seed) || length(seed) != 1L || is.na(seed))
		) {
			stop(
				"simultaneousCIs(): `seed` must be NULL or a single number.",
				call. = FALSE
			)
		}
	}

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
	G <- x$G
	p <- x$p
	se_type <- if (is.null(x$se_type)) "default" else x$se_type
	is_indep <- isTRUE(x$indep_counts_used)
	tes <- beta_hat[treat_inds]

	# --- 3. Resolve cohort-time offsets + first_inds (same helper eventStudy
	#        uses; preserves scattered-cohort support). ---
	offs <- .resolve_event_study_offsets_and_first_inds(x, G = G, T = T_)
	cohort_offsets_int <- offs$cohort_offsets_int
	first_inds <- offs$first_inds

	# twfeCovs estimates a single pooled treatment effect PER COHORT (so
	# `treat_inds` has length G = num_treats, and each effect is already a
	# cohort ATT), not one effect per (g, t) cell. There is therefore no
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
				"single pooled effect per cohort (no per-(g, t) cell or ",
				"event-time structure). Use family = 'cohort' or ",
				"family = 'custom'.",
				call. = FALSE
			)
		}
		first_inds <- seq_len(G)
	}

	# --- 4. Build the K x num_treats psi_tes matrix (row k = the contrast
	#        that picks effect k out of the (g, t) treatment-effect vector). ---
	psi_tes_mat <- .build_psi_tes_for_family(
		family = family,
		contrasts = contrasts,
		G = G,
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
		G = G,
		T = T_,
		cohort_offsets_int = cohort_offsets_int
	)

	# --- 5. Edge case 1: no valid SEs -> stop(). ---
	if (!isTRUE(has_valid_ses)) {
		stop(
			"simultaneousCIs(): the fitted object was constructed with ",
			"calc_ses = FALSE; standard errors are not available. Re-fit ",
			"with q < 1 (FETWFE/BETWFE) and a satisfied rank condition, or -- for ",
			"a `gls = FALSE` high-dimensional fit -- use `debiasedATT()` for the ",
			"overall-ATT standard error.",
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
		# #236: reuse the custom inverted block stored on the fit (NULL for
		# non-custom fits -> byte-identical `fusion_structure` dispatch).
		d_inv_treat <- .gen_inv_treat_block(
			num_treats = num_treats,
			first_inds = first_inds,
			G = G,
			fusion_structure = x$fusion_structure,
			d_inv_treat = x$internal$d_inv_treat
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
	#        zero estimates / zero SEs / pointwise critical value. In fixed-p this
	#        is consistent with selection consistency (the truth is ~ all-zero), so
	#        a message() suffices. In the high-dimensional (p >= NT) regime
	#        selection is NOT consistent and the debiased band center (generally
	#        non-zero) is bypassed here, so an all-zero band must not be read as
	#        "all effects are zero" -- a direct accessor call warns instead. The
	#        silent fit-time precompute passes `warn_degenerate_highdim = FALSE`
	#        and keeps the message() (swallowed by its suppressMessages()), so
	#        fitting stays quiet. (The cleaner fall-through to the desparsified
	#        path is the deferred #304.) ---
	if (length(sel_treat_inds_shifted) == 0L) {
		if (p >= N * T_ && isTRUE(warn_degenerate_highdim)) {
			warning(
				"simultaneousCIs(): the bridge penalty zeroed out every treatment ",
				"effect, so the returned band is degenerate (all zero). In the ",
				"high-dimensional (p >= NT) regime bridge selection is NOT ",
				"consistent, so this does NOT imply the effects are zero, and the ",
				"debiased high-dimensional band center is bypassed (#304); treat ",
				"this result as unreliable.",
				call. = FALSE
			)
		} else {
			message(
				"simultaneousCIs(): the bridge penalty zeroed out every treatment ",
				"effect; estimates and standard errors are all zero."
			)
		}
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
			adjusted_p_values = rep(NA_real_, K),
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

	# --- 8. Shape Psi (p_sel x K): map each effect's per-cell contrast into the
	#        selected support (theta-space for FETWFE, beta-space for OLS family).
	#        Used by both the analytic Sigma assembly and the bootstrap IF. ---
	Psi <- t(psi_tes_mat %*% d_inv_treat_sel)

	# --- Bootstrap method (#142): build + perturb the per-unit IF matrix instead
	#     of the analytic qmvnorm critical value. Runs BEFORE getGramInv() so the
	#     high-dimensional p >= NT regime (singular Gram, no analytic inverse, so
	#     the nodewise/desparsified direction is used) is reachable; returns early
	#     with the same S3 shape (+ method/B/seed/diagnostics fields). ---
	if (identical(method, "bootstrap")) {
		# High-dimensional full `p >= NT`: build the per-effect full theta-space
		# directions `targets = A' a_beta` (A the inverse fusion transform), the
		# input to the full-design desparsified construction. The desparsified path
		# is FETWFE-only (the only estimator with a regularized `p >= NT` fit); a
		# non-fetwfe `p >= NT` fit (e.g. betwfe) leaves `targets` NULL and falls
		# through to the fixed-p selected-support construction -- valid when its
		# selected support is low-dimensional (the usual sparse case). Fixed-p
		# (`p < NT`) likewise leaves `targets` NULL.
		targets <- NULL
		cell_targets <- NULL
		a_att <- NULL
		if (p >= N * T_ && is_fetwfe) {
			A <- genFullInvFusionTransformMat(
				first_inds = first_inds,
				T = T_,
				G = G,
				d = x$d,
				num_treats = num_treats,
				fusion_structure = x$fusion_structure,
				d_inv_treat = x$internal$d_inv_treat
			)
			a_beta_mat <- matrix(0, p, K)
			a_beta_mat[treat_inds, ] <- t(psi_tes_mat)
			targets <- crossprod(A, a_beta_mat) # p x K full theta-space directions
			# Per-cell full theta-space directions (#309): cell_targets[, j] =
			# A' e_{treat_inds[j]}, the direction whose debiased value is the
			# desparsified per-cohort-time effect tau_db_j. Only the event_study
			# propensity channel needs them (its Sigma_2 must be desparsified too).
			if (identical(family, "event_study")) {
				E_cells <- matrix(0, p, num_treats)
				E_cells[cbind(treat_inds, seq_len(num_treats))] <- 1
				cell_targets <- crossprod(A, E_cells) # p x num_treats
			}
			# Overall-ATT theta-space direction (cohort_probs-weighted), the direction
			# the shared CV penalty constant is selected on (#295 D2: one lambda_c for
			# point + band). For the common consecutive-cohort layout this equals
			# debiasedATT()'s `a_th` exactly (same A / cohort_probs / first_inds), so
			# the CV'd constant -- hence the band center -- matches debiasedATT() under
			# `lambda_c = "cv"`. (For scattered adoption times the two use different
			# `first_inds` -- getFirstInds() in debiasedATT() vs the offset resolver
			# here -- but high-dim debiasedATT() is unsupported there anyway: it errors
			# at its ATT-identity guard, so no shared-constant divergence can ship and
			# the band simply CVs on its own correct direction.)
			a_beta_att <- numeric(p)
			cohort_of_treat <- rep(seq_len(G), times = (T_ - 1):(T_ - G))
			for (g in seq_len(G)) {
				idx <- treat_inds[cohort_of_treat == g]
				a_beta_att[idx] <- x$cohort_probs[g] / length(idx)
			}
			a_att <- as.numeric(crossprod(A, a_beta_att))
		}
		# event_study carries the non-zero cohort-probability (propensity)
		# variance channel Sigma_2: build the SAME per-effect Jacobian list the
		# analytic path uses (Step 10 below) so the bootstrap propensity IF and the
		# analytic Sigma_2 are provably consistent. NULL for the other families
		# (Sigma_2 = 0), which keeps their bootstrap byte-identical to Phase 1/2.
		J_list <- if (identical(family, "event_study")) {
			.build_j_list_for_family(
				family = family,
				K = K,
				G = G,
				T = T_,
				num_treats = num_treats,
				cohort_offsets_int = cohort_offsets_int,
				first_inds = first_inds,
				cohort_probs_overall = cohort_probs_overall,
				d_inv_treat_sel = d_inv_treat_sel
			)
		} else {
			NULL
		}
		# Cell-space cohort-weight Jacobians (#309): the SAME coefficients as J_list
		# but on the identity per-cell map (`d_inv_treat_sel = I`), so the high-dim
		# propensity channel can re-pool the debiased per-cell effects as
		# A_db[, k] = j_cells[[k]] %*% tau_db. Only built for high-dim event_study
		# (cell_targets non-NULL); NULL otherwise -> fixed-p post-selection path.
		j_cells <- if (
			identical(family, "event_study") && !is.null(cell_targets)
		) {
			.build_j_list_for_family(
				family = family,
				K = K,
				G = G,
				T = T_,
				num_treats = num_treats,
				cohort_offsets_int = cohort_offsets_int,
				first_inds = first_inds,
				cohort_probs_overall = cohort_probs_overall,
				d_inv_treat_sel = diag(num_treats)
			)
		} else {
			NULL
		}
		return(.simultaneous_cis_bootstrap(
			family = family,
			alpha = alpha,
			B = B,
			seed = seed,
			multiplier = multiplier,
			lambda_c = lambda_c,
			riesz_max_iter = riesz_max_iter,
			riesz_tol = riesz_tol,
			X_final = X_final,
			y_final = y_final,
			N = N,
			T = T_,
			treat_inds = treat_inds,
			sel_feat_inds = sel_feat_inds,
			sel_treat_inds_shifted = sel_treat_inds_shifted,
			Psi = Psi,
			K = K,
			estimates = estimates,
			effect_labels = effect_labels,
			pointwise_crit = pointwise_crit,
			bonferroni_crit = bonferroni_crit,
			targets = targets,
			a_att = a_att,
			J_list = J_list,
			theta_sel = theta_sel,
			cohort_probs_overall = cohort_probs_overall,
			G = G,
			cell_targets = cell_targets,
			j_cells = j_cells
		))
	}

	# --- 9. (analytic only) getGramInv() + (if cluster) the cluster-robust
	#        sandwich, reusing the eventStudy machinery. ---
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
			"invertible; the analytic method's assumptions are not satisfied. ",
			"For a high-dimensional (p >= NT) design, use method = 'bootstrap'.",
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
		G = G,
		T = T_,
		num_treats = num_treats,
		cohort_offsets_int = cohort_offsets_int,
		first_inds = first_inds,
		cohort_probs_overall = cohort_probs_overall,
		d_inv_treat_sel = d_inv_treat_sel
	)
	Sigma_pi_hat <- .multinomial_cov(cohort_probs_overall[1:G])
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
			# integrate; the per-effect critical value is exact, and the
			# single-step max-T adjusted p-value reduces to the pointwise Wald
			# p-value for the lone non-degenerate effect (NA for degenerate ones).
			crit <- pointwise_crit
			adjusted_p_values <- rep(NA_real_, K)
			nd_one <- which(nondeg)
			if (length(nd_one) == 1L) {
				adjusted_p_values[nd_one] <- 2 *
					stats::pnorm(-abs(estimates[nd_one] / ses[nd_one]))
			}
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
			# preserved by `.with_preserved_rng()`, which save/restores the
			# caller's .Random.seed around a fixed internal set.seed(1L) (the
			# same helper getBetaCV() uses; #195).
			rng_out <- .with_preserved_rng(1L, {
				qmv <- mvtnorm::qmvnorm(
					p = 1 - alpha,
					corr = rho,
					tail = "both.tails",
					algorithm = mvtnorm::GenzBretz()
				)
				# Single-step max-T adjusted p-values over the
				# non-degenerate sub-family -- the exact dual of the qmvnorm
				# band, in the same RNG-protected region.
				list(
					crit = qmv$quantile,
					adjusted_p_values = .maxt_adjusted_p_nd(
						estimates,
						ses,
						nondeg,
						rho
					)
				)
			})
			crit <- rng_out$crit
			adjusted_p_values <- rng_out$adjusted_p_values
		}
	} else {
		# Conservative same-data: the Cauchy-Schwarz upper bound does not
		# generalize to a K x K matrix. Fall back to Bonferroni-corrected
		# pointwise CIs on the conservative scalar SEs.
		v1 <- pmax(diag(Sigma_1), 0)
		v2 <- pmax(diag(Sigma_2), 0)
		ses <- .cauchy_schwarz_se(v1, v2)
		crit <- if (K == 1L) pointwise_crit else bonferroni_crit
		# Bonferroni-adjusted p-values -- the dual of the Bonferroni band --
		# computed from the same Cauchy-Schwarz ses; degenerate (se = 0)
		# effects -> NA.
		pw_cons <- ifelse(
			ses > 0,
			2 * stats::pnorm(-abs(estimates / ses)),
			NA_real_
		)
		adjusted_p_values <- pmin(1, K * pw_cons)
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
		adjusted_p_values = adjusted_p_values,
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

# .maxt_adjusted_p_nd
#' @title Single-step max-T multiplicity-adjusted p-values
#' @description Computes the single-step max-T (family-wise) adjusted p-value for
#'   each non-degenerate effect in a family, the exact dual of the simultaneous
#'   band's `mvtnorm::qmvnorm()` critical value. For effect `k` with standardized
#'   statistic `z_k = estimate_k / se_k`, the adjusted p-value is
#'   `1 - P(max_j |Z_j| <= |z_k|)` for `Z ~ N(0, rho)`, evaluated via
#'   `mvtnorm::pmvnorm()` over the **same** non-degenerate correlation matrix the
#'   band uses (Hothorn, Bretz & Westfall 2008; `multcomp`-style single step).
#' @param estimates Numeric vector of length `K`; the family's point estimates.
#' @param ses Numeric vector of length `K`; the per-effect standard errors
#'   (`sqrt(diag(Sigma))`).
#' @param nondeg Logical vector of length `K`; `TRUE` for non-degenerate
#'   (positive-variance) effects, with `sum(nondeg) >= 2`.
#' @param rho The `sum(nondeg) x sum(nondeg)` correlation matrix
#'   `cov2cor(Sigma[nondeg, nondeg])` (the band's `rho`).
#' @return Numeric vector of length `K`: the adjusted p-value for each
#'   non-degenerate effect, `NA_real_` for degenerate effects.
#' @details Must be called inside the caller's `.Random.seed` save/restore region
#'   (it consumes RNG via `mvtnorm::GenzBretz()`). The result is clamped to
#'   `[0, 1]` to absorb the quasi-Monte-Carlo integration noise.
#' @keywords internal
#' @noRd
.maxt_adjusted_p_nd <- function(estimates, ses, nondeg, rho) {
	K <- length(estimates)
	adj <- rep(NA_real_, K)
	nd_idx <- which(nondeg)
	m <- length(nd_idx)
	z <- abs(estimates[nd_idx] / ses[nd_idx])
	for (j in seq_len(m)) {
		prob <- mvtnorm::pmvnorm(
			lower = rep(-z[j], m),
			upper = rep(z[j], m),
			corr = rho,
			algorithm = mvtnorm::GenzBretz()
		)
		adj[nd_idx[j]] <- min(1, max(0, 1 - as.numeric(prob)))
	}
	adj
}

# .cauchy_schwarz_se
#' @title Conservative (Cauchy-Schwarz) combined standard error
#' @description Combines the two per-effect variance components into the
#'   conservative standard error used by the `se_type = "conservative"` band:
#'   `v1` is the first-order / OLS-sandwich term (`diag(Sigma_1)`) and `v2` is
#'   the cohort-probability term (`diag(Sigma_2)`). With the cross-covariance
#'   left unrestricted, the Cauchy-Schwarz upper bound on
#'   `Var(component_1 + component_2)` is `v1 + v2 + 2 * sqrt(v1 * v2)`, whose
#'   square root equals `sqrt(v1) + sqrt(v2)` (the sum of the two component
#'   SEs). This is the per-effect band analogue of the conservative overall-ATT
#'   SE assembled in `getTeResults*()`.
#' @param v1 Numeric vector (length `K`, entries `>= 0`); the first variance
#'   component per effect (`diag(Sigma_1)`).
#' @param v2 Numeric vector (length `K`, entries `>= 0`); the second variance
#'   component per effect (`diag(Sigma_2)`).
#' @return Numeric vector (length `K`): the conservative combined SE per effect,
#'   `sqrt(v1 + v2 + 2 * sqrt(v1 * v2))` (equivalently `sqrt(v1) + sqrt(v2)`).
#' @keywords internal
#' @noRd
.cauchy_schwarz_se <- function(v1, v2) {
	sqrt(v1 + v2 + 2 * sqrt(v1 * v2))
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
#'   for all four classes; twfeCovs gained an `alpha` slot in #204). Passed
#'   explicitly so all four callers agree.
#' @param has_valid_ses Logical; `res$calc_ses` from the core (read off the
#'   classed object by the caller).
#' @return A list with elements `ci_low`, `ci_high`, and `adjusted_p_values`
#'   (numeric vectors of length `G`, aligned to `x$catt_df` cohort order), or
#'   `NULL` to signal "leave `catt_df` unchanged".
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
				has_valid_ses = TRUE,
				# Fit-time precompute stays silent: keep the high-dim degenerate
				# notice a message() (swallowed by suppressMessages above) rather
				# than a warning(). The user gets the warning by calling
				# simultaneousCIs() / eventStudy() directly. See #304.
				warn_degenerate_highdim = FALSE
			)
		),
		error = function(e) NULL
	)
	if (is.null(sci)) {
		return(NULL)
	}
	# Positional alignment: `.build_psi_tes_for_family(family = "cohort")`
	# iterates g = 1:G in cohort-block order, the SAME order
	# `getCohortATTsFinal()` builds `catt_df`. So `sci$ci` is already
	# row-aligned to `catt_df`. The `effect` labels differ
	# ("Cohort <offset>" vs `catt_df$cohort` = `c_names`), so rely on
	# position, not labels; assert the row count defensively.
	if (nrow(sci$ci) != nrow(x$catt_df)) {
		return(NULL)
	}
	list(
		ci_low = sci$ci$simultaneous_ci_low,
		ci_high = sci$ci$simultaneous_ci_high,
		adjusted_p_values = sci$adjusted_p_values
	)
}

#' @title Apply the fit's `ci_type` to its cohort-family bounds (fit-time)
#' @description Internal finalizer for the `ci_type` default (#197). Given a
#'   fully-classed estimator object, if `ci_type == "simultaneous"` it
#'   overwrites `catt_df`'s `ci_low`/`ci_high` with the cohort-family
#'   simultaneous band (via `.apply_simultaneous_catt_band()`), re-validates,
#'   and returns the updated object. For `ci_type == "pointwise"` (or when the
#'   band degrades to `NULL`) it is a no-op pass-through. Hoisted out of the
#'   four entry-point tails to avoid a 4-site copy (WORKFLOW_LESSONS section 14).
#' @param out A fully-classed `fetwfe`/`etwfe`/`betwfe`/`twfeCovs` object.
#' @param alpha Numeric; the alpha the fit used (read from the fit's `alpha` slot).
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
		# Overwrite the two interval-bound columns AND p_value with the
		# simultaneous-band duals: under ci_type = "simultaneous" the p_value
		# becomes the single-step max-T multiplicity-adjusted p-value that
		# matches the band (#200). se / selected / estimate are untouched.
		# (The catt_df S3 layer only intercepts the OLD Title-Case names;
		# ci_low / ci_high / p_value fall through to `$<-.data.frame`.)
		cd <- out$catt_df
		cd$ci_low <- band$ci_low
		cd$ci_high <- band$ci_high
		cd$p_value <- band$adjusted_p_values
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
#'   `(g, t)` treatment-effect vector (`beta_hat[treat_inds]`) that defines
#'   effect `k`. The point estimate for effect `k` is
#'   `psi_tes_mat[k, ] %*% tes`.
#' @param family,contrasts,G,T,num_treats,cohort_offsets_int,first_inds,cohort_probs_overall
#'   See `.simultaneous_cis_impl()`.
#' @return A numeric `K x num_treats` matrix.
#' @keywords internal
#' @noRd
.build_psi_tes_for_family <- function(
	family,
	contrasts,
	G,
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
		# Effect g = cohort g's ATT = mean(tes over cohort g's block).
		# ORDER-INVARIANT (#258): rows are built in g = 1:G order, matching
		# getCohortATTsFinal()'s catt_df row order, so
		# .apply_simultaneous_catt_band() can positionally row-bind the resulting
		# band to catt_df. Do NOT reorder without making that binding label-based.
		psi_tes <- matrix(0, nrow = G, ncol = num_treats)
		for (g in 1:G) {
			block <- .cohort_block_inds(g, G, first_inds, num_treats)
			psi_tes[g, block] <- 1 / length(block)
		}
		return(psi_tes)
	}

	if (family == "all_post_treatment") {
		# One effect per (g, t) cell: the identity contrast.
		return(diag(num_treats))
	}

	# family == "event_study": one effect per post-treatment event time
	# e = 0, ..., T - 2. Effect e pools the cells (g, e) over the valid cohort
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
#'   `family = "custom"`, `Sigma_2` is set to zero on the assumption that the
#'   user's contrasts operate directly on the fixed per-cell effects; note this
#'   is anti-conservative (not conservative) for a contrast that pools across
#'   cohorts in a probability-weighted way, which does carry cohort-probability
#'   variance.
#' @param family,K,G,T,num_treats,cohort_offsets_int,first_inds,cohort_probs_overall,d_inv_treat_sel
#'   See `.simultaneous_cis_impl()`.
#' @return A length-K list of numeric matrices (each `G x ncol(d_inv_treat_sel)`).
#' @keywords internal
#' @noRd
.build_j_list_for_family <- function(
	family,
	K,
	G,
	T,
	num_treats,
	cohort_offsets_int,
	first_inds,
	cohort_probs_overall,
	d_inv_treat_sel
) {
	p_sel <- ncol(d_inv_treat_sel)
	zero_jac <- matrix(0, nrow = G, ncol = p_sel)

	if (family == "event_study") {
		event_times <- 0:(T - 2L)
		return(lapply(seq_len(K), function(kk) {
			e <- event_times[kk]
			V_e <- which(cohort_offsets_int <= T - e)
			.build_jacobian(
				cohort_probs_overall = cohort_probs_overall,
				G = G,
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
	# single (g, t) cell (all_post_treatment) do not depend on the estimated
	# cohort probabilities; custom contrasts operate on the fixed per-cell
	# effects. (Mirrors getCohortATTsFinal(), whose per-cohort SEs carry no
	# var_2 term.)
	rep(list(zero_jac), K)
}

#' @title Build the `effect` column labels for a family
#' @param family,K,G,T,cohort_offsets_int See `.simultaneous_cis_impl()`.
#' @return A character vector of length `K`.
#' @keywords internal
#' @noRd
.effect_labels_for_family <- function(
	family,
	K,
	G,
	T,
	cohort_offsets_int
) {
	if (family == "event_study") {
		return(paste0("e", 0:(T - 2L)))
	}
	if (family == "cohort") {
		return(paste0("Cohort ", cohort_offsets_int[seq_len(G)]))
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
