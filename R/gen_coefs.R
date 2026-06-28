# Coefficient generators and treatment-effect truth extraction for the
# simulation pipeline: `genCoefs()`, its core (`genCoefsCore()`), and
# `getTes()` / `getActualCohortTes()`. Moved from R/gen_funcs.R in
# 1.9.25.

#' @title Validate the targeted-sparsity arguments of `genCoefs()` (#332)
#' @description Shared, RNG-free validation of `n_signal_cohorts` /
#'   `treat_base_levels` for `genCoefs()` and `genCoefsCore()`. Both `NULL`
#'   (the default) is the uniform-`density` mode and is accepted silently. The
#'   numeric non-degeneracy of the resulting cohort effects (`att != 0`,
#'   `V2 > 0`) is enforced separately by a build-time assertion in
#'   `genCoefsCore()` once `beta` is assembled.
#' @param n_signal_cohorts `NULL`, or a length-1 value coercible to a positive
#'   integer in `1..(G - 1)` (the count of signal cohorts beyond the baseline).
#' @param treat_base_levels `NULL`, or a finite numeric vector of length `G`.
#' @param G Integer; the resolved cohort count.
#' @return `invisible(NULL)`; errors via `stop()` on an invalid spec.
#' @keywords internal
#' @noRd
.validate_targeted_sparsity <- function(
	n_signal_cohorts,
	treat_base_levels,
	G
) {
	if (is.null(n_signal_cohorts) && is.null(treat_base_levels)) {
		return(invisible(NULL)) # uniform-`density` mode
	}
	if (!is.null(n_signal_cohorts) && !is.null(treat_base_levels)) {
		stop(
			"genCoefs(): supply at most one of `n_signal_cohorts` or ",
			"`treat_base_levels` (the count convenience OR the explicit ",
			"per-cohort base-level vector), not both.",
			call. = FALSE
		)
	}
	if (G < 2) {
		stop(
			"genCoefs(): targeted sparsity requires G >= 2 (cross-cohort ",
			"heterogeneity, V2 > 0, needs at least two cohorts).",
			call. = FALSE
		)
	}
	if (!is.null(n_signal_cohorts)) {
		k <- n_signal_cohorts
		if (
			!is.numeric(k) ||
				length(k) != 1L ||
				!is.finite(k) ||
				k != round(k) ||
				k < 1
		) {
			stop(
				"genCoefs(): `n_signal_cohorts` must be a single positive ",
				"integer.",
				call. = FALSE
			)
		}
		if (k > G - 1) {
			stop(
				"genCoefs(): `n_signal_cohorts` (",
				as.integer(k),
				") must be ",
				"<= G - 1 (",
				G - 1,
				"): at least one cohort must stay at the ",
				"baseline for the cohort effects to be heterogeneous (V2 > 0).",
				call. = FALSE
			)
		}
	} else {
		if (
			!is.numeric(treat_base_levels) ||
				length(treat_base_levels) != G ||
				any(!is.finite(treat_base_levels))
		) {
			stop(
				"genCoefs(): `treat_base_levels` must be a finite numeric vector ",
				"of length G (",
				G,
				"), one fused base level per cohort.",
				call. = FALSE
			)
		}
	}
	invisible(NULL)
}

#' Generate Coefficient Vector for Data Generation
#'
#' This function generates a coefficient vector \code{beta} for simulation studies of the fused
#' extended two-way fixed effects estimator. It returns an S3 object of class
#' \code{"FETWFE_coefs"} containing \code{beta} along with simulation parameters \code{G},
#' \code{T}, and \code{d}. See the simulation studies section of Faletto (2025) for details.
#'
#' Optional arguments \code{assignment_type} and \code{assignment_strength}
#' control whether cohort membership in the simulated panel is drawn
#' marginally (independent of the covariates, the original behavior) or from
#' a covariate-dependent propensity-score model --- either a multinomial-logit
#' or an ordered-logit (proportional-odds) model. The default
#' \code{assignment_type = "marginal"} preserves the pre-1.14.0 behavior
#' byte-identically. See \code{vignette("simulation_vignette", package = "fetwfe")}
#' for worked examples.
#'
#' @param G Optional integer. The number of treated cohorts (treatment is assumed to start in periods 2 to
#'   \code{G + 1}). Defaults to \code{NULL}; supply either \code{G} or the
#'   deprecated alias \code{R} (described below).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms
#'   corresponding to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1]. The probability that any given entry in the initial
#'   coefficient vector \code{theta} is nonzero. \code{density = 1} gives a fully dense
#'   (non-sparse) coefficient vector.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each
#'   nonzero entry is set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a
#'   positive value).
#' @param fusion_structure Character. One of \code{"cohort"} (default) or
#'   \code{"event_study"}. Controls the basis in which the true treatment-effect
#'   coefficients are sparse. Under \code{"cohort"} the sparse transformed
#'   coefficients \code{theta} are mapped back through the default two-way
#'   (cohort) inverse fusion transform (byte-identical to previous behavior);
#'   under \code{"event_study"} they are mapped through the event-study inverse
#'   fusion transform, so the true treatment effects are sparse in the
#'   event-study basis (effects sharing the same time since treatment,
#'   \eqn{e = t - g}, tend to be equal across cohorts). This is the
#'   simulation-side companion to the \code{fusion_structure} argument of
#'   \code{fetwfe()}.
#' @param assignment_type Character. One of \code{"marginal"} (default),
#'   \code{"multinomial"}, or \code{"ordered"}. Selects the data-generating
#'   process for cohort assignment in \code{simulateData()}.
#'   \itemize{
#'     \item \code{"marginal"}: each unit's cohort is drawn uniformly,
#'       independent of its covariates. Original pre-1.14.0 behavior,
#'       preserved byte-identically.
#'     \item \code{"multinomial"}: cohort assignment follows a
#'       multinomial-logit propensity-score model
#'       \eqn{\pi_g(x) = \exp(\gamma_g^\top x) / \sum_{g'} \exp(\gamma_{g'}^\top x)}{pi_g(x) = exp(gamma_g' x) / sum_g' exp(gamma_g'' x)}
#'       with \eqn{\gamma_0 \equiv 0}{gamma_0 = 0} (never-treated reference).
#'       Requires \eqn{d \ge 1}{d >= 1}.
#'     \item \code{"ordered"}: cohort assignment follows a proportional-odds
#'       (ordered-logit) model with cumulative probabilities
#'       \eqn{P(W \le g | x) = \mathrm{plogis}(\alpha_g - \gamma^\top x)}{P(W <= g | x) = plogis(alpha_g - gamma' x)}.
#'       Cutpoints \eqn{\alpha_g}{alpha_g} are chosen so the marginal cohort
#'       probabilities are approximately uniform 1/(G + 1). Requires
#'       \eqn{d \ge 1}{d >= 1}.
#'   }
#' @param assignment_strength Non-negative numeric scalar. Scales the logit
#'   coefficients in the propensity-score model. \code{0} reduces both
#'   non-marginal types to the uniform marginal distribution by construction.
#'   Larger values produce stronger covariate-cohort coupling. Defaults to
#'   \code{1.0}. Ignored when \code{assignment_type = "marginal"}.
#' @param assignment_interactions Optional. A list of length-2 integer
#'   vectors, each naming a pair of covariate indices \eqn{(j, k)} in
#'   \eqn{[1, d]} whose elementwise product \eqn{x_{i,j} \cdot x_{i,k}}
#'   enters the propensity model as an additional column. Self-interactions
#'   \code{c(j, j)} are allowed and yield a quadratic term \eqn{x_j^2}.
#'   Unordered pairs are canonicalized to \code{c(min(j, k), max(j, k))}
#'   and duplicates are silently deduplicated (inspect
#'   \code{coefs$assignment_coefs$interactions} on the returned object to
#'   verify the retained list). The interaction columns enter the
#'   propensity model only; the outcome model continues to use the
#'   original covariates. Defaults to \code{NULL} (no interactions ---
#'   v1.14.0 behavior is preserved byte-identically). Passing \code{list()}
#'   (empty list) is treated as equivalent to \code{NULL} --- no interactions
#'   specified, behavior is identical to the v1.14.0 marginal-cohort path
#'   within \code{multinomial} / \code{ordered} DGPs. Errors when
#'   \code{assignment_type = "marginal"} (the marginal DGP has no
#'   propensity model to augment). New in 1.14.1.
#' @param assignment_interaction_strength Optional non-negative numeric
#'   scalar. Scales the Gaussian draws of the interaction coefficients
#'   independently of \code{assignment_strength}. Defaults to \code{NULL},
#'   which means "fall through to \code{assignment_strength}". Useful when
#'   stress-testing whether the nonlinear-propensity angle alone drives
#'   downstream differences (without simultaneously cranking the linear
#'   angle). Ignored when \code{assignment_interactions = NULL}. New in
#'   1.14.1.
#' @param n_signal_cohorts (Optional) Targeted-sparsity mode (#332). A single
#'   positive integer \code{k} in \code{1..(G - 1)}. When supplied, the
#'   coefficient vector is built \emph{deterministically} with the treatment
#'   signal placed on exactly \code{k} cohorts' fused base levels (and the
#'   covariate / interaction blocks left at zero), giving \code{||theta||_0 = k}
#'   non-zeros with a guaranteed non-zero, \emph{heterogeneous} cohort effect
#'   (overall ATT \eqn{\neq 0} and positive cohort-weight variance \eqn{V_2 > 0}).
#'   This produces the simultaneously sparse and non-degenerate
#'   high-dimensional (\eqn{p \ge NT}) data-generating process the desparsified /
#'   nodewise debiased-ATT coverage study needs, which the default uniform
#'   \code{density} placement cannot (it almost never lands signal on the tiny
#'   treatment-effect coordinate block). Cohort 1 is held at the baseline and
#'   cohorts \code{2..(k + 1)} receive the distinct magnitudes
#'   \code{eff_size * (1, 2, ..., k)} (so \code{eff_size} must be non-zero).
#'   Defaults to \code{NULL}; \code{NULL} (with \code{treat_base_levels = NULL})
#'   is the unchanged, byte-identical uniform-\code{density} behavior. Mutually
#'   exclusive with \code{treat_base_levels}.
#' @param treat_base_levels (Optional) Targeted-sparsity mode (#332), explicit
#'   form. A finite numeric vector of length \code{G} giving each cohort's fused
#'   base level directly (placed on \code{getTreatInds()[getFirstInds()]}); every
#'   other coordinate is left at zero, so \code{||theta||_0} is the number of
#'   non-zero entries. Note these are \emph{fused} base levels, not the per-cohort
#'   ATTs: under \code{fusion_structure = "cohort"} a vector
#'   \code{(b_1, ..., b_G)} maps to the cumulative cohort effects
#'   \code{cumsum(b)} (e.g. \code{c(0, m, 0)} gives cohort ATTs \code{(0, m, m)}),
#'   and \code{"event_study"} uses a different but still heterogeneity-preserving
#'   map. The result must be non-degenerate and heterogeneous (\eqn{V_2 > 0}); an
#'   all-zero or constant-effect specification is rejected. Defaults to
#'   \code{NULL}. Mutually exclusive with \code{n_signal_cohorts}.
#' @param seed (Optional) Integer. Seed for reproducibility. Three
#'   deterministic offsets share this seed: the main coefficient draw uses
#'   \code{seed}; the assignment coefficients use \code{seed + 1L}; the
#'   Monte Carlo integration in \code{getTes()} uses \code{seed + 2L}.
#'   \code{NA} (or \code{NULL}) means "draw from the ambient random-number
#'   generator" --- no \code{set.seed()} is called, so any preceding
#'   \code{set.seed()} is respected and the result varies across calls.
#' @param verbose Logical. If \code{TRUE}, emit a \code{message()} when
#'   \code{assignment_interactions} canonicalization removes duplicate or
#'   unordered pairs (e.g., when the user passes both \code{c(1, 2)} and
#'   \code{c(2, 1)}). Default \code{FALSE} (silent --- users can verify the
#'   final canonical list via \code{coefs$assignment_coefs$interactions}).
#' @param R Deprecated. The former name for \code{G}; still accepted with a
#'   deprecation warning, and will be removed in a future release. Use
#'   \code{G}.
#'
#' @return An object of class \code{"FETWFE_coefs"}, which is a list containing:
#' \describe{
#'   \item{beta}{A numeric vector representing the full coefficient vector after the inverse fusion
#'      transform.}
#'   \item{theta}{A numeric vector representing the coefficient vector in the transformed feature
#'		space. \code{theta} is a sparse vector, which aligns with an assumption that deviations from the
#'		restrictions encoded in the FETWFE model are sparse. \code{beta} is derived from
#'		\code{theta}.}
#'   \item{fusion_structure}{The fusion structure (\code{"cohort"} or
#'      \code{"event_study"}) used to build the treatment-effect coefficients.}
#'   \item{G}{The provided number of treated cohorts.}
#'   \item{R}{Deprecated alias for \code{G}, retained for backward
#'      compatibility; populated with the same value. Use \code{G}. Will be
#'      removed in a future release.}
#'   \item{T}{The provided number of time periods.}
#'   \item{d}{The provided number of covariates.}
#'   \item{seed}{The provided seed.}
#'   \item{assignment_type}{The selected cohort-assignment DGP
#'         (\code{"marginal"} / \code{"multinomial"} / \code{"ordered"}).
#'         New in 1.14.0.}
#'   \item{assignment_strength}{The scaling factor applied to the assignment
#'         coefficients. New in 1.14.0.}
#'   \item{assignment_interaction_strength}{The scaling factor applied to
#'         the interaction coefficients. \code{NULL} when no interactions
#'         were specified or when the user passed \code{NULL} (the
#'         fall-through default). New in 1.14.1.}
#'   \item{assignment_coefs}{\code{NULL} when
#'         \code{assignment_type = "marginal"}; otherwise a list with elements
#'         \code{type}, \code{strength}, \code{coefs} (the gamma matrix or
#'         vector), and (for ordered) \code{cutpoints}. Starting in 1.14.1,
#'         \code{assignment_coefs} also carries the sub-slots
#'         \code{interactions} (the canonicalized + deduplicated list of
#'         pairs, or \code{NULL}), \code{delta} (the interaction coefficient
#'         matrix for multinomial or vector for ordered, or \code{NULL}),
#'         and \code{interaction_strength} (the effective scaling factor
#'         used for the \code{delta} draws, or \code{NULL} when no
#'         interactions). New in 1.14.0; \code{interactions}, \code{delta},
#'         and \code{interaction_strength} sub-slots new in 1.14.1.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = G + (T - 1) + d + dG + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = G + (T - 1) + d + dG + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times G - \frac{G(G+1)}{2}}{num_treats = T * G - G(G+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries
#'   occurring with probability \code{density}. Nonzero entries are set to \code{eff_size} or
#'   \code{-eff_size} (with a 60\% chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion
#'   transform to \code{theta} using internal routines:
#'   \code{genBackwardsInvFusionTransformMat()} for the fixed-effect blocks and,
#'   for the treatment-effect block, \code{genInvTwoWayFusionTransformMat()} when
#'   \code{fusion_structure = "cohort"} or \code{genInvEventStudyFusionTransformMat()}
#'   when \code{fusion_structure = "event_study"}.
#' }
#'
#' The multinomial-logit and proportional-odds reference DGPs are the
#' canonical parametric propensity-score models named in Faletto (2025)
#' line 1016; the propensity-weighted population-truth aggregation matches
#' Eq. \code{att.estimator.weighted} (line 837).
#'
#' Note on the random-number generator: passing an explicit numeric
#' \code{seed} calls \code{set.seed(seed)} internally and \strong{leaves the
#' global RNG advanced} after the call returns. This is deliberate --- it
#' makes a simulation both reproducible (the same \code{seed} always yields
#' the same coefficients) and varying (drawing data afterwards consumes the
#' advanced stream). To draw from / preserve the ambient random stream
#' instead --- without calling \code{set.seed()} --- pass \code{seed = NA}
#' (or \code{seed = NULL}).
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @seealso \code{\link{fetwfe}}, whose \code{fusion_structure} argument this
#'   mirrors on the estimation side;
#'   \code{vignette("fusion_structure_vignette", package = "fetwfe")} for the
#'   cohort-vs-event-study distinction, and
#'   \code{vignette("simulation_vignette", package = "fetwfe")} for the full
#'   simulation pipeline.
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients
#'   coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#'   # Simulate data using the coefficients
#'   sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
#'
#'   # Event-study-sparse truth: treatment effects that share the same time
#'   # since treatment are fused across cohorts (the simulation-side companion
#'   # to fetwfe()'s fusion_structure = "event_study").
#'   coefs_es <- genCoefs(
#'     G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
#'     fusion_structure = "event_study", seed = 123
#'   )
#'
#'   # Covariate-dependent cohort assignment: multinomial-logit DGP
#'   coefs_mn <- genCoefs(
#'     G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
#'     assignment_type = "multinomial", assignment_strength = 1.0,
#'     seed = 123
#'   )
#'   sim_mn <- simulateData(coefs_mn, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)
#'
#'   # Covariate-dependent cohort assignment with nonlinear propensity
#'   # (multinomial-logit + a single x1*x2 interaction term in the propensity
#'   # model only; outcome model continues to use plain X):
#'   coefs_int <- genCoefs(
#'     G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
#'     assignment_type = "multinomial",
#'     assignment_interactions = list(c(1, 2)),
#'     assignment_interaction_strength = 1.5,
#'     seed = 123
#'   )
#' }
#'
#' @export
genCoefs <- function(
	G = NULL,
	T,
	d,
	density,
	eff_size,
	fusion_structure = c("cohort", "event_study"),
	assignment_type = c("marginal", "multinomial", "ordered"),
	assignment_strength = 1.0,
	assignment_interactions = NULL,
	assignment_interaction_strength = NULL,
	n_signal_cohorts = NULL,
	treat_base_levels = NULL,
	seed = NULL,
	verbose = FALSE,
	R = NULL
) {
	# Resolve the canonical cohort count (G), mapping the deprecated `R`
	# alias and warning if it is supplied (#41). The body below uses `G`
	# directly.
	G <- .resolve_cohort_count_arg(G, R, "genCoefs")

	# Validate the targeted-sparsity arguments (#332) up front (pure checks,
	# consume no RNG). Both NULL => uniform-`density` mode, byte-identical.
	.validate_targeted_sparsity(n_signal_cohorts, treat_base_levels, G)

	# Check that T is a numeric scalar and at least 2.
	if (!is.numeric(T) || length(T) != 1 || T < 2) {
		stop("T must be a numeric value greater than or equal to 2")
	}

	# Check that G is a numeric scalar and at least 1.
	if (!is.numeric(G) || length(G) != 1 || G < 1) {
		stop(
			"G must be a numeric value greater than or equal to 1 (at least one treated cohort)"
		)
	}

	# Check that G does not exceed T - 1.
	if (G > T - 1) {
		stop("G must be less than or equal to T - 1")
	}

	# Check that d is a numeric scalar and is non-negative.
	if (!is.numeric(d) || length(d) != 1 || d < 0) {
		stop("d must be a non-negative numeric value")
	}

	# Check that density is a numeric scalar in (0, 1] (1 = fully dense, non-sparse).
	if (
		!is.numeric(density) ||
			length(density) != 1 ||
			density <= 0 ||
			density > 1
	) {
		stop("density must be numeric, greater than 0 and at most 1")
	}

	# Check that eff_size is numeric.
	if (!is.numeric(eff_size) || length(eff_size) != 1) {
		stop("eff_size must be a numeric value")
	}

	fusion_structure <- match.arg(fusion_structure)
	assignment_type <- match.arg(assignment_type)
	.validate_strength_arg(assignment_strength, "assignment_strength")
	if (assignment_type != "marginal" && d < 1) {
		stop(sprintf(
			"assignment_type = '%s' requires d >= 1 (at least one covariate)",
			assignment_type
		))
	}

	# Interactions × marginal cross-check
	if (!is.null(assignment_interactions) && assignment_type == "marginal") {
		stop(paste0(
			"assignment_interactions can only be specified when assignment_type ",
			"is 'multinomial' or 'ordered'; the marginal DGP does not use a ",
			"propensity model."
		))
	}

	# Interactions structural validation + canonicalize + dedupe
	if (!is.null(assignment_interactions)) {
		if (!is.list(assignment_interactions)) {
			stop(
				"assignment_interactions must be a list of length-2 integer-pair vectors (or NULL)"
			)
		}
		canon <- vector("list", length(assignment_interactions))
		for (k in seq_along(assignment_interactions)) {
			p <- assignment_interactions[[k]]
			if (!is.numeric(p) || length(p) != 2L) {
				stop(sprintf(
					"assignment_interactions[[%d]] must be a length-2 numeric vector; got length %d",
					k,
					length(p)
				))
			}
			j <- as.integer(p[1L])
			k_idx <- as.integer(p[2L])
			if (anyNA(c(j, k_idx))) {
				stop(sprintf(
					"assignment_interactions[[%d]] coerced to NA; must be integer indices",
					k
				))
			}
			if (j < 1L || j > d || k_idx < 1L || k_idx > d) {
				stop(sprintf(
					"assignment_interactions[[%d]] = c(%d, %d) has an index outside [1, d = %d]",
					k,
					j,
					k_idx,
					d
				))
			}
			canon[[k]] <- c(min(j, k_idx), max(j, k_idx))
		}
		# Dedupe. Silent by default (per Decision Log Q1: users verify via
		# coefs$assignment_coefs$interactions inspection); when verbose =
		# TRUE the canonicalize/dedup step emits a message naming the
		# count of removed pairs.
		deduped <- list()
		for (k in seq_along(canon)) {
			already <- FALSE
			for (existing in deduped) {
				if (identical(canon[[k]], existing)) {
					already <- TRUE
					break
				}
			}
			if (!already) {
				deduped[[length(deduped) + 1L]] <- canon[[k]]
			}
		}
		if (verbose && length(canon) > length(deduped)) {
			message(sprintf(
				"Deduplicated %d assignment_interactions pair(s) after canonicalization.",
				length(canon) - length(deduped)
			))
		}
		assignment_interactions <- deduped
		# Note: assignment_interactions = list() (empty list, e.g. from a
		# user passing list() directly) is treated as equivalent to NULL
		# downstream — see Decision Log RC2 maintainer override.
	}

	# Interaction strength validation
	.validate_strength_arg(
		assignment_interaction_strength,
		"assignment_interaction_strength",
		allow_null = TRUE
	)

	stopifnot(G >= 1)
	stopifnot(T >= 2)
	stopifnot(G <= T - 1)

	core_obj <- genCoefsCore(
		G = G,
		T = T,
		d = d,
		density = density,
		eff_size = eff_size,
		fusion_structure = fusion_structure,
		n_signal_cohorts = n_signal_cohorts,
		treat_base_levels = treat_base_levels,
		seed = seed
	)
	if (is.null(core_obj$beta)) {
		stop(
			"Internal error: genCoefsCore() did not return expected components."
		)
	}

	# Draw assignment coefficients under a deterministic seed offset so
	# the propensity-coefs stream is independent of the main beta draw
	# but reproducible from the same `seed`.
	if (assignment_type == "marginal") {
		assignment_coefs <- NULL
	} else {
		assignment_seed <- if (is.null(seed)) NULL else seed + 1L
		assignment_coefs <- .gen_assignment_coefs(
			G = G,
			d = d,
			type = assignment_type,
			strength = assignment_strength,
			interactions = assignment_interactions,
			interaction_strength = assignment_interaction_strength,
			seed = assignment_seed
		)
	}

	# Create an S3 object of class "FETWFE_coefs"
	obj <- list(
		beta = core_obj$beta,
		theta = core_obj$theta,
		G = G,
		R = G,
		T = T,
		d = d,
		fusion_structure = fusion_structure,
		seed = seed,
		assignment_type = assignment_type,
		assignment_strength = assignment_strength,
		assignment_interaction_strength = assignment_interaction_strength,
		assignment_coefs = assignment_coefs
	)
	class(obj) <- "FETWFE_coefs"
	return(obj)
}


#' Compute True Treatment Effects
#'
#' @description
#' This function extracts the true treatment effects from a full coefficient vector
#' as generated by \code{genCoefs()}. It returns the per-cohort CATTs and an
#' overall ATT. Under the default marginal cohort-assignment DGP, the overall
#' ATT is the equal-weighted mean of the cohort-specific effects. Under the
#' covariate-dependent DGPs introduced in 1.14.0, the overall ATT is a
#' propensity-weighted mean using cohort weights
#' \eqn{E[\pi_g(X)] / \sum_{g' \text{ treated}} E[\pi_{g'}(X)]}{E[pi_g(X)] / sum_{g' treated} E[pi_g'(X)]},
#' matching Faletto (2025) Eq. \code{att.estimator.weighted} (line 837) at the
#' population level. The expected propensities are computed by Monte Carlo
#' integration over the covariate distribution.
#'
#' @param coefs_obj An object of class \code{"FETWFE_coefs"} containing the coefficient vector
#' and simulation parameters.
#'
#' @return An object of class \code{"FETWFE_tes"}, which is a list with the
#' following elements:
#' \describe{
#'   \item{att_true}{A numeric value representing the overall average treatment
#'         effect on the treated. Under the marginal DGP this is the
#'         equal-weighted mean of the cohort-specific effects; under
#'         covariate-dependent DGPs it is the propensity-weighted mean using
#'         \code{cohort_weights}.}
#'   \item{actual_cohort_tes}{A numeric vector of length \code{G} containing the
#'         true cohort-specific treatment effects, calculated by averaging the
#'         coefficients corresponding to the treatment dummies for each cohort.
#'         Intrinsic to \eqn{\beta}{beta}; does not depend on the assignment
#'         DGP.}
#'   \item{cohort_times}{An integer vector of length \code{G} giving the calendar
#'         time period at which each treated cohort first adopts treatment. In
#'         the simulator's convention cohort \code{g} adopts at calendar time
#'         \code{g + 1} (cohort 0 is never-treated).}
#'   \item{cohort_weights}{Numeric vector of length \code{G} summing to 1. Under
#'         the marginal DGP this is uniform \code{1/G}. Under
#'         \code{assignment_type = "multinomial"} or \code{"ordered"} it is
#'         \eqn{E[\pi_g(X)] / \sum_{g' \text{ treated}} E[\pi_{g'}(X)]}{E[pi_g(X)] / sum_{g' treated} E[pi_g'(X)]}.
#'         New in 1.14.0.}
#'   \item{assignment_type}{Character; the cohort-assignment DGP carried over
#'         from \code{coefs_obj} (one of \code{"marginal"},
#'         \code{"multinomial"}, or \code{"ordered"}). Determines whether
#'         \code{cohort_weights} is uniform (marginal) or propensity-weighted.
#'         New in 1.18.1.}
#'   \item{assignment_strength}{Numeric; the assignment-strength scaling carried
#'         over from \code{coefs_obj} (meaningful only when
#'         \code{assignment_type != "marginal"}). \code{NULL} for
#'         \code{FETWFE_coefs} objects saved before 1.14.0. New in 1.18.1.}
#'   \item{G, T, d, seed}{The generating parameters carried over from
#'         \code{coefs_obj} so that \code{print()} and \code{summary()} on the
#'         returned object are self-describing.}
#'   \item{R}{Deprecated alias for \code{G}, retained for backward
#'         compatibility; populated with the same value. Use \code{G}. Will be
#'         removed in a future release.}
#' }
#' Use \code{print()} or \code{summary()} on the returned object for a
#' formatted display.
#'
#' @details
#' The function internally uses auxiliary routines \code{getNumTreats()}, \code{getP()},
#' \code{getFirstInds()}, \code{getTreatInds()}, and \code{getActualCohortTes()} to determine the
#' correct indices of treatment effect coefficients in \code{beta}. The overall treatment effect
#' is computed as a weighted average of the cohort-specific effects (uniform
#' weights under the marginal DGP, propensity weights otherwise).
#'
#' Under non-marginal DGPs, \eqn{E[\pi_g(X)]}{E[pi_g(X)]} is estimated by
#' Monte Carlo integration over the X distribution (Gaussian by default) with
#' \code{M = 10000} draws. The Monte Carlo seed is offset from the main
#' \code{coefs_obj$seed} by \code{+ 2L} per the documented seed-offset
#' convention.
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' \dontrun{
#' # Generate coefficients
#' coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
#'
#' # Compute the true treatment effects:
#' te_results <- getTes(coefs)
#'
#' # Overall average treatment effect on the treated:
#' print(te_results$att_true)
#'
#' # Cohort-specific treatment effects:
#' print(te_results$actual_cohort_tes)
#'
#' # Or use the new print method for a self-describing display:
#' print(te_results)
#'
#' # Propensity-weighted truth under covariate-dependent DGP:
#' coefs_mn <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2,
#'                     assignment_type = "multinomial", assignment_strength = 1.0,
#'                     seed = 42)
#' te_mn <- getTes(coefs_mn)
#' te_mn$att_true       # propensity-weighted overall ATT
#' te_mn$cohort_weights # length G; sums to 1
#' }
#'
#' @export
getTes <- function(coefs_obj) {
	if (!inherits(coefs_obj, "FETWFE_coefs")) {
		stop(
			"coefs_obj must be an object of class 'FETWFE_coefs'; got class: ",
			paste(class(coefs_obj), collapse = ", "),
			call. = FALSE
		)
	}

	# Unpack components from the coefs object
	beta <- coefs_obj$beta
	G <- coefs_obj$G
	T <- coefs_obj$T
	d <- coefs_obj$d

	# Backward-compat: saved-from-v1.13.x FETWFE_coefs objects don't have
	# these slots; coerce a missing $assignment_type to "marginal" so the
	# upgrade path doesn't crash downstream.
	assignment_type <- if (is.null(coefs_obj$assignment_type)) {
		"marginal"
	} else {
		coefs_obj$assignment_type
	}
	# Assignment-strength scaling, carried onto the FETWFE_tes object so
	# print()/summary() can describe the DGP (#189). NULL for FETWFE_coefs
	# objects saved before 1.14.0; only displayed for non-marginal DGPs.
	assignment_strength <- coefs_obj$assignment_strength

	num_treats <- getNumTreats(G = G, T = T)

	p <- getP(G = G, T = T, d = d, num_treats = num_treats)

	stopifnot(length(beta) == p)

	first_inds <- getFirstInds(G = G, T = T)
	treat_inds <- getTreatInds(G = G, T = T, d = d, num_treats = num_treats)

	actual_cohort_tes <- getActualCohortTes(
		G = G,
		first_inds = first_inds,
		treat_inds = treat_inds,
		coefs = beta,
		num_treats = num_treats
	)

	if (assignment_type == "marginal") {
		att_true <- as.numeric(mean(actual_cohort_tes))
		cohort_weights <- rep(1 / G, G)
	} else {
		# Compute E[pi_g(X)] under the propensity-score DGP via MC
		# integration. Reserved-offset seed convention: +2L from main seed.
		# This matches paper Eq. att.estimator.weighted (line 837):
		#   tau_ATT = sum_{g in treated} (N_g / N_tau) * tau_ATT(g)
		# at the population level, where N_g / N_tau converges to
		#   E[pi_g(X)] / sum_{g' treated} E[pi_g'(X)].
		mc_seed <- if (is.null(coefs_obj$seed)) {
			NULL
		} else {
			coefs_obj$seed + 2L
		}
		# Convention A (#254): getTes() is an analysis/truth function, so compute
		# the expected cohort probabilities under the fixed seed `mc_seed` but
		# RESTORE the caller's RNG -- it must not perturb the user's stream. When
		# `mc_seed` is NULL/NA (coefs built without a seed, or with seed = NA) the
		# draw is from the ambient RNG (non-deterministic truth, as before) and is
		# still restored.
		expected_probs <- .with_preserved_rng(
			mc_seed,
			.expected_cohort_probs(
				assignment_coefs = coefs_obj$assignment_coefs,
				d = d,
				distribution = "gaussian",
				M = 10000L,
				seed = NULL
			)
		)
		# expected_probs is length (G + 1); index 1 is never-treated,
		# indices 2..R+1 are treated cohorts.
		treated_probs <- expected_probs[2:(G + 1)]
		cohort_weights <- as.numeric(treated_probs / sum(treated_probs))
		att_true <- as.numeric(sum(cohort_weights * actual_cohort_tes))
	}

	# Cohort adoption times in the simulator's convention: cohort g adopts
	# at calendar time g + 1 (cohort 0 = never-treated, by convention
	# encoded in the panel's `time` integer values 1..T). Stored so
	# downstream tooling (e.g., `tidy.FETWFE_tes`) can label rows with the
	# same scheme that `tidy.<estimator>` uses on a fitted panel.
	cohort_times <- as.integer(seq_len(G) + 1L)

	out <- list(
		att_true = att_true,
		actual_cohort_tes = actual_cohort_tes,
		G = G,
		R = G,
		T = T,
		d = d,
		seed = coefs_obj$seed,
		cohort_times = cohort_times,
		cohort_weights = cohort_weights,
		assignment_type = assignment_type,
		assignment_strength = assignment_strength
	)
	class(out) <- "FETWFE_tes"
	return(out)
}


#' Generate Coefficient Vector for Data Generation (core)
#'
#' This function generates a coefficient vector \code{beta} along with a sparse auxiliary vector
#' \code{theta} for simulation studies of the fused extended two-way fixed effects estimator. The
#' returned \code{beta} is formatted to align with the design matrix created by
#' \code{genRandomData()}, and is a valid input for the \code{beta} argument of that function. The
#' vector \code{theta} is sparse, with nonzero entries occurring with probability \code{density} and
#' scaled by \code{eff_size}. See the simulation studies section of Faletto (2025) for details.
#'
#' @param G Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to
#' \code{G + 1}). Defaults to \code{NULL}; supply either \code{G} or the
#' deprecated alias \code{R} (described below).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms
#' corresponding to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1]. The probability that any given entry in the initial
#' coefficient vector \code{theta} is nonzero. \code{density = 1} gives a fully dense
#' (non-sparse) coefficient vector.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each
#' nonzero entry is set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a
#' positive value).
#' @param fusion_structure Character. One of \code{"cohort"} (default) or
#' \code{"event_study"}. Selects the inverse fusion transform applied to the
#' treatment-effect block, controlling the basis in which the true treatment
#' effects are sparse. \code{"cohort"} is byte-identical to previous behavior;
#' \code{"event_study"} fuses effects at the same time since treatment across
#' cohorts. See \code{genCoefs()} for details.
#' @param n_signal_cohorts,treat_base_levels (Optional) Targeted-sparsity mode
#' (#332). Both \code{NULL} (the default) is the unchanged, byte-identical
#' uniform-\code{density} path. Otherwise the treatment signal is placed
#' deterministically on the per-cohort fused base levels for a sparse,
#' non-degenerate, heterogeneous high-dimensional DGP. See \code{genCoefs()} for
#' the full description; the two are mutually exclusive.
#' @param seed (Optional) Integer. Seed for reproducibility. \code{NA} (or
#' \code{NULL}) means "draw from the ambient random-number generator" --- no
#' \code{set.seed()} is called.
#' @param R Deprecated. The former name for \code{G}; still accepted with a
#' deprecation warning, and will be removed in a future release. Use \code{G}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{beta}}{A numeric vector representing the full coefficient vector after the inverse
#'   fusion transform.}
#'   \item{theta}{A numeric vector representing the coefficient vector in the transformed feature
#'		space. \code{theta} is a sparse vector, which aligns with an assumption that deviations from the
#'		restrictions encoded in the FETWFE model are sparse. \code{beta} is derived from
#'		\code{theta}.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = G + (T - 1) + d + dG + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = G + (T - 1) + d + dG + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times G - \frac{G(G+1)}{2}}{num_treats = T * G - G(G+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries
#' occurring
#'   with probability \code{density}. Nonzero entries are set to \code{eff_size} or \code{-eff_size}
#'   (with a 60\% chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion
#'   transform to \code{theta} using internal routines:
#'   \code{genBackwardsInvFusionTransformMat()} for the fixed-effect blocks and,
#'   for the treatment-effect block, \code{genInvTwoWayFusionTransformMat()} when
#'   \code{fusion_structure = "cohort"} or \code{genInvEventStudyFusionTransformMat()}
#'   when \code{fusion_structure = "event_study"}.
#' }
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' \dontrun{
#'   # Set parameters for the coefficient generation
#'   G <- 3         # Number of treated cohorts
#'   T <- 6         # Total number of time periods
#'   d <- 2         # Number of covariates
#'   density <- 0.1 # Probability that an entry in the initial vector is nonzero
#'   eff_size <- 1.5  # Scaling factor for nonzero coefficients
#'   seed <- 789    # Seed for reproducibility
#'
#'   # Generate coefficients using genCoefsCore()
#'   coefs_core <- genCoefsCore(G = G, T = T, d = d, density = density,
#'   eff_size = eff_size, seed = seed)
#'   beta <- coefs_core$beta
#'   theta <- coefs_core$theta
#'
#'   # For diagnostic purposes, compute the expected length of beta.
#'   # The length p is defined internally as:
#'   #   p = G + (T - 1) + d + d*G + d*(T - 1) + num_treats + num_treats*d,
#'   # where num_treats = T * G - (G*(G+1))/2.
#'   num_treats <- T * G - (G * (G + 1)) / 2
#'   p_expected <- G + (T - 1) + d + d * G + d * (T - 1) + num_treats + num_treats * d
#'
#'   cat("Length of beta:", length(beta), "\nExpected length:", p_expected, "\n")
#' }
#'
#' @export
genCoefsCore <- function(
	G = NULL,
	T,
	d,
	density,
	eff_size,
	fusion_structure = c("cohort", "event_study"),
	n_signal_cohorts = NULL,
	treat_base_levels = NULL,
	seed = NULL,
	R = NULL
) {
	# Resolve the canonical cohort count (G), mapping the deprecated `R`
	# alias and warning if it is supplied (#41). The body below uses `G`
	# directly.
	G <- .resolve_cohort_count_arg(G, R, "genCoefsCore")

	fusion_structure <- match.arg(fusion_structure)

	# Defensive re-validation (genCoefsCore is exported and callable directly);
	# pure checks, no RNG, before `.apply_seed` (#332).
	.validate_targeted_sparsity(n_signal_cohorts, treat_base_levels, G)

	.apply_seed(seed)

	# Check that T is a numeric scalar and at least 2.
	if (!is.numeric(T) || length(T) != 1 || T < 2) {
		stop("T must be a numeric value greater than or equal to 2")
	}

	# Check that G is a numeric scalar and at least 1.
	if (!is.numeric(G) || length(G) != 1 || G < 1) {
		stop(
			"G must be a numeric value greater than or equal to 1 (at least one treated cohort)"
		)
	}

	# Check that G does not exceed T - 1.
	if (G > T - 1) {
		stop("G must be less than or equal to T - 1")
	}

	# Check that d is a numeric scalar and is non-negative.
	if (!is.numeric(d) || length(d) != 1 || d < 0) {
		stop("d must be a non-negative numeric value")
	}

	# Check that density is a numeric scalar in (0, 1] (1 = fully dense, non-sparse).
	if (
		!is.numeric(density) ||
			length(density) != 1 ||
			density <= 0 ||
			density > 1
	) {
		stop("density must be numeric, greater than 0 and at most 1")
	}

	# Check that eff_size is numeric.
	if (!is.numeric(eff_size) || length(eff_size) != 1) {
		stop("eff_size must be a numeric value")
	}

	stopifnot(G >= 1)
	stopifnot(T >= 2)
	stopifnot(G <= T - 1)

	num_treats <- getNumTreats(G = G, T = T)

	p <- getP(G = G, T = T, d = d, num_treats = num_treats)

	theta <- rep(0, p)

	targeted <- !is.null(n_signal_cohorts) || !is.null(treat_base_levels)
	if (!targeted) {
		# ---- Default uniform-`density` path. Kept VERBATIM so that, when the
		#      targeted-sparsity args are unset, the RNG stream (the rbinom at this
		#      line, re-drawn on an all-zero result, then the rfunc sign draw whose
		#      length depends on it) and hence `theta`/`beta` are byte-identical to
		#      pre-#332 behavior. ----
		# Make sure at least one feature is selected
		pass_condition <- FALSE
		while (!pass_condition) {
			theta_inds <- which(as.logical(rbinom(
				n = p,
				size = 1,
				prob = density
			)))
			pass_condition <- length(theta_inds) > 0
		}

		num_coefs <- length(theta_inds)
		# Generate signs of coefficients in transformed space, and bias away from
		# 0.5 (as described in paper)
		signs <- rfunc(num_coefs, prob = 0.6)

		theta[theta_inds] <- eff_size * signs
	} else {
		# ---- Targeted-sparsity path (#332): place a small, heterogeneous signal
		#      on the per-cohort fused BASE levels and leave everything else --
		#      including the covariate/interaction blocks -- at zero, so
		#      `s = ||theta||_0` is exactly the treatment-signal count. FULLY
		#      DETERMINISTIC: no rbinom/rfunc/sample/rnorm here, so this branch
		#      cannot perturb the default RNG stream and is reproducible from its
		#      args alone.
		#
		#      `theta[treat_inds]` are FUSED differences; `treat_inds[first_inds[g]]`
		#      is cohort g's base level. Cohort 1's base (`first_inds[1]`) is the
		#      unique all-ones column of the inverse treatment-fusion transform --
		#      it lifts EVERY cohort's effect equally, contributing to `att` but
		#      0 to the cohort-weight variance V2. The count mode therefore SKIPS
		#      cohort 1 and writes distinct magnitudes on cohorts 2..(k+1), which
		#      load only cohorts {g..G} (a nested staircase) and so produce
		#      heterogeneous cohort effects (V2 > 0). The build-time assertion
		#      after `beta` is assembled enforces the same non-degeneracy for the
		#      free-form `treat_base_levels` vector. See issue #332 (high-dim
		#      p > NT desparsified-ATT coverage study). ----
		ts_first_inds <- getFirstInds(G = G, T = T)
		ts_treat_inds <- getTreatInds(
			G = G,
			T = T,
			d = d,
			num_treats = num_treats
		)
		if (!is.null(n_signal_cohorts)) {
			k <- as.integer(n_signal_cohorts)
			theta[ts_treat_inds[ts_first_inds[2:(k + 1)]]] <- eff_size *
				seq_len(k)
		} else {
			theta[ts_treat_inds[ts_first_inds]] <- treat_base_levels
		}
	}

	# Now we have coefficients that are sparse in the appropriate feature space.
	# The last step is to transform them to the original feature space. Since
	# theta = D %*% beta, beta = solve(D) %*% theta.
	beta <- rep(as.numeric(NA), p)

	beta[1:G] <- genBackwardsInvFusionTransformMat(G) %*% theta[1:G]

	stopifnot(all(is.na(beta[(G + 1):(G + T - 1)])))
	beta[(G + 1):(G + T - 1)] <- genBackwardsInvFusionTransformMat(T - 1) %*%
		theta[(G + 1):(G + T - 1)]

	if (d > 0) {
		# Coefficients corresponding to X don't need to be transformed
		stopifnot(all(is.na(beta[(G + T - 1 + 1):(G + T - 1 + d)])))
		beta[(G + T - 1 + 1):(G + T - 1 + d)] <- theta[
			(G + T - 1 + 1):(G + T - 1 + d)
		]

		# Cohort-X interactions (one cohort at a time, with all interactions for
		# X. So G blocks of size d.)

		for (j in 1:d) {
			first_ind_j <- G + T - 1 + d + j
			last_ind_j <- G + T - 1 + d + (G - 1) * d + j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)

			stopifnot(length(inds_j) == G)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- genBackwardsInvFusionTransformMat(G) %*%
				theta[inds_j]
		}

		stopifnot(all(!is.na(beta[1:(G + T - 1 + d + G * d)])))
		stopifnot(all(is.na(beta[(G + T - 1 + d + G * d + 1):p])))

		# Time-X interactions
		for (j in 1:d) {
			first_ind_j <- G + T - 1 + d + G * d + j
			last_ind_j <- G + T - 1 + d + G * d + (T - 2) * d + j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)
			stopifnot(length(inds_j) == T - 1)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- genBackwardsInvFusionTransformMat(T - 1) %*%
				theta[inds_j]
		}

		stopifnot(all(!is.na(beta[1:(G + T - 1 + d + G * d + (T - 1) * d)])))
		stopifnot(all(is.na(beta[(G + T - 1 + d + G * d + (T - 1) * d + 1):p])))
	}

	# Base treatment effects: need to identify indices of first treatment
	# effect for each cohort
	first_inds <- getFirstInds(G = G, T = T)

	treat_inds <- getTreatInds(G = G, T = T, d = d, num_treats = num_treats)

	stopifnot(all(is.na(beta[treat_inds])))

	beta[treat_inds] <- .gen_inv_treat_block(
		num_treats,
		first_inds,
		G,
		fusion_structure
	) %*%
		theta[treat_inds]

	stopifnot(all(
		!is.na(beta[1:(G + T - 1 + d + G * d + (T - 1) * d + num_treats)])
	))

	if (d > 0) {
		stopifnot(all(is.na(beta[
			(G + T - 1 + d + G * d + (T - 1) * d + num_treats + 1):p
		])))

		# Treatment effect-X interactions
		for (j in 1:d) {
			first_ind_j <- G + T - 1 + d + G * d + (T - 1) * d + num_treats + j
			last_ind_j <- G +
				T -
				1 +
				d +
				G * d +
				(T - 1) * d +
				num_treats +
				(num_treats - 1) * d +
				j

			inds_j <- seq(first_ind_j, last_ind_j, by = d)

			stopifnot(length(inds_j) == num_treats)
			stopifnot(all(is.na(beta[inds_j])))

			beta[inds_j] <- .gen_inv_treat_block(
				num_treats,
				first_inds,
				G,
				fusion_structure
			) %*%
				theta[inds_j]
		}
	}

	stopifnot(all(!is.na(beta)))

	# Targeted-sparsity non-degeneracy contract (#332): in targeted mode the point
	# of the DGP is a NON-degenerate, HETEROGENEOUS cohort effect. The two guards
	# below catch the two DEGENERATE failure modes on the assembled `beta`
	# (intrinsic to it, DGP-weight-independent): all-zero cohort effects (=> overall
	# att = 0) and all-equal cohort effects (=> V2 = 0). Between them they reject
	# `eff_size = 0`, an all-zero `treat_base_levels`, and the homogeneous trap
	# (e.g. signal only on cohort 1's all-ones base, or any constant
	# `treat_base_levels`) -- for either fusion_structure -- without the caller
	# having to reason about the fused map. Note they do NOT prove `att != 0` for an
	# arbitrary mixed-sign free-form `treat_base_levels` (equal-and-opposite levels
	# could cancel to att ~ 0 while keeping V2 > 0); the roxygen only promises
	# V2 > 0 for the explicit vector, which is the caller's responsibility on that
	# escape hatch. (The count mode cannot reach either failure given the validated
	# 1 <= k <= G-1, but the check is cheap insurance.)
	if (targeted) {
		catt_built <- getActualCohortTes(
			G = G,
			first_inds = first_inds,
			treat_inds = treat_inds,
			coefs = beta,
			num_treats = num_treats
		)
		if (all(catt_built == 0)) {
			stop(
				"genCoefs() targeted mode produced all-zero cohort effects ",
				"(att == 0). Supply a nonzero `eff_size` (count mode) or a ",
				"`treat_base_levels` with a nonzero entry.",
				call. = FALSE
			)
		}
		if (length(unique(catt_built)) == 1L) {
			stop(
				"genCoefs() targeted mode produced homogeneous cohort effects ",
				"(V2 == 0): all cohorts have the same effect, so there is no ",
				"cohort-weight variance. Place signal on at least one non-baseline ",
				"cohort with a distinct level (count mode does this automatically; ",
				"for `treat_base_levels` use entries that are not all equal).",
				call. = FALSE
			)
		}
	}

	# Confirm beta satisfies input requirements of genRandomData() (make
	# up values for N, sig_eps_sq, and sig_eps_c_sq that meet requirements)

	testGenRandomDataInputs(
		beta = beta,
		G = G,
		T = T,
		d = d,
		N = G + 1,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	return(list(beta = beta, theta = theta))
}


#' Calculate True Cohort Average Treatment Effects from Coefficients
#'
#' Given a full coefficient vector and information about treatment effect indices,
#' this function calculates the true average treatment effect for each cohort by
#' averaging the relevant treatment effect coefficients.
#'
#' @param G Integer. Number of treated cohorts.
#' @param first_inds Integer vector. `first_inds[g]` is the index (within the
#'   block of treatment effect coefficients) of the first treatment effect for cohort `g`.
#' @param treat_inds Integer vector. Indices in the full `coefs` vector that
#'   correspond to the block of all treatment effect coefficients.
#' @param coefs Numeric vector. The full true coefficient vector \eqn{\beta}.
#' @param num_treats Integer. Total number of treatment effect parameters.
#'
#' @return A numeric vector of length G, where the g-th element is the true
#'   average treatment effect for cohort g.
#' @keywords internal
#' @noRd
getActualCohortTes <- function(G, first_inds, treat_inds, coefs, num_treats) {
	actual_cohort_tes <- rep(as.numeric(NA), G)

	for (g in 1:G) {
		inds_g <- .cohort_block_inds(g, G, first_inds, num_treats)
		actual_cohort_tes[g] <- mean(coefs[treat_inds][inds_g])
	}

	return(actual_cohort_tes)
}
