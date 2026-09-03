# Two-tier diagnostic floor for cluster-sandwich quadratic forms.
#
# TWO FAMILIES LIVE HERE, and the asymmetry between them is deliberate.
#
#   * `.floor_cluster_quad()` -- the SCALAR helper (#139, 2024). Its
#     `warning()` / `stop()` raise ORDINARY conditions, deliberately NOT
#     classed. Every one of its call sites is provably outside
#     `.fit_band_for_family()`'s `tryCatch(error = function(e) NULL)` -- a
#     transitive call-graph walk from `.simultaneous_cis_impl()` intersects
#     the SCALAR helper's callers in the empty set -- so there is nothing
#     for a class to key on and classing them would buy nothing (#470).
#     (Deliberately narrower than "the floor-calling functions": the
#     vectorized family below IS reachable from there, which is the whole
#     reason its conditions are classed.)
#
#     THAT CLAIM IS LOAD-BEARING AND IS PINNED, by A11 in
#     `tests/testthat/test-cluster_floor.R` (#476). If it stops holding, a
#     scalar `stop()` raised inside the protected region is swallowed by
#     `error = function(e) NULL` and the band degrades to the POINTWISE one
#     under a `[simultaneous 95% CI]` header -- #470's own defect, at the
#     sites #139 was written for. A11 going red is an instruction to class
#     these two conditions and teach `.fit_band_for_family()` to re-raise
#     them, not to relax the assertion.
#
#     It is also NOT re-expressed on `.floor_psd_diag_core()` below, even
#     though the two share their tiering shape. The reason is the classing
#     asymmetry above, and it stands alone: the core raises CLASSED
#     conditions, so delegating to it would class these two as well and hand
#     `.fit_band_for_family()` a pair it has no site to key on. It would also
#     restate them in the core's wording, which says "on the covariance
#     diagonal" -- true of the vectorized family's argument, false of a
#     scalar quadratic form.
#
#     A COST ARGUMENT USED TO STAND HERE TOO -- that a re-expression "would
#     be a test change in two files for no behavioral gain" -- and it is
#     deleted rather than reworded, because the cost measures at ZERO (#476).
#     A re-expression changing BOTH messages and BOTH condition classes left
#     `test-cluster_floor.R` and `test-namespace-inspect-463.R` green (177
#     assertions before, 177 after) and the whole suite byte-identical: the
#     `test-namespace-inspect-463.R` pin is a round trip that moves with the
#     body, and the unit block greps substrings the new wording keeps. Do not
#     reintroduce a cost claim about those two files without re-measuring.
#
#     The resulting same-file duplication is bounded by a check rather than by
#     intent, but only partly: `test-cluster_floor.R`'s A5b is an EXACT set
#     over every namespace function mentioning `err_threshold` /
#     `warn_threshold`, so a third copy THAT NAMES EITHER THRESHOLD goes
#     red automatically. A third copy with the constants written inline
#     (`if (any(v < -1))`) names neither, has no call site for A1/A4/A9 to
#     collect, and is not in A5c's pin set -- measured: the whole guardrail
#     stays green. That spelling is the more likely accident, not the less.
#
#   * `.floor_cluster_quad_diag()` / `.floor_variance_diag()` -- the
#     VECTORIZED family (#470), sharing `.floor_psd_diag_core()`. Their
#     conditions ARE classed
#     (`fetwfe_negative_variance_floored` /
#     `fetwfe_negative_variance_catastrophic`), because they are the only
#     conditions in this file that cross `.fit_band_for_family()`: that
#     function keys on those two classes to capture both tiers inside its
#     protected region and re-raise them below it. Without the classes the
#     warning is converted by `options(warn = 2)`, swallowed by the
#     `error = function(e) NULL` handler, and the band silently degrades to
#     the POINTWISE one under a `[simultaneous 95% CI]` header -- the #433
#     wrong-answer shape. The classes are therefore part of the interface,
#     not an implementation detail.
#
# Mathematically the cluster-sandwich quadratic forms wrapped by this helper
# (`t(psi) %*% sandwich %*% psi`, sums of outer products) are PSD by
# construction, so any negative value in well-conditioned data is a
# floating-point artifact at machine epsilon (~`1e-15`). Large negatives
# would indicate a bug — a broken PSD invariant, a sign error, or a
# numerical breakdown — and the pre-existing `max(q, 0)` floor would
# silently absorb them, producing `SE = 0` with no signal to the user.
#
# This helper layers a two-tier diagnostic on top of the floor (#139):
#
#   q <  -1     : stop()    -- catastrophic; outside any plausible FP-noise
#                              range or realistic SE^2 magnitude in DiD
#                              applications.
#   q <  -1e-10 : warning() -- clearly outside FP noise; surfaces the
#                              anomaly without breaking the fit.
#   q in [-1e-10, 0] : floor silently to 0 (FP cancellation in well-
#                              conditioned data).
#
# On well-conditioned data the behavior is unchanged from the prior bare
# `max(q, 0)`; the warning/error only surfaces when the quadratic form
# is in a range that genuinely indicates trouble.
#
# Inputs that don't match the contract of a single non-NA numeric scalar
# (length != 1, NA, non-numeric) are passed through unchanged so the
# downstream code can handle them as it would have today.
#
# @param q numeric scalar -- the quadratic form value to floor.
# @param site character scalar -- short site identifier (function name,
#   optionally with `/sub_context` suffix) that gets surfaced in the
#   diagnostic message so a developer can locate the firing site.
# @param err_threshold numeric scalar (default -1); `q < err_threshold`
#   fires `stop()`.
# @param warn_threshold numeric scalar (default -1e-10); `q < warn_threshold`
#   (but `>= err_threshold`) fires `warning()`.
# @return The floored quadratic form: `max(q, 0)` on numeric scalars; `q`
#   unchanged on pass-through inputs.
#
# @keywords internal
#' @noRd
.floor_cluster_quad <- function(
	q,
	site,
	err_threshold = -1,
	warn_threshold = -1e-10
) {
	# Pass-through if the value isn't a single non-NA numeric -- let
	# downstream code handle those cases as it would today.
	if (length(q) != 1L || !is.numeric(q) || is.na(q)) {
		return(q)
	}
	if (q < err_threshold) {
		stop(
			"Cluster-sandwich quadratic form is catastrophically negative (",
			signif(q, 3),
			") at site '",
			site,
			"'; this indicates a bug or severe numerical breakdown. ",
			"Please file an issue at ",
			"https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
	if (q < warn_threshold) {
		warning(
			"Negative cluster-sandwich quadratic form (",
			signif(q, 3),
			") clipped to 0 at site '",
			site,
			"'; if the magnitude is non-trivial this may indicate ",
			"numerical instability or a regression.",
			call. = FALSE
		)
	}
	max(q, 0)
}


# .floor_psd_diag_core
#
# The shared core of the vectorized floor family (#470). Takes the VECTOR of
# covariance-diagonal entries -- never a scalar quadratic form, which is
# `.floor_cluster_quad()`'s job -- floors it at zero, and raises at most ONE
# aggregated condition per call, naming the offending indices and the most
# negative value. One condition per call rather than one per entry: an
# event-study family can carry `K` in the tens, and `K` warnings from a single
# call would bury the signal they exist to raise.
#
# `NA` entries are neither diagnosed nor floored, and this needs no special
# case: `pmax(v, 0)` already propagates `NA` and `which(v < threshold)`
# already drops it, so today's behavior on `NA` is preserved exactly.
#
# Both conditions are classed. See the file header for why these are classed
# and the scalar helper's are not.
#
# @param v numeric vector -- the covariance diagonal to floor. Anything
#   non-numeric is passed through unchanged, matching the scalar helper's
#   spirit.
# @param site character scalar -- short site identifier (function name,
#   optionally with `/sub_context` suffix) surfaced in the diagnostic message
#   so a developer can locate the firing site.
# @param subject character scalar -- the LOWERCASE noun phrase naming what
#   went negative ("cluster-sandwich quadratic form" for the sandwich
#   wrapper, "variance" for the model-based one). The warning tier uses it
#   verbatim; the error tier opens with it, so the core sentence-cases it
#   there. Each wrapper therefore supplies exactly one lowercase phrase.
# @param err_threshold numeric scalar (default -1); any entry
#   `< err_threshold` fires `stop()`.
# @param warn_threshold numeric scalar (default -1e-10); any entry
#   `< warn_threshold`, with none below `err_threshold`, fires `warning()`.
# @return `pmax(v, 0)` on numeric input; `v` unchanged on pass-through.
#
# @keywords internal
#' @noRd
.floor_psd_diag_core <- function(
	v,
	site,
	subject,
	err_threshold = -1,
	warn_threshold = -1e-10
) {
	# Pass-through if the value isn't numeric -- let downstream code handle
	# it as it would today.
	if (!is.numeric(v)) {
		return(v)
	}
	err_idx <- which(v < err_threshold)
	warn_idx <- which(v < warn_threshold)
	if (length(err_idx) > 0L || length(warn_idx) > 0L) {
		fatal <- length(err_idx) > 0L
		idx <- if (fatal) err_idx else warn_idx
		# Cap the rendered index list. `K` can be in the tens, and a caller
		# that captures the condition should not be handed a multi-kilobyte
		# message -- #431 shipped a half-megabyte one exactly that way.
		shown <- if (length(idx) > 10L) {
			paste0(paste(idx[seq_len(10L)], collapse = ", "), ", ...")
		} else {
			paste(idx, collapse = ", ")
		}
		# ONE spelling of the index list in both tiers -- `indices` with a
		# comma-separated list, even when there is a single entry -- so the
		# two messages cannot drift apart.
		where <- paste0(
			" at site '",
			site,
			"': ",
			length(idx),
			" of ",
			length(v),
			" entries "
		)
		worst <- signif(min(v[idx]), 3)
		if (fatal) {
			# Sentence-cased here, inline, because `subject` opens this
			# message and the wrappers supply a lowercase phrase. Inline on
			# purpose, as is the capping above: factoring either into a
			# fourth helper would have to grow `test-cluster_floor.R`'s A5c
			# pin set with it.
			subject_cap <- paste0(
				toupper(substring(subject, 1L, 1L)),
				substring(subject, 2L)
			)
			# The `ci_type = "pointwise"` clause names the escape hatch the
			# NEWS bullet advertises: this tier blocks a fit outright, and an
			# error that blocks a fit should say how to get one. Phrased as
			# "to obtain a fit without the simultaneous band" rather than "pass
			# ci_type = ..." because it must be TRUE ON EVERY ROUTE, and the
			# core serves more than the fit-time one -- on a direct
			# `simultaneousCIs()` call `ci_type` is not that call's remedy,
			# but it is still how the user obtains a usable fit. It is a
			# constant here rather than a per-caller formal: `.floor_psd_diag_core()`'s
			# formals are pinned by name in `test-cluster_floor.R`'s A5c, and a
			# formal added only to vary this sentence would have to grow that
			# pin for no behavioral gain.
			stop(structure(
				list(
					message = paste0(
						subject_cap,
						" on the covariance diagonal is catastrophically ",
						"negative",
						where,
						"below ",
						err_threshold,
						" (indices ",
						shown,
						"; most negative ",
						worst,
						"). This indicates a bug or severe numerical ",
						"breakdown. To obtain a fit without the simultaneous ",
						"band, refit with ci_type = \"pointwise\". ",
						"Please file an issue at ",
						"https://github.com/gregfaletto/fetwfePackage/issues."
					),
					call = NULL
				),
				class = c(
					"fetwfe_negative_variance_catastrophic",
					"error",
					"condition"
				)
			))
		}
		warning(structure(
			list(
				message = paste0(
					"Negative ",
					subject,
					" on the covariance diagonal",
					where,
					"clipped to 0 (indices ",
					shown,
					"; most negative ",
					worst,
					"). If the magnitude is non-trivial this may indicate ",
					"numerical instability or a regression."
				),
				call = NULL
			),
			class = c(
				"fetwfe_negative_variance_floored",
				"warning",
				"condition"
			)
		))
	}
	pmax(v, 0)
}


# .floor_cluster_quad_diag
#
# Matrix in, matrix out with its DIAGONAL floored and diagnosed (#470). The
# K x K cluster-sandwich quadratic form `t(Psi) %*% sandwich %*% Psi` is PSD in
# exact arithmetic, so a negative diagonal entry is either FP cancellation or a
# broken invariant. Off-diagonals are left untouched: they can legitimately
# take either sign.
#
# Takes the quadratic form LITERALLY as its first argument, so the form and its
# floor stay written together in one expression -- the only arrangement
# `tests/testthat/test-cluster_floor.R` can verify, since a form hoisted into
# its own statement and a silent revert are the same shape to a lexical check.
# Exactly two formals, and it names no threshold, so that guardrail's A5a
# arity rule and A5b threshold rule apply to it unchanged.
#
# @param M numeric matrix -- the cluster-sandwich quadratic form. A
#   non-two-dimensional argument is a programming error and raises an ordinary
#   `stop()`; see the guard below for why it is not a pass-through and why the
#   condition is deliberately unclassed.
# @param site character scalar -- short site identifier for the message.
# @return `M` with its diagonal floored at zero.
#
# @keywords internal
#' @noRd
.floor_cluster_quad_diag <- function(M, site) {
	# A guard is needed at all because `diag()` of a length-one numeric BUILDS
	# an identity matrix rather than reading a diagonal -- `diag(5)` is 5 x 5 --
	# so a scalar reaching here would be "floored" into something else entirely.
	#
	# It is a `stop()` rather than a pass-through because nothing legitimately
	# hands this function a non-matrix: its only caller passes
	# `t(Psi_full) %*% sandwich_full %*% Psi_full`. A silent pass-through is
	# strictly worse than either alternative, and #476 measured why: a one-token
	# swap of `.floor_cluster_quad(` to `.floor_cluster_quad_diag(` at the
	# SCALAR `.compute_att_var1()` site -- which hands its floor an
	# `as.numeric()` of the same triple product -- reverted both the #139
	# diagnostic AND the underlying #84-item-9 `max(q, 0)` floor, returning a
	# NEGATIVE `att_var_1` and an understated overall-ATT standard error, with
	# the entire suite byte-identically green. `test-cluster_floor.R`'s A10 now
	# pins which floor function each site calls; this `stop()` is the other half.
	#
	# PLAIN `stop()`, deliberately NOT one of the two classed conditions above:
	# this is a bug report about the caller, not a variance diagnostic, so
	# `.fit_band_for_family()` must neither capture it nor re-raise it below its
	# `tryCatch()`. It degrades to `NULL` there like any other ordinary error.
	#
	# The numeric test lives in `.floor_psd_diag_core()` and not here so the
	# vectorized family makes that decision in exactly one place. Here it would
	# also be redundant: `diag()` of a non-numeric matrix comes back
	# non-numeric, the core passes it through unchanged, and
	# `diag(M) <- diag(M)` rewrites the same values, so a character matrix
	# returns byte-identically either way (pinned in the unit block of
	# `test-cluster_floor.R`). An earlier draft of this comment justified the
	# placement by claiming `diag()` on a `Matrix`-classed argument returns an
	# ordinary numeric vector -- measured FALSE, and deleted: only
	# `Matrix::bdiag` is imported, so `diag` resolves to `base::diag`, which
	# ERRORS on a `ddiMatrix` (`long vectors not supported yet`).
	if (length(dim(M)) != 2L) {
		stop(
			".floor_cluster_quad_diag() requires a two-dimensional matrix ",
			"but was given a ",
			class(M)[1],
			" of length ",
			length(M),
			" at site '",
			site,
			"'; this is a bug in fetwfe, not a property of your data. ",
			"Please file an issue at ",
			"https://github.com/gregfaletto/fetwfePackage/issues.",
			call. = FALSE
		)
	}
	diag(M) <- .floor_psd_diag_core(
		diag(M),
		site,
		"cluster-sandwich quadratic form"
	)
	M
}


# .floor_variance_diag
#
# Vector in, vector out, floored and diagnosed (#470). The sibling of
# `.floor_cluster_quad_diag()` for the MODEL-BASED covariance diagonals in
# `.simultaneous_cis_impl()`, which are not sandwich quantities -- so the
# subject is "variance" and the #139 wording would be a false statement about
# them. The two share `.floor_psd_diag_core()`, and therefore both condition
# classes, so one handler pair in `.fit_band_for_family()` covers the family.
#
# Deliberately absent from `test-cluster_floor.R`'s recognized floor-function
# set: its call sites are not sandwich forms, so that guardrail's A9 ("every
# floor call is handed a real quadratic form") would report
# `.simultaneous_cis_impl` as a decoy. See that file's `WHAT GETS THROUGH`
# list.
#
# @param v numeric vector -- a covariance diagonal.
# @param site character scalar -- short site identifier for the message.
# @return `pmax(v, 0)`.
#
# @keywords internal
#' @noRd
.floor_variance_diag <- function(v, site) {
	.floor_psd_diag_core(v, site, "variance")
}
