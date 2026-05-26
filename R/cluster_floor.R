# Two-tier diagnostic floor for cluster-sandwich quadratic forms.
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
