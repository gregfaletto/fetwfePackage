# High-dimensional (p >= NT) debiasing direction for debiasedATT() (#31).
#
# The fixed-p debiased path builds the debiasing direction as the exact/ridged
# inverse v = (Sigma_hat + tiny*I)^{-1} a, which is invalid once p >= NT
# (Sigma_hat = X'X/n is singular and the tiny ridge yields an unstable v). In the
# high-dimensional regime v is instead the nodewise (desparsified) Riesz
# representer of the high-dimensional FETWFE theory, computed as the
# l1-penalized solve below. Ported verbatim from the validated reference
# `simulations/highdim_functions.R::riesz_lasso` in the FETWFE paper repo.

#' Nodewise (desparsified) debiasing direction for the high-dimensional regime
#'
#' @description
#' Solves `v_hat = argmin_v (1/2) v' Sig v - a' v + lambda * ||v||_1` by
#' coordinate descent. By the KKT conditions the minimizer satisfies
#' `|| Sig %*% v_hat - a ||_inf <= lambda` --- the relaxed-inverse feasibility the
#' high-dimensional FETWFE theory remainder bound uses --- and as `lambda -> 0`
#' with `Sig` nonsingular it returns the exact inverse `Sig^{-1} a` (the fixed-p
#' boundary). `O(p^2)` per sweep via a rank-1 maintenance of `Sig %*% v`.
#'
#' @param Sig Numeric matrix; the (whitened theta-space) Gram `X'X / n`.
#' @param a Numeric vector; the target direction (length `ncol(Sig)`).
#' @param lambda Numeric scalar; the nodewise l1 penalty.
#' @param max_iter Integer; maximum coordinate-descent sweeps.
#' @param tol Numeric; convergence tolerance on the max coordinate change.
#' @return `v` with attributes `feasibility` (`= ||Sig v - a||_inf`), `iters`,
#'   and `converged` (`FALSE` => `lambda` likely too small / more iterations
#'   needed).
#' @keywords internal
#' @noRd
riesz_lasso <- function(Sig, a, lambda, max_iter = 5000L, tol = 1e-9) {
	p <- length(a)
	v <- numeric(p)
	Sv <- numeric(p) # Sig %*% v, maintained incrementally
	dgSig <- diag(Sig)
	it <- 0L
	converged <- FALSE
	for (it in seq_len(max_iter)) {
		max_change <- 0
		for (j in seq_len(p)) {
			if (dgSig[j] <= 0) {
				# All-zero (or non-positive-diagonal) Gram column: leave v_j = 0.
				# Its residual |Sv - a|_j = |a_j| still enters `feasibility` below --
				# see the note there; this is intentional, not an oversight.
				next
			}
			cj <- Sv[j] - dgSig[j] * v[j] # (Sig v)_j without the j-th term
			rho <- a[j] - cj
			vj_new <- sign(rho) * max(0, abs(rho) - lambda) / dgSig[j] # soft-threshold
			dv <- vj_new - v[j]
			if (dv != 0) {
				Sv <- Sv + Sig[, j] * dv # rank-1 update of Sig %*% v
				v[j] <- vj_new
				adv <- abs(dv)
				if (adv > max_change) {
					max_change <- adv
				}
			}
		}
		if (max_change < tol) {
			converged <- TRUE
			break
		}
	}
	# Feasibility certificate ||Sig v - a||_inf, the KKT residual. Semantics note
	# (#325): the max is taken over ALL p coordinates, INCLUDING any skipped
	# zero-diagonal columns (whose residual is the full |a_j|). This is deliberate:
	# a zero Gram column with a nonzero target means the Riesz representer genuinely
	# does not exist in that coordinate, so the constraint truly cannot be met and
	# reporting infeasibility is correct. Do NOT "fix" this by restricting the max
	# to active (`dgSig > 0`) coordinates -- that would mask a real infeasibility and
	# let an undefined direction pass the gate. Unreachable from the real pipeline
	# (the full design's `min(diag(Sig))` is bounded well above 0), so this guards a
	# theoretical corner, not an observed one.
	attr(v, "feasibility") <- max(abs(Sv - a)) # ||Sig v - a||_inf
	attr(v, "iters") <- it
	attr(v, "converged") <- converged # FALSE => lambda likely too small
	v
}

#' Theory-scaled nodewise penalty `lambda_node = c * scale * sqrt(log p / N)`
#'
#' @description The high-dimensional FETWFE theory rate scale. The effective
#'   sample size is the number of clusters `N` (units), not `n = NT`. `scale`
#'   rescales to the units of the constraint `||Sig v - a||_inf` (the accessor
#'   passes `scale = max(abs(a))`).
#' @keywords internal
#' @noRd
lambda_node_default <- function(p, N, const = 1.0, scale = 1.0) {
	const * scale * sqrt(log(p) / N)
}

# Relative tolerance for the high-dim nodewise KKT feasibility certificate. The
# desparsified L1 direction BINDS its constraint by construction
# (||Sig v - a||_inf = lambda_node at the L1 optimum), so a converged solver
# overshoots only by coordinate-descent round-off (~1e-9..1e-8 relative). The
# feasibility slack must therefore floor ABOVE solver precision -- it is
# deliberately DECOUPLED from the convergence control `riesz_tol`, which
# previously doubled as the slack at 1e-9 (below precision), spuriously flagging
# legitimately-binding directions as infeasible (gregfaletto/fetwfe#88). A
# genuine infeasibility (lambda_c too small / non-convergence) overshoots by
# >> 1e-6, so 1e-6 cleanly separates round-off from a real constraint violation.
.RIESZ_FEASIBILITY_RTOL <- 1e-6

#' @title KKT feasibility test for a high-dim nodewise direction
#' @description `TRUE` iff the relaxed-inverse certificate holds to within the
#'   feasibility tolerance: `feasibility <= lambda_node * (1 + rtol)`. Vectorized.
#'   The single source for every nodewise feasibility check -- the `debiasedATT()`
#'   / `simultaneousCIs()` warnings, the `lambda_c = "cv"` grid gate, and the
#'   high-dim print diagnostics -- so they always agree. See `.RIESZ_FEASIBILITY_RTOL`.
#' @param feasibility Numeric (scalar or vector); `||Sig v - a||_inf`.
#' @param lambda_node Numeric (scalar or vector); the nodewise penalty scale.
#' @param rtol Relative feasibility tolerance; defaults to `.RIESZ_FEASIBILITY_RTOL`.
#' @return Logical of the recycled length of `feasibility` / `lambda_node`.
#' @keywords internal
#' @noRd
.riesz_feasible <- function(
	feasibility,
	lambda_node,
	rtol = .RIESZ_FEASIBILITY_RTOL
) {
	feasibility <= lambda_node * (1 + rtol)
}

#' Cross-validate the high-dimensional nodewise penalty constant `lambda_c`
#'
#' @description
#' Selects the leading constant `lambda_c` of the nodewise penalty
#' `lambda_node = lambda_c * max(|a|) * sqrt(log p / N)` **per fit** by
#' cross-validating the unit-level Riesz loss `L(v) = 0.5 v' Sig v - a' v`,
#' restricted to the **KKT-feasible region** (the grid points where
#' `riesz_lasso()` converges and `||Sig v - a||_inf <= lambda_node`). This is the
#' desparsified-lasso / auto-DML penalty-selection standard (van de Geer; hdi;
#' Chernozhukov-Newey-Singh), replacing the fixed theory-scale `lambda_c = 1.0`
#' (#295, Decision D2). The theory scale `1.0` is the documented **fallback** when
#' no grid point is feasible.
#'
#' The constant is in pure-constant units (the `max(|a|)` constraint scale and the
#' `sqrt(log p / N)` rate are factored out by `lambda_node_default()`), so the
#' same selected constant transfers across effect directions -- the D2 "one
#' `lambda_c` for the point estimate AND the simultaneous band" requirement
#' (`debiasedATT()` CV-selects on the overall-ATT direction; the band reuses the
#' constant, scaling each effect's `lambda_node_k` by its own `max(|a_k|)`).
#'
#' **Determinism.** The only RNG draw is the unit-fold assignment, seeded from the
#' data (`as.integer(N * T)`, the `.fit_q1_nuisance()` convention) under
#' `.with_preserved_rng()`; `riesz_lasso()` is itself deterministic. So the CV is
#' reproducible and identical at the `debiasedATT()` and `simultaneousCIs()` call
#' sites (same `Sig` / `X` / `N` / `T` / direction).
#'
#' @param Sig_full Numeric `p x p`; the full-data Gram `crossprod(X) / n`.
#' @param a Numeric length `p`; the target (theta-space) direction.
#' @param X Numeric `n x p` (`n = N*T`); the unit-major balanced design (for the
#'   fold-specific Grams).
#' @param N_units,T Integer; units and time periods.
#' @param mult_grid Numeric; multipliers of the theory scale to search (the grid
#'   spans below and above `1.0` so the CV sees the feasibility edge).
#' @param nfolds Integer; CV folds over units.
#' @param riesz_max_iter,riesz_tol Passed to `riesz_lasso()` for the fold solves.
#' @param gate_max_iter Integer; iteration cap for the full-data feasibility gate
#'   ONLY (the gate just classifies feasible/infeasible, so a smaller cap is
#'   cheap); floored at `riesz_max_iter` if that is smaller.
#' @return A list: `lambda_c` (the selected constant, or `1.0` on fallback),
#'   `fallback` (logical), `cv_loss` / `feasible` (length-`mult_grid`), `mult_grid`.
#' @keywords internal
#' @noRd
.cv_lambda_node <- function(
	Sig_full,
	a,
	X,
	N_units,
	T,
	mult_grid = c(0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0),
	nfolds = 5L,
	riesz_max_iter = 5000L,
	riesz_tol = 1e-9,
	gate_max_iter = 500L
) {
	p <- length(a)
	# Theory anchor; the grid is `mult_grid * lam0` (pure-constant units).
	lam0 <- lambda_node_default(
		p = p,
		N = N_units,
		const = 1.0,
		scale = max(abs(a))
	)
	gate_iter <- min(gate_max_iter, riesz_max_iter)

	# --- Feasibility gate (full data): only KKT-feasible grid points compete. ---
	# Deliberately conservative: a grid point qualifies only if the reduced-budget
	# (`gate_iter`) solve BOTH converges AND certifies `||Sig v - a||_inf <= lambda`.
	# The deployed solve (debiasedATT() / .build_regression_if_highdim()) runs the
	# full `riesz_max_iter` budget and accepts feasibility regardless of `converged`,
	# so the gate can in principle exclude a lambda that becomes feasible only after
	# more iterations, biasing selection toward larger (better-conditioned) constants.
	# That is the safe direction for a pre-filter (CV competes only among confidently
	# feasible candidates) and was not observed to bite; revisit under the #88
	# coverage work if the feasible edge must extend to smaller lambda.
	lambdas <- mult_grid * lam0
	feasible <- vapply(
		lambdas,
		function(lam) {
			vf <- riesz_lasso(
				Sig_full,
				a,
				lam,
				max_iter = gate_iter,
				tol = riesz_tol
			)
			isTRUE(attr(vf, "converged")) &&
				.riesz_feasible(attr(vf, "feasibility"), lam)
		},
		logical(1)
	)
	feas_idx <- which(feasible)
	if (length(feas_idx) == 0L) {
		# No feasible grid point -> documented theory-scale fallback.
		return(list(
			lambda_c = 1.0,
			fallback = TRUE,
			cv_loss = rep(NA_real_, length(mult_grid)),
			feasible = feasible,
			mult_grid = mult_grid
		))
	}

	# --- Unit-fold CV of the Riesz loss over the feasible grid points. ---
	# Data-derived fold seed (the only RNG draw; RNG state preserved). Unit-major
	# balanced panel: unit i occupies rows ((i-1)*T + 1):(i*T). Two configs sharing
	# `N*T` share this fold seed, which is harmless: the seed only fixes a partition
	# of units into folds, so reuse just means an identical (still valid) split.
	cv_seed <- as.integer(min(as.numeric(N_units) * T, .Machine$integer.max))
	fold <- .with_preserved_rng(
		cv_seed,
		sample(rep(seq_len(nfolds), length.out = N_units))
	)
	unit_of_row <- rep(seq_len(N_units), each = T)
	cv_loss <- rep(NA_real_, length(mult_grid))
	for (gi in feas_idx) {
		lam <- lambdas[gi]
		fold_losses <- vapply(
			seq_len(nfolds),
			function(f) {
				test_rows <- unit_of_row %in% which(fold == f)
				X_tr <- X[!test_rows, , drop = FALSE]
				X_te <- X[test_rows, , drop = FALSE]
				Sig_tr <- crossprod(X_tr) / nrow(X_tr)
				Sig_te <- crossprod(X_te) / nrow(X_te)
				v_tr <- riesz_lasso(
					Sig_tr,
					a,
					lam,
					max_iter = riesz_max_iter,
					tol = riesz_tol
				)
				# Held-out Riesz loss 0.5 v' Sig_te v - a' v.
				0.5 *
					as.numeric(crossprod(v_tr, Sig_te %*% v_tr)) -
					sum(a * v_tr)
			},
			numeric(1)
		)
		cv_loss[gi] <- mean(fold_losses)
	}
	best <- feas_idx[which.min(cv_loss[feas_idx])]
	list(
		lambda_c = mult_grid[best],
		fallback = FALSE,
		cv_loss = cv_loss,
		feasible = feasible,
		mult_grid = mult_grid
	)
}
