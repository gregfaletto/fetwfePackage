# High-dimensional (p >= NT) debiasing direction for debiasedATT() (#31).
#
# The fixed-p debiased path builds the debiasing direction as the exact/ridged
# inverse v = (Sigma_hat + tiny*I)^{-1} a, which is invalid once p >= NT
# (Sigma_hat = X'X/n is singular and the tiny ridge yields an unstable v). In the
# high-dimensional regime v is instead the nodewise (desparsified) Riesz
# representer of paper Theorem `debiased.highdim.thm`, computed as the
# l1-penalized solve below. Ported verbatim from the validated reference
# `simulations/highdim_functions.R::riesz_lasso` in the FETWFE paper repo.

#' Nodewise (desparsified) debiasing direction for the high-dimensional regime
#'
#' @description
#' Solves `v_hat = argmin_v (1/2) v' Sig v - a' v + lambda * ||v||_1` by
#' coordinate descent. By the KKT conditions the minimizer satisfies
#' `|| Sig %*% v_hat - a ||_inf <= lambda` --- the relaxed-inverse feasibility the
#' Theorem `debiased.highdim.thm` remainder bound uses --- and as `lambda -> 0`
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
				next # all-zero column: leave v_j = 0
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
	attr(v, "feasibility") <- max(abs(Sv - a)) # ||Sig v - a||_inf
	attr(v, "iters") <- it
	attr(v, "converged") <- converged # FALSE => lambda likely too small
	v
}

#' Theory-scaled nodewise penalty `lambda_node = c * scale * sqrt(log p / N)`
#'
#' @description The Theorem `debiased.highdim.thm` rate scale. The effective
#'   sample size is the number of clusters `N` (units), not `n = NT`. `scale`
#'   rescales to the units of the constraint `||Sig v - a||_inf` (the accessor
#'   passes `scale = max(abs(a))`).
#' @keywords internal
#' @noRd
lambda_node_default <- function(p, N, c = 1.0, scale = 1.0) {
	c * scale * sqrt(log(p) / N)
}
