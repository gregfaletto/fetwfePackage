# Regression tests for #165 / hardened in #178: the GLS whitening step in
# `.estimate_variance_and_gls()` (R/gls_machinery.R) uses the block-apply
# identity (I_N kron A) %*% vec_{T,N}(M) = vec_{T,N}(A %*% M) instead of
# materialising the (N*T) x (N*T) Kronecker matrix.
#
# Block 1 (sanity prologue): the algebraic identity `(I_N kron A) %*% vec(M)
# = vec(A %*% M)` holds for a small numeric fixture. Documents the math the
# optimisation relies on; this assertion is necessary but not sufficient.
#
# Block 2 (load-bearing): directly call `fetwfe:::.estimate_variance_and_gls()`
# and assert its `y_gls` / `X_gls` match an explicit-Kronecker reference
# computed inline, to round-off (`tolerance = 1e-14`, not 0 -- the two sides
# are different floating-point operation orders, so bit-identity is a property
# of the BLAS rather than of correctness; see #427). A future refactor that
# silently breaks the production reshape — e.g. uses the wrong variance
# component in `A`'s scale, drops the `matrix(..., nrow = T)` reshape, or
# any other change that makes the transform not equal to (I_N kron A) %*% y
# — fails this assertion. The pre-#178 test re-derived both sides inline
# and so passed regardless of what the production code did.
#
# Block 3 (sig_eps_c_sq = 0 edge case): same production-call shape as
# Block 2, exercising the degenerate-covariance branch where the closed
# form collapses to `(1/sqrt(sig_eps_sq)) * I_T`. This block keeps
# `tolerance = 0` where Block 2 no longer does, and that is deliberate: on
# this branch `A` is exactly diagonal, so both operation orders are bit-equal
# on any BLAS. See the comment at the assertions for why.

library(testthat)

test_that("block-apply identity `(I_N kron A) %*% vec(M) = vec(A %*% M)` holds (sanity prologue)", {
	# The mathematical identity that justifies the block-apply form. Holds
	# for any conformable `A` (T x T), `y` (N*T-vector in (unit, time)
	# order), `X` (N*T x p in the same order). This is a sanity check on
	# the math, NOT a check on the production code (the production code's
	# correctness is locked by Block 2 below).
	set.seed(1L)
	N <- 7L
	T <- 4L
	p <- 3L
	A <- matrix(rnorm(T * T), nrow = T, ncol = T)
	y <- rnorm(N * T)
	X <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)

	# Explicit Kronecker form.
	y_kron <- as.vector(kronecker(diag(N), A) %*% y)
	X_kron <- kronecker(diag(N), A) %*% X

	# Block-apply form.
	y_block <- as.vector(A %*% matrix(y, nrow = T))
	X_block <- matrix(A %*% matrix(X, nrow = T), nrow = N * T, ncol = p)

	# NOT `tolerance = 0`. The identity is exact in exact arithmetic, but the two
	# sides are different sequences of floating-point operations -- an (N*T)x(N*T)
	# matrix-vector product versus a TxT by Tx(N or p) product -- so they agree
	# only up to round-off, and which way the last bits fall is a property of the
	# BLAS. `tolerance = 0` here passed on macOS/Accelerate and Windows and failed
	# on all four Linux jobs (#427), where the observed disagreement was 1-4 ULPs:
	# max |diff| 4.44e-16 on values of order 1, mean 1.8e-16 to 1.96e-16. The
	# tolerance below is ~22x that maximum and two orders under the 100x cap.
	#
	# The band was measured, not assumed: on this fixture `tolerance = 1e-14`
	# passes a uniform relative perturbation of 1e-14 and fails at 5e-14, and a
	# *single* element wrong by more than ~1e-14 fails -- `all.equal.numeric`
	# defaults to `countEQ = FALSE`, so the mean is taken over only the differing
	# entries and a localized error is not averaged away across the vector. Noise
	# 4.4e-16, threshold 1e-14, any real defect at least twelve orders above.
	expect_equal(y_block, y_kron, tolerance = 1e-14)
	expect_equal(
		X_block,
		matrix(X_kron, nrow = N * T, ncol = p),
		tolerance = 1e-14
	)
})

test_that(".estimate_variance_and_gls() matches explicit Kronecker GLS form (load-bearing)", {
	# Tiny fixture. Variance components supplied so the REML path in
	# `estOmegaSqrtInv()` doesn't fire; we're testing the GLS transform,
	# not the variance-component estimator.
	set.seed(42L)
	N <- 12L
	T <- 5L
	p <- 7L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	# X_ints is the input to the REML estimator and only consulted when
	# variance components are NA; any conformable matrix works here.
	X_ints <- X_mod
	sig_eps_sq <- 2.3
	sig_eps_c_sq <- 1.1

	res <- fetwfe:::.estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		T = T,
		p = p,
		verbose = FALSE
	)

	# Explicit-Kronecker reference. Mirrors `.estimate_variance_and_gls()`'s
	# closed form for `Omega_sqrt_inv` (in `R/gls_machinery.R`) and its scaling
	# `A <- sqrt(sig_eps_sq) * Omega_sqrt_inv`.
	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) *
		(diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv
	big_kron <- kronecker(diag(N), A)
	y_gls_ref <- big_kron %*% y # (N*T) x 1 matrix
	X_gls_ref <- big_kron %*% X_mod # (N*T) x p matrix

	# Production return shapes per `.estimate_variance_and_gls()`: `y_gls` is
	# (N*T) x 1 matrix (NOT a length-N*T vector); `X_gls` is (N*T) x p
	# matrix. `ignore_attr = TRUE` because `kronecker(...)` returns
	# matrices without dimnames whereas the production return may carry
	# names.
	# `tolerance = 1e-14`, not 0, for the reason given in the prologue above: the
	# production path applies `A` block-wise while the reference builds the full
	# Kronecker product, so the two agree only to round-off. Observed on Linux
	# (#427): max |diff| 4.44e-16, mean 1.05e-16, on values of order 1.
	expect_equal(res$y_gls, y_gls_ref, tolerance = 1e-14, ignore_attr = TRUE)
	expect_equal(res$X_gls, X_gls_ref, tolerance = 1e-14, ignore_attr = TRUE)
	expect_identical(res$sig_eps_sq, sig_eps_sq)
	expect_identical(res$sig_eps_c_sq, sig_eps_c_sq)
})

test_that(".estimate_variance_and_gls() handles sig_eps_c_sq = 0 (edge case)", {
	# Degenerate-covariance branch: when `sig_eps_c_sq = 0`, the closed
	# form for `Omega_sqrt_inv` reduces to `(1/sqrt(sig_eps_sq)) * I_T`
	# (per the comment in `.estimate_variance_and_gls()`). Verify the production
	# function still matches the explicit Kronecker form on this branch.
	set.seed(7L)
	N <- 30L
	T <- 4L
	p <- 5L
	y <- rnorm(N * T)
	X_mod <- matrix(rnorm(N * T * p), nrow = N * T, ncol = p)
	X_ints <- X_mod
	sig_eps_sq <- 1.5
	sig_eps_c_sq <- 0

	res <- fetwfe:::.estimate_variance_and_gls(
		y = y,
		X_ints = X_ints,
		X_mod = X_mod,
		sig_eps_sq = sig_eps_sq,
		sig_eps_c_sq = sig_eps_c_sq,
		N = N,
		T = T,
		p = p,
		verbose = FALSE
	)

	J_over_T <- matrix(1 / T, nrow = T, ncol = T)
	Omega_sqrt_inv <- (1 / sqrt(sig_eps_sq)) *
		(diag(T) - J_over_T) +
		(1 / sqrt(sig_eps_sq + T * sig_eps_c_sq)) * J_over_T
	A <- sqrt(sig_eps_sq) * Omega_sqrt_inv
	big_kron <- kronecker(diag(N), A)

	# `tolerance = 0` here is deliberate, not an oversight left behind when
	# Blocks 1-2 were widened (#427). On this branch `sig_eps_c_sq = 0` makes
	# `sig_eps_sq + T * sig_eps_c_sq` bit-identical to `sig_eps_sq`, so the two
	# `1 / sqrt(...)` coefficients above are the same double and `A`'s
	# off-diagonals cancel to *exactly* zero -- verified at T = 4, 5, 7. `A` is
	# then exactly diagonal, every cross term in both operation orders is exactly
	# zero, and `a + 0 = a` under any summation order. Byte-equality here is
	# guaranteed by the branch condition rather than by the BLAS, so it is a real
	# invariant worth pinning exactly.
	expect_equal(res$y_gls, big_kron %*% y, tolerance = 0, ignore_attr = TRUE)
	expect_equal(
		res$X_gls,
		big_kron %*% X_mod,
		tolerance = 0,
		ignore_attr = TRUE
	)
})
