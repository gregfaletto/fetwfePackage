# Regression tests for #205. getGramInv() must degrade gracefully -- emit the
# "not invertible" warning and return list(gram_inv = NA, calc_ses = FALSE) --
# rather than abort with an uncaught LAPACK "system is computationally singular"
# error (or silently invert a garbage near-singular Gram), when the
# selected-feature Gram is (near-)singular. Pre-fix it guarded on an absolute
# 1e-16 floor (below .Machine$double.eps) and then called an UNGUARDED solve();
# the fix uses an rcond-aware relative tolerance and wraps solve() in a tryCatch
# routing to the same graceful fallback.

# Build an N*T x p design whose column-centered Gram, (1/(NT)) X'X, has (to high
# precision) the requested eigenvalues, by scaling an orthonormal, already-
# centered basis. This lets us place the smallest Gram eigenvalue exactly in the
# rcond regime that distinguishes the fixed guard from the old 1e-16 floor --
# the off-diagonal coupling to the tiny eigendirection is ~1e-24, so the small
# eigenvalue is well-resolved (not lost in floating-point noise).
.make_X_for_eigvals <- function(N, T, eigvals) {
	n <- N * T
	p <- length(eigvals)
	set.seed(205)
	M <- scale(matrix(stats::rnorm(n * p), n, p), center = TRUE, scale = FALSE)
	Q <- qr.Q(qr(M))[, seq_len(p), drop = FALSE] # orthonormal, centered basis
	Q %*% diag(sqrt(n * eigvals), p, p)
}

test_that("getGramInv() does not abort on a near-singular Gram (#205)", {
	N <- 60L
	T <- 2L
	# Smallest Gram eigenvalue ~1.5e-16: ABOVE the old absolute 1e-16 floor (so
	# the old guard missed it), but rcond ~1.5e-16 is BELOW solve()'s tolerance
	# (.Machine$double.eps ~2.2e-16), so the old unguarded solve() aborted with an
	# uncaught LAPACK error. The fix's rcond-aware guard catches it first.
	X_final <- .make_X_for_eigvals(N, T, c(1, 1, 1.5e-16))

	# On the OLD code this call either threw (solve abort) or returned
	# calc_ses = TRUE (silent garbage inverse) -- both fail the assertions below.
	res <- expect_warning(
		getGramInv(
			N = N,
			T = T,
			X_final = X_final,
			treat_inds = 1L,
			num_treats = 1L,
			calc_ses = TRUE
		),
		"not invertible"
	)
	expect_false(res$calc_ses)
	expect_true(length(res$gram_inv) == 1L && is.na(res$gram_inv))
})

test_that("getGramInv() degrades gracefully on an exactly-singular Gram (#205)", {
	N <- 20L
	T <- 2L
	set.seed(205)
	z1 <- stats::rnorm(N * T)
	z2 <- stats::rnorm(N * T)
	X_final <- cbind(z1, z1, z2) # duplicated column -> exactly rank-deficient

	res <- expect_warning(
		getGramInv(
			N = N,
			T = T,
			X_final = X_final,
			treat_inds = 1L,
			num_treats = 1L,
			calc_ses = TRUE
		),
		"not invertible"
	)
	expect_false(res$calc_ses)
	expect_true(length(res$gram_inv) == 1L && is.na(res$gram_inv))
})
