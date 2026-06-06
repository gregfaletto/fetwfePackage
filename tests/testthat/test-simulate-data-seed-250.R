library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Issue #250: simulateData()'s `seed` argument and the flipped default.
#
# As of fetwfe 1.24.0 simulateData() no longer reuses coefs_obj$seed. The new
# contract:
#   - seed = NULL (default): draw from the ambient RNG AND emit a warning whose
#     text contains "ambient".
#   - seed = NA: draw from the ambient RNG silently (caller-controlled RNG).
#   - seed = <numeric scalar>: set.seed(seed) -> reproducible panel.
#   - anything else: stop("seed must be ...").
#
# These default-path calls (the bare simulateData(coefs, ...) without an
# explicit seed) are the ONLY live simulateData() calls in the suite permitted
# to omit a seed; every other live call is pinned for byte-exactness. The
# mutation check on R/gen_data.R (the default flip + the warning) is owned by
# the parent; this file just locks the observable contract green.
# ------------------------------------------------------------------------------

# Small, fast fixture reused across the tests below.
.seed250_coefs <- function() {
	genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 123)
}
.SEED250_N <- 60L

# ------------------------------------------------------------------------------
# 1) Default (seed = NULL) flips to the ambient RNG and warns.
# ------------------------------------------------------------------------------
test_that("simulateData default draws from the ambient RNG and warns (#250)", {
	coefs <- .seed250_coefs()

	# The default path emits a warning mentioning the ambient generator.
	expect_warning(
		simulateData(
			coefs,
			N = .SEED250_N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		),
		"ambient"
	)

	# Two default draws use the ambient generator, so consecutive panels differ
	# (the removed reuse behavior would have made them byte-identical). This is
	# the mutation-bearing assertion: reverting the default flip makes it FAIL.
	a <- suppressWarnings(simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	))
	b <- suppressWarnings(simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	))
	expect_false(identical(a$y, b$y))
})

# ------------------------------------------------------------------------------
# 2) seed = NA is silent and honors the caller's RNG state.
# ------------------------------------------------------------------------------
test_that("simulateData seed = NA is silent and respects the ambient RNG (#250)", {
	coefs <- .seed250_coefs()

	# No warning when the ambient generator is opted into explicitly.
	expect_warning(
		simulateData(
			coefs,
			N = .SEED250_N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = NA
		),
		NA
	)

	# Same outer set.seed() -> identical panel (caller controls the RNG).
	set.seed(20250250L)
	a <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
	set.seed(20250250L)
	b <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
	expect_identical(a$y, b$y)

	# A different outer seed -> a different panel.
	set.seed(99999L)
	c <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
	expect_false(identical(a$y, c$y))
})

# ------------------------------------------------------------------------------
# 3) A numeric seed is silent and reproducible.
# ------------------------------------------------------------------------------
test_that("simulateData seed = <int> is silent and reproducible (#250)", {
	coefs <- .seed250_coefs()

	# No warning when an explicit integer seed is supplied.
	expect_warning(
		simulateData(
			coefs,
			N = .SEED250_N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 456
		),
		NA
	)

	# The same integer seed always yields the same panel, independent of the
	# ambient RNG state.
	set.seed(1L)
	a <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 456
	)
	set.seed(2L)
	b <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 456
	)
	expect_identical(a$y, b$y)
	expect_identical(a$X, b$X)
	expect_identical(a$assignments, b$assignments)

	# A different integer seed yields a different panel.
	d <- simulateData(
		coefs,
		N = .SEED250_N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 789
	)
	expect_false(identical(a$y, d$y))
})

# ------------------------------------------------------------------------------
# 4) Validation: non-numeric / non-scalar seeds error.
# ------------------------------------------------------------------------------
test_that("simulateData rejects invalid seed values (#250)", {
	coefs <- .seed250_coefs()

	expect_error(
		simulateData(
			coefs,
			N = .SEED250_N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = "x"
		),
		"seed must be"
	)
	expect_error(
		simulateData(
			coefs,
			N = .SEED250_N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = c(1, 2)
		),
		"seed must be"
	)
})
