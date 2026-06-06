library(testthat)
library(fetwfe)

# ==============================================================================
# Seed / RNG conventions (#254).
#
# Two roles, two conventions:
#   A (estimator nuisance): fixed seed + restore the caller's RNG via
#     .with_preserved_rng() -- getTes()'s truth integration follows this.
#   B (data generation): ambient / advance the caller's stream --
#     genCoefs()/genCoefsCore()/simulateDataCore() follow this, and `seed = NA`
#     is a silent-ambient value (harmonized with simulateData()).
#
# All getTes() coverage uses a NON-MARGINAL (multinomial) DGP so the call
# routes through .expected_cohort_probs() (the formerly-perturbing path); a
# marginal coefs object would be vacuous here.
# ==============================================================================

# Shared non-marginal fixture (multinomial => non-marginal => routes through
# .expected_cohort_probs()). seed = 42 makes the truth deterministic.
coefs_nm <- genCoefs(
	G = 3,
	T = 5,
	d = 2,
	density = 0.5,
	eff_size = 2,
	assignment_type = "multinomial",
	seed = 42
)

# ------------------------------------------------------------------------------
# (3a) BYTE-EXACT: the non-marginal att_true is pinned. tolerance = 1e-6 catches
# RNG-ordering drift (even a 1-draw shift moves att_true by ~2.5e-4, far above
# 1e-6; a larger reseed change moves it more) while staying robust to platform
# floating-point. The exact value is 1.822751376679849.
# ------------------------------------------------------------------------------
test_that("getTes() att_true is byte-exact on the non-marginal fixture", {
	expect_equal(getTes(coefs_nm)$att_true, 1.822751, tolerance = 1e-6)
})

# ------------------------------------------------------------------------------
# (3b) NON-PERTURBATION (seeded caller): getTes() must not advance the caller's
# .Random.seed (Convention A). This is the core bug-fix regression.
# ------------------------------------------------------------------------------
test_that("getTes() does not perturb the caller's RNG (seeded caller)", {
	set.seed(1)
	before <- .Random.seed
	invisible(getTes(coefs_nm))
	expect_identical(.Random.seed, before)
})

# ------------------------------------------------------------------------------
# (3c) NA-seed coefs: the cross-item blocker regression. A coefs object whose
# $seed is NA must let getTes() run (mc_seed becomes NA -> .with_preserved_rng
# -> .apply_seed treats NA as ambient) AND leave .Random.seed unchanged. Because
# genCoefs(seed = NA) is now legal (item 2b), the freshly-built variant must run
# too.
# ------------------------------------------------------------------------------
test_that("getTes() runs on an NA-seed coefs object and preserves the RNG", {
	coefs_na <- coefs_nm
	coefs_na$seed <- NA

	expect_no_error(getTes(coefs_na))

	set.seed(7)
	before <- .Random.seed
	invisible(getTes(coefs_na))
	expect_identical(.Random.seed, before)
})

test_that("getTes() runs on a freshly genCoefs(seed = NA) non-marginal coefs", {
	expect_no_error(
		getTes(genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			assignment_type = "multinomial",
			seed = NA
		))
	)
})

# ------------------------------------------------------------------------------
# (3d) NULL-seed coefs: getTes() runs and leaves .Random.seed unchanged. No
# att_true pin -- the truth is drawn from the ambient stream (non-deterministic).
# ------------------------------------------------------------------------------
test_that("getTes() runs on a NULL-seed coefs object and preserves the RNG", {
	coefs_null <- coefs_nm
	coefs_null$seed <- NULL

	expect_no_error(getTes(coefs_null))

	set.seed(3)
	before <- .Random.seed
	invisible(getTes(coefs_null))
	expect_identical(.Random.seed, before)
})

# ------------------------------------------------------------------------------
# Item 2: `seed = NA` is a silent-ambient value across the data-gen API, and
# honors an outer set.seed() (Convention B: ambient draw advances the caller's
# stream). genCoefs / genCoefsCore / simulateDataCore must run without error and
# be reproducible under the same outer seed.
# ------------------------------------------------------------------------------
test_that("genCoefs(seed = NA) runs silently and honors an outer set.seed()", {
	expect_no_error(
		genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = NA)
	)

	set.seed(99)
	a <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = NA)
	set.seed(99)
	b <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = NA)
	expect_identical(a$beta, b$beta)
	expect_identical(a$theta, b$theta)
})

test_that("genCoefsCore(seed = NA) runs silently and honors an outer set.seed()", {
	expect_no_error(
		genCoefsCore(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = NA
		)
	)

	set.seed(123)
	a <- genCoefsCore(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = NA
	)
	set.seed(123)
	b <- genCoefsCore(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = NA
	)
	expect_identical(a$beta, b$beta)
	expect_identical(a$theta, b$theta)
})

test_that("simulateDataCore(seed = NA) runs silently and honors an outer set.seed()", {
	beta <- genCoefsCore(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)$beta

	expect_no_error(
		simulateDataCore(
			N = 30,
			T = 5,
			G = 3,
			d = 2,
			sig_eps_sq = 1,
			sig_eps_c_sq = 1,
			beta = beta,
			seed = NA
		)
	)

	set.seed(55)
	a <- simulateDataCore(
		N = 30,
		T = 5,
		G = 3,
		d = 2,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1,
		beta = beta,
		seed = NA
	)
	set.seed(55)
	b <- simulateDataCore(
		N = 30,
		T = 5,
		G = 3,
		d = 2,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1,
		beta = beta,
		seed = NA
	)
	expect_identical(a$y, b$y)
	expect_identical(a$X, b$X)
})

# ------------------------------------------------------------------------------
# Item 2: invalid seeds (non-numeric, length > 1) error with the canonical
# message from .apply_seed().
# ------------------------------------------------------------------------------
test_that("the data-gen functions reject invalid seeds with the canonical message", {
	expect_error(
		genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = "x"),
		regexp = "seed must be"
	)
	expect_error(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = c(1, 2)
		),
		regexp = "seed must be"
	)
	expect_error(
		genCoefsCore(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = "x"
		),
		regexp = "seed must be"
	)
	expect_error(
		genCoefsCore(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = c(1, 2)
		),
		regexp = "seed must be"
	)

	beta <- genCoefsCore(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)$beta
	expect_error(
		simulateDataCore(
			N = 30,
			T = 5,
			G = 3,
			d = 2,
			sig_eps_sq = 1,
			sig_eps_c_sq = 1,
			beta = beta,
			seed = "x"
		),
		regexp = "seed must be"
	)
	expect_error(
		simulateDataCore(
			N = 30,
			T = 5,
			G = 3,
			d = 2,
			sig_eps_sq = 1,
			sig_eps_c_sq = 1,
			beta = beta,
			seed = c(1, 2)
		),
		regexp = "seed must be"
	)
})
