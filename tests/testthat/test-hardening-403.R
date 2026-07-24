library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #403 defense-in-depth hardening batch. Each test reproduces a latent trap that
# is NOT reachable as a live bug through today's shipped entry points, so it
# drives the internal function directly with a constructed edge input. Each is
# red-green: it FAILS on the pre-fix line and PASSES after (mutation-checked
# during development). Item 2 (CV ridge pseudo-rows) is resolved by documenting
# the acceptance -> behavior-neutral, no test.
# ------------------------------------------------------------------------------

# --- Item 1: getBetaBIC() warns + floors on an interpolating lambda -----------
test_that("getBetaBIC() warns and floors when a lambda interpolates (mse ~ 0) (#403)", {
	N <- 10L
	T <- 2L
	p <- 3L
	set.seed(403)
	X_mod <- matrix(stats::rnorm(N * T * p), nrow = N * T, ncol = p)
	beta_true <- c(1.5, -0.5, 2)
	y <- as.numeric(X_mod %*% beta_true) # exact fit -> one lambda interpolates
	# fit$beta is (p + 1) x n_lambda: col 1 interpolates (eta = 0, beta_true);
	# col 2 is the null model (eta = mean(y), zero slopes).
	fit <- list(beta = cbind(c(0, beta_true), c(mean(y), 0, 0, 0)))
	expect_warning(
		getBetaBIC(
			fit,
			N = N,
			T = T,
			p = p,
			X_mod = X_mod,
			y = y,
			scale_center = rep(0, p),
			scale_scale = rep(1, p)
		),
		"interpolates"
	)
})

# --- Item 3: getGramInv() keeps a 1x1 matrix with one selected treat feature ---
test_that("getGramInv() returns a matrix (not a scalar) with one selected treatment feature (#403)", {
	N <- 40L
	T <- 2L
	set.seed(403)
	# Well-conditioned 3-column design so getGramInv succeeds (calc_ses = TRUE).
	X_final <- matrix(stats::rnorm(N * T * 3), nrow = N * T, ncol = 3)
	res <- getGramInv(
		N = N,
		T = T,
		X_final = X_final,
		treat_inds = 1L, # column 1 is the single treatment feature
		num_treats = 1L,
		calc_ses = TRUE
	)
	expect_true(res$calc_ses)
	# Pre-fix: gram_inv[TRUE-once, TRUE-once] drops to a scalar (is.matrix FALSE),
	# so the following nrow()/ncol() guards compare NULL and pass vacuously.
	expect_true(is.matrix(res$gram_inv))
	expect_identical(dim(res$gram_inv), c(1L, 1L))
})

# --- Item 4: .build_regression_if() gives an actionable singular-Gram error ----
test_that(".build_regression_if() routes a singular centered Gram to an actionable error (#403)", {
	N <- 20L
	T <- 2L
	n <- N * T
	set.seed(403)
	z <- stats::rnorm(n)
	X_sel <- cbind(z, z, stats::rnorm(n)) # duplicated column -> singular Sig
	y <- stats::rnorm(n)
	Psi_full <- matrix(stats::rnorm(3 * 2), nrow = 3, ncol = 2)
	# Pre-fix: bare LAPACK "system is computationally singular" (does not match).
	expect_error(
		.build_regression_if(
			X_sel = X_sel,
			y = y,
			N = N,
			T = T,
			Psi_full = Psi_full
		),
		"not invertible"
	)
})

# --- Item 5: .build_propensity_if() count guard actually bites -----------------
test_that(".build_propensity_if() rejects non-integral N * cohort_probs (#403)", {
	G <- 3L
	N <- 100L
	T <- 5L
	K <- 2L
	set.seed(403)
	A <- matrix(stats::rnorm(G * K), nrow = G, ncol = K)
	# Non-integral N * pi_hat: N * 0.234 = 23.4. Pre-fix the tautological
	# `sum(n_g) + n_never != N` never fires; the integrality check does.
	bad <- c(0.234, 0.3, 0.4)
	expect_error(
		.build_propensity_if(
			cohort_probs_overall = bad,
			G = G,
			N = N,
			T = T,
			A = A
		),
		"could not recover integer cohort counts"
	)
	# Control: integral N * pi_hat must NOT trip the guard.
	good <- c(0.2, 0.3, 0.4)
	expect_error(
		.build_propensity_if(
			cohort_probs_overall = good,
			G = G,
			N = N,
			T = T,
			A = A
		),
		NA
	)
})

# --- Item 6: augment() accepts a fit whose covs slot is NULL (d = 0 / legacy) --
test_that("augment() does not require the unused covs slot (#403)", {
	cf <- genCoefs(G = 3, T = 4, d = 0, density = 0.5, eff_size = 2, seed = 403)
	sim <- simulateData(
		cf,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 403
	)
	fit <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		verbose = FALSE
	)
	# A d = 0 fit stores covs = character(0); a hand-built / legacy object can
	# legitimately carry covs = NULL (the C8 validator allows it). Use
	# `fit["covs"] <- list(NULL)` so the slot stays PRESENT with a NULL value
	# (`fit$covs <- NULL` would delete the slot and trip the class validator).
	fit["covs"] <- list(NULL)
	stopifnot(is.null(fit$covs), "covs" %in% names(fit))
	# Pre-fix: augment() errors "missing slots covs ... earlier dev build".
	au <- augment(fit, data = sim$pdata)
	expect_s3_class(au, "data.frame")
	expect_true(all(c(".fitted", ".resid") %in% names(au)))
})

# --- Item 7: .expected_cohort_probs() validates distribution before seeding ----
test_that(".expected_cohort_probs() does not touch the RNG on an invalid distribution (#403)", {
	assignment_coefs <- matrix(0, nrow = 3, ncol = 4) # unused before the error
	set.seed(999)
	before <- .Random.seed
	# Pre-fix: .apply_seed(seed = 1) runs BEFORE the distribution is validated, so
	# the ambient RNG is advanced before the error. Post-fix: match.arg() errors
	# first, leaving .Random.seed untouched.
	expect_error(
		.expected_cohort_probs(
			assignment_coefs = assignment_coefs,
			d = 2L,
			distribution = "not_a_distribution",
			seed = 1L
		)
	)
	expect_identical(.Random.seed, before)
})
