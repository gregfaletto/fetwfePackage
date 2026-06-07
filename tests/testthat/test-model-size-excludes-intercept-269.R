library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #269: the reported model-size diagnostics (lambda.max/min/star_model_size) count
# the number of selected *features* and EXCLUDE the always-present intercept (which
# is row 1 of the coefficient vector). This matches glmnet's `df` convention and the
# @param wording ("selects close to 0 / p features"). Under the old convention each
# size counted the intercept too, so it was one larger.
#
# The discriminating assertion is lambda_star_model_size == (# nonzero SLOPE coefs):
# under the old convention it would equal (# nonzero slopes) + 1.
# ------------------------------------------------------------------------------
test_that("reported model sizes count features and exclude the intercept (#269)", {
	skip_if_not_installed("lme4")

	cf <- genCoefs(G = 4, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 269)
	sim <- simulateData(
		cf,
		N = 160,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 269
	)
	base <- list(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		verbose = FALSE
	)

	check_invariant <- function(res) {
		p <- res$p
		# Coefficient vector is length p+1: intercept at [1], the p slope features at
		# [2:(p+1)]. lambda_star_model_size must equal the nonzero slope count only.
		n_slopes <- sum(res$internal$theta_hat[2:(p + 1)] != 0)
		expect_identical(
			as.integer(res$lambda_star_model_size),
			as.integer(n_slopes)
		)
		# Feature-count path bounds: 0 <= max <= star <= min <= p (was <= p + 1 when
		# the intercept was counted).
		expect_gte(res$lambda.max_model_size, 0L)
		expect_lte(res$lambda.max_model_size, res$lambda_star_model_size)
		expect_lte(res$lambda_star_model_size, res$lambda.min_model_size)
		expect_lte(res$lambda.min_model_size, p)
	}

	# Exercises both selection paths' size computation (BIC: getBetaBIC; CV: getBetaCV)
	# plus the shared .lambda_path_diagnostics() for lambda.max/min sizes.
	check_invariant(do.call(fetwfe, c(base, list(lambda_selection = "bic"))))
	check_invariant(do.call(fetwfe, c(base, list(lambda_selection = "cv"))))
})
