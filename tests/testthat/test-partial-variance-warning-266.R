library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #266: supplying exactly ONE of sig_eps_sq / sig_eps_c_sq used to silently
# discard it -- REML re-estimates BOTH noise variances jointly, so the supplied
# value was overwritten with no signal. .estimate_variance_and_gls() now warns
# when exactly one is supplied. The guard is shared by all four estimators; etwfe
# (pure OLS) is the fastest path through it.
# ------------------------------------------------------------------------------
test_that("supplying exactly one noise variance warns; both/neither is silent (#266)", {
	skip_if_not_installed("lme4") # the NA path estimates via REML (estOmegaSqrtInv)

	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	base <- list(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs
	)
	warns_of <- function(extra) {
		w <- character(0)
		withCallingHandlers(
			suppressMessages(do.call(etwfe, c(base, extra))),
			warning = function(c) {
				w <<- c(w, conditionMessage(c))
				invokeRestart("muffleWarning")
			}
		)
		w
	}

	# Partial provision (either direction) -> warns and names the discarded value.
	w_sq <- warns_of(list(sig_eps_sq = 1, sig_eps_c_sq = NA))
	expect_true(any(grepl("Only one of", w_sq, fixed = TRUE)))
	expect_true(any(grepl("sig_eps_sq = 1", w_sq, fixed = TRUE)))

	w_c <- warns_of(list(sig_eps_sq = NA, sig_eps_c_sq = 0.5))
	expect_true(any(grepl("Only one of", w_c, fixed = TRUE)))
	expect_true(any(grepl("sig_eps_c_sq = 0.5", w_c, fixed = TRUE)))

	# Both supplied OR both omitted -> NO partial-provision warning.
	expect_false(any(grepl(
		"Only one of",
		warns_of(list(sig_eps_sq = 1, sig_eps_c_sq = 0.5)),
		fixed = TRUE
	)))
	expect_false(any(grepl(
		"Only one of",
		warns_of(list(sig_eps_sq = NA, sig_eps_c_sq = NA)),
		fixed = TRUE
	)))
})
