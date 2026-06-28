library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #343 Item 1: the internal core validator check_etwfe_core_inputs() rejects
# sig_eps_sq = 0, matching the public validator (.collect_etwfe_input_violations,
# the #185 fix -- 0 makes 1/sqrt(sig_eps_sq) non-finite in the GLS transform).
# This is defense-in-depth: the public validator gates first on every exported
# path, so 0 cannot reach the core via the public API, but the two validators
# should agree. sig_eps_c_sq = 0 stays allowed (no unit-level random effects).
# ------------------------------------------------------------------------------
test_that("check_etwfe_core_inputs() rejects sig_eps_sq = 0, allows sig_eps_c_sq = 0 (#343)", {
	mk <- function(sig_eps_sq = 1, sig_eps_c_sq = 1) {
		fetwfe:::check_etwfe_core_inputs(
			in_sample_counts = c("0" = 2L, "2" = 3L),
			N = 5L,
			T = 3L,
			sig_eps_sq = sig_eps_sq,
			sig_eps_c_sq = sig_eps_c_sq,
			indep_counts = NA,
			verbose = FALSE,
			alpha = 0.05,
			add_ridge = FALSE
		)
	}
	expect_no_error(mk()) # valid: positive idiosyncratic variance
	expect_error(mk(sig_eps_sq = 0), "sig_eps_sq") # now rejected (was >= 0)
	expect_no_error(mk(sig_eps_c_sq = 0)) # 0 unit-level variance still allowed
})

# ------------------------------------------------------------------------------
# #343 Item 4c: is_twfe_covs x add_ridge x indep_counts interact cleanly. These
# three are independent in the code (the collapse sets p before the ridge rows
# are appended, and indep_counts only feeds the cohort-probability / ATT branch),
# but the combined path was untested. Pin that it runs and yields a finite ATT.
# ------------------------------------------------------------------------------
test_that("twfeCovs() runs with add_ridge + indep_counts together (#343)", {
	cf <- genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 2, seed = 343)
	sim <- simulateData(
		cf,
		N = 90,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 343
	)
	# the simulated object supplies indep_counts; add_ridge = TRUE turns on the
	# ridge-row augmentation on the collapsed design.
	res <- twfeCovsWithSimulatedData(sim, add_ridge = TRUE)
	expect_s3_class(res, "twfeCovs")
	expect_true(is.finite(res$att_hat))
	expect_true(isTRUE(res$indep_counts_used))

	# Prove add_ridge is actually active on this combined path (not silently
	# dropped): the ridge-augmented fit must differ from the unaugmented one.
	res_noridge <- twfeCovsWithSimulatedData(sim, add_ridge = FALSE)
	expect_true(is.finite(res_noridge$att_hat))
	expect_false(isTRUE(all.equal(res$att_hat, res_noridge$att_hat)))
})
