# Behavior-preserving guardrail for the etwfe_core() / twfeCovs_core()
# unification (#327). The existing `_snaps` lock only PRINTED output; this pins
# the full NUMERIC output of etwfe() / twfeCovs() (which exercise the two OLS
# cores end-to-end) across a fixture matrix, so the refactor that extracts the
# shared `.ols_estimator_core()` reproduces the pre-refactor output to numeric
# tolerance, and any future drift between the two estimators is caught. (The
# refactor is bit-identical by construction -- a pure reorganization with no
# arithmetic change -- but we assert at testthat's default tolerance, which is
# tight enough to catch any substantive regression while staying robust across
# BLAS/platforms.)
#
# The fingerprint flattens every numeric output the cores produce (overall ATT --
# including the asymptotically-exact two-sample ATT when indep_counts is supplied
# (the c7_indep configs), which exercises the indep getTeResultsOLS() call -- the
# full catt_df, and the coefficient vector) and reduces it to collision-resistant
# aggregates; the position-weighted sum (`wsum`) catches reorderings/swaps that
# plain sum/sumsq would miss. Goldens captured from the pre-refactor code (commit
# f1c8cca).

.fp_327 <- function(fit) {
	v <- c(
		fit$att_hat,
		fit$att_se,
		fit$att_p_value,
		as.numeric(fit$catt_hats),
		as.numeric(fit$catt_ses),
		as.numeric(unlist(
			fit$catt_df[c("estimate", "se", "ci_low", "ci_high", "p_value")]
		)),
		as.numeric(fit$beta_hat)
	)
	v <- v[is.finite(v)]
	c(
		n = length(v),
		sum = sum(v),
		sumsq = sum(v^2),
		sumabs = sum(abs(v)),
		wsum = sum(v * seq_along(v)),
		mx = max(v),
		mn = min(v)
	)
}

.mk_327 <- function(
	estimator,
	G,
	T,
	d,
	se_type = "default",
	add_ridge = FALSE,
	reml = FALSE,
	indep = FALSE,
	seed = 1
) {
	co <- genCoefs(
		G = G,
		T = T,
		d = d,
		density = 0.5,
		eff_size = 2,
		seed = seed
	)
	dat <- simulateData(
		co,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = seed
	)
	args <- list(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		se_type = se_type,
		add_ridge = add_ridge,
		verbose = FALSE
	)
	if (!reml) {
		args$sig_eps_sq <- 1
		args$sig_eps_c_sq <- 0.5
	}
	# Exercise the asymptotically-exact two-sample ATT path: indep_counts routes
	# the overall ATT through the indep getTeResultsOLS() call, a branch the other
	# six configs never reach. The refactor reorganized this branch, so pin its
	# output to lock byte-identity there too.
	if (indep) {
		args$indep_counts <- dat$indep_counts
	}
	do.call(estimator, args)
}

.CONFIGS_327 <- list(
	c1_base = list(G = 3, T = 5, d = 2),
	c2_cluster = list(G = 3, T = 5, d = 2, se_type = "cluster"),
	c3_ridge = list(G = 3, T = 5, d = 2, add_ridge = TRUE),
	c4_d0 = list(G = 3, T = 5, d = 0),
	c5_bigGT = list(G = 4, T = 6, d = 3),
	c6_reml = list(G = 3, T = 5, d = 2, reml = TRUE),
	# Same DGP as c1_base but with indep_counts supplied, so the overall ATT comes
	# from the indep (two-sample) getTeResultsOLS() call -- a branch otherwise
	# unexercised by this guardrail.
	c7_indep = list(G = 3, T = 5, d = 2, indep = TRUE)
)

# Goldens: fingerprint of each (estimator, config) from the pre-#327 code.
.GOLDEN_327 <- list(
	etwfe_c1_base = c(
		n = 74,
		sum = 60.8352899168,
		sumsq = 304.06524078,
		sumabs = 108.225547072,
		wsum = 2231.05132864,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c2_cluster = c(
		n = 74,
		sum = 60.2500429937,
		sumsq = 302.901520801,
		sumabs = 107.200328027,
		wsum = 2222.88197939,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c3_ridge = c(
		n = 74,
		sum = 60.8355894105,
		sumsq = 304.068391429,
		sumabs = 108.226076316,
		wsum = 2231.06287584,
		mx = 6.14711965305,
		mn = -4.35014739268
	),
	etwfe_c4_d0 = c(
		n = 40,
		sum = -38.7019017011,
		sumsq = 113.245976802,
		sumabs = 47.2982242601,
		wsum = -917.585716081,
		mx = 0.837538966036,
		mn = -4.3645520168
	),
	etwfe_c5_bigGT = c(
		n = 126,
		sum = 103.329393934,
		sumsq = 503.04966703,
		sumabs = 188.023020045,
		wsum = 7023.93349848,
		mx = 5.05575657097,
		mn = -6.08136842374
	),
	etwfe_c6_reml = c(
		n = 74,
		sum = 61.0043302714,
		sumsq = 304.480412608,
		sumabs = 108.533286057,
		wsum = 2233.64525524,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c7_indep = c(
		n = 74,
		sum = 60.7968737847,
		sumsq = 304.0921750954,
		sumabs = 108.1871309400,
		wsum = 2230.6565019494,
		mx = 6.1470871639,
		mn = -4.3501244277
	),
	twfeCovs_c1_base = c(
		n = 36,
		sum = 40.3928807217,
		sumsq = 143.179199422,
		sumabs = 40.684445326,
		wsum = 1055.74546402,
		mx = 8.48022246831,
		mn = -0.145782302154
	),
	twfeCovs_c2_cluster = c(
		n = 36,
		sum = 45.6229901083,
		sumsq = 161.68019933,
		sumabs = 51.4796001752,
		wsum = 1137.28023272,
		mx = 8.48022246831,
		mn = -1.64492081469
	),
	twfeCovs_c3_ridge = c(
		n = 36,
		sum = 40.3929796973,
		sumsq = 143.179930206,
		sumabs = 40.684542513,
		wsum = 1055.74808929,
		mx = 8.48024443005,
		mn = -0.145781407869
	),
	twfeCovs_c4_d0 = c(
		n = 34,
		sum = -35.0606655433,
		sumsq = 81.3634195152,
		sumabs = 42.0840786104,
		wsum = -620.146197574,
		mx = 0.874312552911,
		mn = -2.86568612327
	),
	twfeCovs_c5_bigGT = c(
		n = 47,
		sum = 30.8793423015,
		sumsq = 190.00469572,
		sumabs = 67.5089511377,
		wsum = 596.753874573,
		mx = 5.11653367631,
		mn = -6.09899877872
	),
	twfeCovs_c6_reml = c(
		n = 36,
		sum = 40.4838258143,
		sumsq = 143.323584631,
		sumabs = 40.8177266326,
		wsum = 1057.30820486,
		mx = 8.48022246831,
		mn = -0.166950409144
	),
	twfeCovs_c7_indep = c(
		n = 36,
		sum = 40.3357038891,
		sumsq = 143.1233086215,
		sumabs = 40.6272684934,
		wsum = 1055.6800476111,
		mx = 8.4802224683,
		mn = -0.1457823022
	)
)

test_that("etwfe() / twfeCovs() numeric output matches the pre-#327 golden to numeric tolerance (#327)", {
	skip_on_cran()
	for (est in c("etwfe", "twfeCovs")) {
		for (nm in names(.CONFIGS_327)) {
			cfg <- .CONFIGS_327[[nm]]
			fit <- do.call(
				.mk_327,
				c(list(estimator = get(est)), cfg)
			)
			key <- paste0(est, "_", nm)
			# Guard that an indep config actually consumed indep_counts -- else
			# the overall ATT would silently fall back to the in-sample branch and
			# the indep coverage this config exists for would be illusory.
			if (isTRUE(cfg$indep)) {
				expect_true(isTRUE(fit$indep_counts_used), info = key)
			}
			expect_equal(.fp_327(fit), .GOLDEN_327[[key]], info = key)
		}
	}
})
