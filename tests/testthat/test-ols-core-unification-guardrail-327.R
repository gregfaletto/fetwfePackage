# Behavior-preserving guardrail for the etwfe_core() / twfeCovs_core()
# unification (#327). The existing `_snaps` lock only PRINTED output; this pins
# the NUMERIC output of etwfe() / twfeCovs() (which exercise the two OLS
# cores end-to-end) across a fixture matrix.
#
# The fingerprint flattens these numeric outputs of the cores (overall ATT --
# including the asymptotically-exact two-sample ATT when indep_counts is supplied
# (the c7_indep configs), which exercises the indep getTeResultsOLS() call -- the
# full catt_df, and the coefficient vector) and reduces it to collision-resistant
# aggregates; the position-weighted sum (`wsum`) catches reorderings/swaps that
# plain sum/sumsq would miss.

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
	# output.
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

# Goldens: fingerprint of each (estimator, config).
.GOLDEN_327 <- list(
	etwfe_c1_base = c(
		n = 74,
		sum = 60.8352734697,
		sumsq = 304.0652407664,
		sumabs = 108.2255306250,
		wsum = 2231.0509503580,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c2_cluster = c(
		n = 74,
		sum = 60.2500425622,
		sumsq = 302.9015208008,
		sumabs = 107.2003275954,
		wsum = 2222.8819694671,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c3_ridge = c(
		n = 74,
		sum = 60.8354731328,
		sumsq = 304.0673411971,
		sumabs = 108.2258834551,
		wsum = 2231.0586485055,
		mx = 6.1471088233,
		mn = -4.3501397377
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
		sum = 103.3293928000,
		sumsq = 503.0496670296,
		sumabs = 188.0230189116,
		wsum = 7023.9334644682,
		mx = 5.05575657097,
		mn = -6.08136842374
	),
	etwfe_c6_reml = c(
		n = 74,
		sum = 61.0043063298,
		sumsq = 304.4804125696,
		sumabs = 108.5332621155,
		wsum = 2233.6447045813,
		mx = 6.14708716389,
		mn = -4.35012442767
	),
	etwfe_c7_indep = c(
		n = 74,
		sum = 60.7968573377,
		sumsq = 304.0921750818,
		sumabs = 108.1871144929,
		wsum = 2230.6561236666,
		mx = 6.1470871639,
		mn = -4.3501244277
	),
	twfeCovs_c1_base = c(
		n = 36,
		sum = 28.2113851932,
		sumsq = 174.5039823141,
		sumabs = 48.4404221105,
		wsum = 889.6701501748,
		mx = 8.5910974771,
		mn = -2.0282300035
	),
	twfeCovs_c2_cluster = c(
		n = 36,
		sum = 34.8652761713,
		sumsq = 221.8957842605,
		sumabs = 64.6894233615,
		wsum = 986.9970631622,
		mx = 8.5910974771,
		mn = -3.1887787766
	),
	twfeCovs_c3_ridge = c(
		n = 36,
		sum = 28.2114271933,
		sumsq = 174.5045643108,
		sumabs = 48.4404921913,
		wsum = 889.6715520479,
		mx = 8.5911123145,
		mn = -2.0282324251
	),
	twfeCovs_c4_d0 = c(
		n = 34,
		sum = -29.9224949548,
		sumsq = 96.4450430598,
		sumabs = 39.8082659448,
		wsum = -562.3335794533,
		mx = 1.0096134340,
		mn = -4.2026152336
	),
	twfeCovs_c5_bigGT = c(
		n = 47,
		sum = 20.9958565766,
		sumsq = 151.4808010364,
		sumabs = 58.7333446913,
		wsum = 478.5060668888,
		mx = 4.8163044763,
		mn = -5.3625646504
	),
	twfeCovs_c6_reml = c(
		n = 36,
		sum = 28.3713380297,
		sumsq = 174.8449365058,
		sumabs = 48.7203315275,
		wsum = 892.4322530116,
		mx = 8.5910974771,
		mn = -2.0548339806
	),
	twfeCovs_c7_indep = c(
		n = 36,
		sum = 28.5235274621,
		sumsq = 174.3981994637,
		sumabs = 48.4130630864,
		wsum = 890.2587409090,
		mx = 8.5910974771,
		mn = -2.0282300035
	)
)

test_that("etwfe() / twfeCovs() numeric output matches its golden fingerprint to numeric tolerance (#327)", {
	skip_on_cran()
	skip_if_not_installed("lme4")
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
