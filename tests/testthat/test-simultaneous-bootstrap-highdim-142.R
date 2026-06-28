# Tests for the high-dimensional (full p >= NT) regime of the simultaneousCIs()
# multiplier bootstrap (#142 Phase 2). There the analytic Gram inverse need not
# exist, so the bootstrap uses the FULL-design desparsified construction
# (debiasedATT()'s #31 nodewise directions generalized to K effects -- uniformly
# valid, NOT post-selection). EXPERIMENTAL. The simulator can make full p > NT
# panels since #293 (small cohorts).

.make_highdim_sim_fit <- function() {
	# G=3, T=5, d=22 => full p = 390 > NT = 200 at N = 40.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
}

hd_fit <- .make_highdim_sim_fit()

test_that("the high-dim fixture really is full p >= NT with valid SEs", {
	expect_gte(ncol(hd_fit$internal$X_final), nrow(hd_fit$internal$X_final))
	expect_true(isTRUE(hd_fit$internal$calc_ses))
})

test_that("the bootstrap runs in the high-dim regime and returns nodewise diagnostics", {
	for (fam in c("cohort", "all_post_treatment")) {
		bo <- simultaneousCIs(
			hd_fit,
			family = fam,
			method = "bootstrap",
			B = 1000,
			seed = 1
		)
		expect_s3_class(bo, "simultaneous_cis")
		expect_identical(bo$regime, "high-dimensional")
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_true(all(
			c("feasibility", "converged", "lambda_node") %in% names(bo)
		))
		# valid studentized sup-t: critical value in [pointwise, Bonferroni]
		expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
		expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	}
})

test_that("at the default lambda_c the nodewise directions converge and are feasible", {
	bo <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_true(all(bo$converged))
	expect_true(all(bo$feasibility <= bo$lambda_node * (1 + 1e-9)))
})

test_that("too-small lambda_c yields infeasible directions and warns (experimental)", {
	expect_warning(
		simultaneousCIs(
			hd_fit,
			family = "cohort",
			method = "bootstrap",
			B = 500,
			seed = 1,
			lambda_c = 0.05
		),
		"feasibility|experimental"
	)
})

test_that("high-dim bootstrap is deterministic given a seed", {
	a <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 800,
		seed = 7
	)
	b <- simultaneousCIs(
		hd_fit,
		family = "cohort",
		method = "bootstrap",
		B = 800,
		seed = 7
	)
	expect_identical(a$critical_value, b$critical_value)
})

test_that("event_study bootstrap uses the desparsified path in the high-dim regime (#299)", {
	# #299: a full p >= NT event_study fit routes through the desparsified
	# `targets` path (the regression-only guard was lifted), so it carries BOTH the
	# full-design desparsified regression channel and the propensity channel
	# `F_pi`, with `regime == "high-dimensional"` + nodewise diagnostics -- not the
	# fixed-p selected-support fallback it used before #299.
	bo <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 500,
		seed = 1
	)
	expect_s3_class(bo, "simultaneous_cis")
	expect_identical(bo$regime, "high-dimensional")
	# nodewise diagnostics attached (the desparsified-path signature)
	expect_true(all(
		c("feasibility", "converged", "lambda_node") %in% names(bo)
	))
	expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
})

test_that("high-dim event_study bootstrap is deterministic given a seed", {
	a <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 800,
		seed = 11
	)
	b <- simultaneousCIs(
		hd_fit,
		family = "event_study",
		method = "bootstrap",
		B = 800,
		seed = 11
	)
	expect_identical(a$critical_value, b$critical_value)
	expect_identical(a$ci$estimate, b$ci$estimate)
	expect_identical(a$ci$simultaneous_ci_low, b$ci$simultaneous_ci_low)
})

test_that("too-small lambda_c warns for high-dim event_study (experimental)", {
	expect_warning(
		simultaneousCIs(
			hd_fit,
			family = "event_study",
			method = "bootstrap",
			B = 500,
			seed = 1,
			lambda_c = 0.05
		),
		"feasibility|experimental"
	)
})

# ------------------------------------------------------------------------------
# Debiasing anchor (#299, mutation-checkable cross-check vs the validated
# debiasedATT()). The high-dim band center is the DEBIASED estimate
# `estimates + colSums(F_reg)/(N*T)` (Theorem 6.6). For a 1 x num_treats `custom`
# contrast that matches debiasedATT()'s overall-ATT direction (cohort g loads its
# treat_inds cells with cohort_probs[g] / (#cells in cohort g)), that center must
# equal debiasedATT(fit)$att. RHS is the independently-validated single-effect
# debiased ATT; LHS is the simultaneousCIs high-dim center -- an independent path,
# not a round-trip. Mutation-checkable against the `colSums(F_mat)/(N*T)` line.
# ------------------------------------------------------------------------------
test_that("high-dim debiased band center equals debiasedATT()$att for the overall-ATT contrast", {
	G <- hd_fit$G
	T_ <- hd_fit$T
	treat_inds <- hd_fit$treat_inds
	num_treats <- length(treat_inds)
	cohort_probs <- hd_fit$cohort_probs
	# cohort of each treatment cell: T-1 cells for cohort 1, ..., T-G for cohort G
	cohort_of_treat <- rep(seq_len(G), times = (T_ - 1):(T_ - G))
	contrast <- numeric(num_treats)
	for (g in seq_len(G)) {
		idx <- which(cohort_of_treat == g)
		contrast[idx] <- cohort_probs[g] / length(idx)
	}
	# overall-ATT direction: weights sum to sum(cohort_probs[1:G])
	expect_equal(sum(contrast), sum(cohort_probs[seq_len(G)]))
	contrast_mat <- matrix(contrast, nrow = 1L)

	sc <- simultaneousCIs(
		hd_fit,
		family = "custom",
		contrasts = contrast_mat,
		method = "bootstrap",
		B = 200,
		seed = 1
	)
	expect_identical(sc$regime, "high-dimensional")
	db <- debiasedATT(hd_fit)
	# the load-bearing cross-check: equal to ~1e-9 (observed |diff| ~4e-16).
	expect_equal(sc$ci$estimate, db$att, tolerance = 1e-9)
})

# ------------------------------------------------------------------------------
# Per-effect q=1 plug-in pin (#303, post-exec review). The cross-check above pins
# ONE direction (via debiasedATT). This pins a MULTI-effect (K = 2) custom family
# against a fully INDEPENDENT reconstruction: the q=1 nuisance plug-in
# `sum(a_theta * theta_q1)` plus the manual nodewise correction
# `mean((X v) resid_q1)` -- built here with `riesz_lasso` + a hand-rolled score,
# NOT the implementation's `.build_regression_if_highdim`. Mutation-checkable:
# reverting the high-dim band center from the q=1 plug-in (`estimates_q1`) to the
# bridge `estimates` shifts these per-effect centers and fails the match.
# ------------------------------------------------------------------------------
test_that("high-dim custom per-effect centers match an independent q=1 reconstruction (#303)", {
	X <- hd_fit$internal$X_final
	y <- as.numeric(hd_fit$internal$y_final)
	n <- nrow(X)
	p <- ncol(X)
	Tt <- hd_fit$T
	G <- hd_fit$G
	ti <- hd_fit$treat_inds
	nt <- length(ti)
	A <- genFullInvFusionTransformMat(
		first_inds = getFirstInds(G = G, T = Tt),
		T = Tt,
		G = G,
		d = hd_fit$d,
		num_treats = nt,
		fusion_structure = hd_fit$fusion_structure,
		d_inv_treat = hd_fit$internal$d_inv_treat
	)
	theta_q1 <- fetwfe:::.fit_q1_nuisance(X, y, n / Tt, Tt)
	resid <- as.numeric(y - theta_q1[1] - X %*% theta_q1[-1])
	Sig <- crossprod(X) / n
	# two distinct multi-effect contrasts (K = 2): first and last treatment cell.
	C <- rbind(c(1, rep(0, nt - 1)), c(rep(0, nt - 1), 1))
	recon <- apply(C, 1, function(cc) {
		ab <- numeric(p)
		ab[ti] <- cc
		ath <- as.numeric(crossprod(A, ab))
		ln <- lambda_node_default(
			p = p,
			N = n / Tt,
			const = 1.0,
			scale = max(abs(ath))
		)
		v <- riesz_lasso(Sig, ath, ln)
		sum(ath * theta_q1[-1]) + mean((X %*% v) * resid)
	})
	sc <- simultaneousCIs(
		hd_fit,
		family = "custom",
		contrasts = C,
		method = "bootstrap",
		B = 100,
		seed = 1
	)
	expect_identical(sc$regime, "high-dimensional")
	# centers are B/seed-independent (plug-in + deterministic correction).
	expect_equal(sc$ci$estimate, recon, tolerance = 1e-9)
})

# ------------------------------------------------------------------------------
# Band / adjusted-p duality under debiasing (#299, plan-review item 1). The
# adjusted-p `t_stat` MUST use the debiased center `estimates_used`, not the bridge
# `estimates`, or the band and its adjusted p-values disagree on a bridge-zeroed
# effect. Fixture: seed = 2 (same G=3,T=5,d=22,N=40 high-dim DGP), where the
# bridge zeroes SOME-but-not-all event-study effects (bridge e3 == 0 exactly, its
# debiased center is non-zero). The mutation-checkable signal: e3's adjusted
# p-value is finite and < 1 here (debiased |center|/se ~1.5); reverting the
# `t_stat` to abs(estimates) would make e3's t_stat abs(0)/se = 0 and its adjusted
# p ~1. NOTE: seed = 6 instead zeroes ALL effects and trips the upstream all-zero
# early-exit (which now warns in the high-dim regime -- see the dedicated test
# below), so it never reaches this duality path.
# ------------------------------------------------------------------------------
test_that("high-dim event_study band/adjusted-p duality holds with the debiased center (#299)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 2
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 2
	)
	dat$indep_counts <- NA
	fit2 <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	expect_gte(ncol(fit2$internal$X_final), nrow(fit2$internal$X_final))
	# the fixture genuinely zeroes SOME (not all) treatment cells in the bridge,
	# so the path reaches the bootstrap AND has a bridge-zero effect.
	n_sel <- sum(fit2$internal$theta_hat[-1][fit2$treat_inds] != 0)
	expect_gt(n_sel, 0L)
	expect_lt(n_sel, length(fit2$treat_inds))

	bo <- simultaneousCIs(
		fit2,
		family = "event_study",
		method = "bootstrap",
		B = 4000,
		seed = 1
	)
	expect_identical(bo$regime, "high-dimensional")
	# bridge ES estimate for the last event time is exactly 0; the debiased center
	# is not -- the configuration that makes the t_stat source decisive.
	an <- simultaneousCIs(fit2, family = "event_study", method = "analytic")
	bridge_zero <- which(abs(an$ci$estimate) < 1e-10)
	expect_gt(length(bridge_zero), 0L)
	expect_true(any(abs(bo$ci$estimate[bridge_zero]) > 1e-6))

	# Exact dual: an effect's simultaneous band excludes 0 iff its adjusted p-value
	# is < alpha (uses the debiased center on BOTH sides).
	excl0 <- (bo$ci$simultaneous_ci_low > 0) | (bo$ci$simultaneous_ci_high < 0)
	sig <- bo$adjusted_p_values < bo$alpha
	expect_identical(excl0, sig)

	# Mutation anchor: a bridge-zero effect's adjusted p-value is driven by its
	# DEBIASED center, so it is < 1 (here ~0.28); under the bug (abs(estimates))
	# it would be ~1 (mean(boot_max >= 0)).
	expect_true(all(bo$adjusted_p_values[bridge_zero] < 0.95))
})

# ------------------------------------------------------------------------------
# High-dim all-zero degenerate band WARNS, not a silent message (#299 review /
# #304). When the bridge zeroes EVERY treatment effect at p >= NT, the upstream
# degenerate early-exit returns an all-zero band BEFORE the debiased bootstrap
# path. In the high-dim regime selection is NOT consistent and the debiased
# center can be non-zero, so that all-zero band must not be read as "all effects
# are zero" -- the early-exit warns() (instead of the fixed-p message()). Fixture:
# seed = 6 (same high-dim DGP) zeroes all 9 treatment cells. Mutation-checkable:
# reverting the early-exit to a plain message() (or dropping its `p >= N * T_`
# guard) drops the warning and fails expect_warning(). Fixed-p behavior (a
# message(), no warning) is unchanged -- covered by test-degenerate-jacobian-225.R.
# ------------------------------------------------------------------------------
test_that("high-dim all-zero-bridge bootstrap returns an INFORMATIVE band, not the degenerate one (#304)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 6
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 6
	)
	dat$indep_counts <- NA
	fit6 <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	# Fixture invariants: genuinely high-dim AND the bridge zeroed EVERY treat cell.
	expect_gte(ncol(fit6$internal$X_final), nrow(fit6$internal$X_final))
	expect_true(all(fit6$internal$theta_hat[-1][fit6$treat_inds] == 0))

	# #304: the high-dim bootstrap now FALLS THROUGH to the desparsified band -- the
	# debiased centers (Theorem 6.6) are non-zero and informative -- instead of the
	# all-zero degenerate early-exit. So NO "zeroed out every" degenerate warning,
	# and the band is informative.
	# Capture warnings: assert NO degenerate ("zeroed out every") warning fires (the
	# #304 fall-through is taken, not the all-zero early-exit), while tolerating the
	# standing high-dim experimental-feasibility warning (orthogonal; fit-specific).
	w <- character(0)
	bo <- withCallingHandlers(
		simultaneousCIs(
			fit6,
			family = "cohort",
			method = "bootstrap",
			B = 300,
			seed = 1
		),
		warning = function(cond) {
			w <<- c(w, conditionMessage(cond))
			invokeRestart("muffleWarning")
		}
	)
	expect_false(any(grepl("zeroed out every|degenerate", w)))
	expect_identical(bo$regime, "high-dimensional")
	expect_false(all(bo$ci$estimate == 0)) # informative debiased centers
	ses <- (bo$ci$simultaneous_ci_high - bo$ci$estimate) / bo$critical_value
	expect_true(all(is.finite(ses)) && all(ses > 0))
	expect_true(any(
		bo$ci$simultaneous_ci_low != bo$ci$simultaneous_ci_high
	)) # non-degenerate

	# event_study (the headline Theorem 6.6 family, which builds a J_list) must ALSO
	# fall through to an informative band, not crash on the empty selected support.
	wes <- character(0)
	boES <- withCallingHandlers(
		simultaneousCIs(
			fit6,
			family = "event_study",
			method = "bootstrap",
			B = 200,
			seed = 1
		),
		warning = function(cond) {
			wes <<- c(wes, conditionMessage(cond))
			invokeRestart("muffleWarning")
		}
	)
	expect_identical(boES$regime, "high-dimensional")
	expect_false(all(boES$ci$estimate == 0))
	expect_false(any(grepl("zeroed out every|degenerate", wes)))

	# #303 identity: the overall-ATT band center == debiasedATT()$att at lambda_c=1.0.
	G <- fit6$G
	ti <- fit6$treat_inds
	Tt <- fit6$T
	sizes <- (Tt - 1):(Tt - G)
	cot <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(length(ti))
	for (g in seq_len(G)) {
		a_beta[which(cot == g)] <- fit6$cohort_probs[g] / sum(cot == g)
	}
	sc <- simultaneousCIs(
		fit6,
		family = "custom",
		contrasts = matrix(a_beta, nrow = 1),
		method = "bootstrap",
		B = 100,
		seed = 1
	)
	expect_equal(sc$ci$estimate, debiasedATT(fit6)$att, tolerance = 1e-9)
})

test_that("gls=FALSE high-dim all-zero-bridge: calc_ses=FALSE, analytic errors, bootstrap informative (#304)", {
	# Same seed=6 all-zero-bridge fixture, but gls=FALSE (no supplied variances).
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 6
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 6
	)
	dat$indep_counts <- NA
	fitF <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
	expect_true(all(fitF$internal$theta_hat[-1][fitF$treat_inds] == 0)) # all zeroed
	# Part B: the degenerate gls=FALSE fit now reports calc_ses=FALSE (was TRUE) with
	# NA SE (so the C1 SE-consistency validator passes at fit construction).
	expect_false(isTRUE(fitF$internal$calc_ses))
	expect_true(is.na(fitF$att_se))
	# Part B reconciliation: analytic now errors at the calc_ses gate (was: returned
	# the degenerate band because calc_ses was wrongly TRUE).
	expect_error(
		simultaneousCIs(fitF, family = "cohort", method = "analytic"),
		"calc_ses = FALSE"
	)
	# Part A: bootstrap still falls through to the informative desparsified band
	# (no degenerate warning; the experimental-feasibility warning is tolerated).
	wF <- character(0)
	boF <- withCallingHandlers(
		simultaneousCIs(
			fitF,
			family = "cohort",
			method = "bootstrap",
			B = 200,
			seed = 1
		),
		warning = function(cond) {
			wF <<- c(wF, conditionMessage(cond))
			invokeRestart("muffleWarning")
		}
	)
	expect_false(any(grepl("zeroed out every|degenerate", wF)))
	expect_identical(boF$regime, "high-dimensional")
	expect_false(all(boF$ci$estimate == 0))
})

# ------------------------------------------------------------------------------
# Non-fetwfe (betwfe) high-dim bootstrap falls back to the fixed-p band, NOT an
# error (#305 review, Major). The desparsified `targets` path is fetwfe-only; a
# non-fetwfe p >= NT fit must leave `targets` NULL and route through the fixed-p
# selected-support construction (regime == "fixed-p") -- for BOTH the propensity-
# channel `event_study` family and a regression-channel family (`cohort`).
# Mutation-checkable: reverting the dispatch guard `if (p >= N*T_ && is_fetwfe)`
# to `if (p >= N*T_)` reinstates the `stop("...only for fetwfe() fits")`, making
# both calls error and failing the expect_s3_class / regime assertions here.
# ------------------------------------------------------------------------------
test_that("non-fetwfe (betwfe) high-dim bootstrap uses the fixed-p band, not an error (#305)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	bfit <- betwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	# fixture invariant: genuinely high-dim (full p >= NT)
	expect_gte(ncol(bfit$internal$X_final), nrow(bfit$internal$X_final))

	for (fam in c("event_study", "cohort")) {
		# #308: the fallback band substantially under-covers (the selected support
		# here is non-empty, so this is NOT the #304 all-zero degenerate path), so a
		# warning fires; the band is still returned (preserving #305).
		expect_warning(
			bo <- simultaneousCIs(
				bfit,
				family = fam,
				method = "bootstrap",
				B = 200,
				seed = 1
			),
			"under-?cover|unreliable"
		)
		expect_s3_class(bo, "simultaneous_cis")
		# load-bearing: non-fetwfe HD routes through the fixed-p selected-support
		# path (NOT the fetwfe-only desparsified path, NOT a stop()).
		expect_identical(bo$regime, "fixed-p")
		expect_false("feasibility" %in% names(bo)) # no nodewise diagnostics
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_true(all(is.finite(bo$ci$simultaneous_ci_high)))
	}
})
