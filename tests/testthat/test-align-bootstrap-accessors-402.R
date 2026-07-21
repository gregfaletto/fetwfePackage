# Tests for #402: align the two bootstrap accessors' contracts.
#   Item 1 -- debiasedATT(method = "bootstrap") now supports indep_counts
#             (two-sample) fits (the scalar gate was lifted; the band already ran
#             the identical construction). The `debiasedATT` errors-on-indep test
#             lived in test-debiased-att-wild-bootstrap-360.R and is now inverted
#             there; the substantive correctness proof lives here.
#   Item 2  -- simultaneousCIs() gains `cv_time_budget` (the #384 wall-clock
#             backstop debiasedATT() already offered) and surfaces `$lambda_cv`.

# ---- Fixtures (built once; the fetwfe fits are the expensive part) ---------

# A two-sample fit: fetwfeWithSimulatedData() unconditionally forwards the
# simulator's indep_counts, so every such fit is two-sample.
.mk402_two_sample_fit <- function() {
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 3)
	sim <- simulateData(
		cf,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 3
	)
	fetwfeWithSimulatedData(sim, q = 0.5, verbose = FALSE)
}

# A single-sample high-dimensional fit: G=3, T=5, d=22 => full p = 390 > NT = 200
# at N = 40, so the band takes the desparsified nodewise construction where
# `lambda_c` / `cv_time_budget` are live (they are inert at p < NT).
.mk402_highdim_fit <- function() {
	cf <- genCoefs(G = 3, T = 5, d = 22, density = 0.5, eff_size = 2, seed = 1)
	dat <- simulateData(
		cf,
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

ts_fit_402 <- .mk402_two_sample_fit()
hd_fit_402 <- .mk402_highdim_fit()

# ---- Item 1: two-sample wild bootstrap ------------------------------------

test_that("debiasedATT(method='bootstrap') runs on two-sample (indep_counts) fits and returns valid two-sample inference (#402)", {
	expect_true(isTRUE(ts_fit_402$indep_counts_used)) # guard: genuinely two-sample

	alpha <- 0.05
	# (i) Runs without error -- SUBSTANTIVE, not a smoke test: the internal V2
	#     influence-function anchor stop()s unless the per-unit two-sample
	#     propensity IF reproduces the two-sample `var_weight` to 1e-6, so success
	#     PROVES the two-sample IF decomposition the lifted gate's message wrongly
	#     claimed was impossible (empirically the anchor matches to ~1e-18).
	res <- debiasedATT(
		ts_fit_402,
		method = "bootstrap",
		seed = 1,
		alpha = alpha
	)
	expect_s3_class(res, "debiased_att")
	expect_identical(res$method, "bootstrap")

	# (ii) The studentized bootstrap crit value honours the #363 Gaussian floor and
	#      drives the interval; `se` is the two-channel two-sample SE.
	expect_gte(res$crit_value, qnorm(1 - alpha / 2))
	expect_equal(
		res$ci_high - res$att,
		res$crit_value * res$se,
		tolerance = 1e-8
	)
	expect_equal(
		res$att - res$ci_low,
		res$crit_value * res$se,
		tolerance = 1e-8
	)
	expect_true(is.finite(res$var_weight) && res$var_weight > 0)
	expect_equal(res$se, sqrt(res$var_reg + res$var_weight), tolerance = 1e-10)

	# (iii) Floored, so the bootstrap interval is never narrower than analytic Wald.
	ana <- debiasedATT(ts_fit_402, method = "analytic", alpha = alpha)
	expect_gte((res$ci_high - res$ci_low) - (ana$ci_high - ana$ci_low), -1e-9)

	# (iv) Defense-in-depth, INDEPENDENT of the production runtime anchor: rebuild
	#      the per-unit two-sample propensity IF from the fit's public slots (bridge
	#      catt / att_hat, marginal cohort_probs_overall) and confirm sum(psi2^2)/n^2
	#      reproduces the two-sample cohort-weight variance the bootstrap uses -- the
	#      exact `.build_propensity_if()` construction the wild bootstrap runs. So even
	#      if the internal anchor's tolerance were loosened, an IF mismatch is caught.
	G <- ts_fit_402$G
	cpo <- ts_fit_402$cohort_probs_overall
	S <- sum(cpo[seq_len(G)])
	a_att_G <- (ts_fit_402$catt_df$estimate - ts_fit_402$att_hat) / S
	psi2 <- as.numeric(fetwfe:::.build_propensity_if(
		cohort_probs_overall = cpo,
		G = G,
		N = ts_fit_402$N,
		T = ts_fit_402$T,
		A = matrix(a_att_G, nrow = G, ncol = 1L)
	))
	n <- ts_fit_402$N * ts_fit_402$T
	expect_equal(sum(psi2^2) / n^2, res$var_weight, tolerance = 1e-6)
})

# ---- Item 2a: cv_time_budget validation + threading -----------------------

test_that("simultaneousCIs() validates cv_time_budget (#402)", {
	for (bad in list(-1, 0, NA_real_, c(1, 2), "x")) {
		expect_error(
			simultaneousCIs(ts_fit_402, cv_time_budget = bad),
			"cv_time_budget"
		)
	}
	# The default Inf and any positive number are accepted (validated even on the
	# analytic path, where the budget is inert -- symmetric with debiasedATT()).
	expect_s3_class(
		simultaneousCIs(ts_fit_402, family = "cohort", cv_time_budget = Inf),
		"simultaneous_cis"
	)
})

test_that("cv_time_budget threads from simultaneousCIs() into the band's CV node, and $lambda_cv is surfaced (#402)", {
	# A vanishing budget leaves the CV fold sweep incomplete -> theory-scale
	# fallback with a warning. That the warning fires at all PROVES the budget
	# threaded end-to-end (generic -> method -> impl -> bootstrap -> .cv_lambda_node);
	# drop any link and no warning fires (mutation-checked).
	ws <- character(0)
	sc <- withCallingHandlers(
		simultaneousCIs(
			hd_fit_402,
			family = "event_study",
			method = "bootstrap",
			lambda_c = "cv",
			cv_time_budget = 1e-6,
			seed = 1
		),
		warning = function(w) {
			ws <<- c(ws, conditionMessage(w))
			invokeRestart("muffleWarning")
		}
	)
	expect_true(any(grepl("cv_time_budget", ws)))
	expect_s3_class(sc, "simultaneous_cis")
	expect_identical(sc$regime, "high-dimensional")

	# Item 2b: with lambda_c = "cv" the band surfaces $lambda_cv (the CV diagnostics),
	# matching debiasedATT(); the 1e-6 budget fired, so it records that fallback.
	expect_true("lambda_cv" %in% names(sc))
	expect_setequal(
		names(sc$lambda_cv),
		c(
			"cv_loss",
			"feasible",
			"mult_grid",
			"fallback",
			"fallback_reason",
			"bailed"
		)
	)
	expect_true(isTRUE(sc$lambda_cv$fallback))
	expect_identical(sc$lambda_cv$fallback_reason, "time")
})

# ---- Item 2b: $lambda_cv absent when no CV runs ---------------------------

test_that("simultaneousCIs() omits $lambda_cv on a fixed-lambda_c band (#402)", {
	sc_fixed <- simultaneousCIs(
		hd_fit_402,
		family = "event_study",
		method = "bootstrap",
		lambda_c = 1.0,
		seed = 1
	)
	expect_false("lambda_cv" %in% names(sc_fixed))
})
