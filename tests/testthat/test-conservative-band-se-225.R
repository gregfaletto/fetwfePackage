library(testthat)
library(fetwfe)

# Issue #225: pin the conservative simultaneous-CI standard-error formula and
# broaden coverage of the recently-added inference surface (#192/#197/#200).
#
# The shipped formula (`.cauchy_schwarz_se()` in R/simultaneous_cis.R, used by
# the `se_type = "conservative"` band at the Sigma_1/Sigma_2 combine step) is
# correct; these tests are drift sentinels.
#
# Why a new test was needed: the pre-existing Test 6 in test-maxt-pvalues-200.R
# recovers the band SE from the band's OWN pointwise CI and then re-derives the
# adjusted p-values from that same SE -- a round-trip that any (even wrong) SE
# satisfies. It also runs under family = "cohort", where each effect is a single
# cohort so the second variance component Sigma_2 = 0 and the conservative SE is
# byte-identical to the tight SE. So it pins the Bonferroni-dual p-value but
# CANNOT catch a regression in the conservative SE formula. These do.

# --- Shared deterministic fixture: one simulated panel; fits built once. ---
# G = 4, T = 6 => the event-study family has cohort-pooling event times, so the
# second variance component Sigma_2 (hence v2) is strictly positive there.
.cs_sim_225 <- local({
	sc <- genCoefs(G = 4, T = 6, d = 2, density = 0.5, eff_size = 1.5, seed = 7)
	simulateData(sc, N = 120, sig_eps_sq = 2, sig_eps_c_sq = 1, seed = 7)
})

.cs_fit_225 <- function(se_type) {
	fetwfe(
		pdata = .cs_sim_225$pdata,
		time_var = .cs_sim_225$time_var,
		unit_var = .cs_sim_225$unit_var,
		treatment = .cs_sim_225$treatment,
		covs = .cs_sim_225$covs,
		response = .cs_sim_225$response,
		sig_eps_sq = .cs_sim_225$sig_eps_sq,
		sig_eps_c_sq = .cs_sim_225$sig_eps_c_sq,
		q = 0.5,
		se_type = se_type,
		verbose = FALSE
	)
}

.cs_fits_225 <- suppressMessages(list(
	default = .cs_fit_225("default"),
	conservative = .cs_fit_225("conservative"),
	cluster = .cs_fit_225("cluster")
))

# Recover the per-effect SE from a simultaneous_cis object's pointwise band:
# pointwise_ci = estimate +/- pointwise_critical_value * se.
.cs_se_from_band <- function(sci) {
	(sci$ci$pointwise_ci_high - sci$ci$estimate) / sci$pointwise_critical_value
}

# ---------------------------------------------------------------------------
# 1. Unit test of the conservative-SE formula in isolation, against
#    hand-computed values. This is the drift sentinel: any change to the
#    formula (dropping the cross term, a wrong coefficient, a scale factor)
#    makes it fail.
# ---------------------------------------------------------------------------
test_that(".cauchy_schwarz_se is sqrt(v1)+sqrt(v2) and matches hand values (#225)", {
	v1 <- c(4, 0, 1, 9)
	v2 <- c(9, 5, 0, 16)
	got <- fetwfe:::.cauchy_schwarz_se(v1, v2)

	# The Cauchy-Schwarz upper bound sqrt(v1 + v2 + 2*sqrt(v1*v2)) is exactly
	# the sum of the two component SEs sqrt(v1) + sqrt(v2).
	expect_equal(got, sqrt(v1) + sqrt(v2))
	# Hand-computed: (2+3, 0+sqrt(5), 1+0, 3+4).
	expect_equal(got, c(5, sqrt(5), 1, 7))

	# Degenerate components: one term zero collapses to the other's SE.
	expect_equal(fetwfe:::.cauchy_schwarz_se(c(4, 0), c(0, 9)), c(2, 3))
	expect_equal(fetwfe:::.cauchy_schwarz_se(0, 0), 0)

	# It strictly exceeds the TIGHT combined SE sqrt(v1 + v2) exactly where both
	# components are positive (the cross term 2*sqrt(v1*v2) > 0), and equals it
	# where a component is zero. Dropping the cross term is the regression this
	# guards against.
	tight <- sqrt(v1 + v2)
	both <- v1 > 0 & v2 > 0
	expect_true(all(got[both] > tight[both]))
	expect_equal(got[!both], tight[!both])
})

# ---------------------------------------------------------------------------
# 2. Integration: the conservative event-study band is actually WIRED to the
#    formula on a v2 > 0 family. se_type changes only the variance combine (not
#    selection or point estimates), so the conservative and tight bands are
#    row-aligned over identical v1, v2: tight = sqrt(v1 + v2), conservative =
#    sqrt(v1 + v2 + 2*sqrt(v1*v2)). The conservative band must therefore be
#    >= the tight band everywhere and STRICTLY wider for >= 1 cohort-pooling
#    (v2 > 0) event time.
# ---------------------------------------------------------------------------
test_that("conservative event-study band exceeds the tight band where v2>0 (#225)", {
	skip_if_not_installed("mvtnorm")
	sci_t <- simultaneousCIs(.cs_fits_225$default, family = "event_study")
	sci_c <- suppressMessages(
		simultaneousCIs(.cs_fits_225$conservative, family = "event_study")
	)

	# Same selection + point estimates under both se_types -> row-aligned.
	expect_equal(sci_c$ci$estimate, sci_t$ci$estimate)

	se_t <- .cs_se_from_band(sci_t)
	se_c <- .cs_se_from_band(sci_c)
	fin <- is.finite(se_t) & is.finite(se_c) & se_t > 0
	expect_gte(sum(fin), 2L)

	expect_true(all(se_c[fin] >= se_t[fin] - 1e-9))
	expect_true(any(se_c[fin] > se_t[fin] + 1e-6))
})

# ---------------------------------------------------------------------------
# 3. family = "custom" is well-formed under cluster AND conservative se_type
#    (previously only the default/tight custom path was covered). Row 1 pools
#    every treatment-effect cell; row 2 picks a single cell.
# ---------------------------------------------------------------------------
test_that("simultaneousCIs(family='custom') is well-formed under cluster + conservative (#225)", {
	skip_if_not_installed("mvtnorm")
	for (st in c("cluster", "conservative")) {
		fit <- .cs_fits_225[[st]]
		nt <- length(fit$treat_inds)
		C <- rbind(rep(1 / nt, nt), c(1, rep(0, nt - 1)))

		sci <- suppressMessages(
			simultaneousCIs(fit, family = "custom", contrasts = C)
		)
		expect_identical(sci$K, 2L)
		# API contract: estimates are the contrast applied to the (g, t) effects.
		expect_equal(
			sci$ci$estimate,
			as.numeric(C %*% fit$beta_hat[fit$treat_inds]),
			info = st
		)

		se <- .cs_se_from_band(sci)
		expect_true(all(is.finite(se) & se >= 0), info = st)

		# Simultaneous band is at least as wide as the pointwise band (K = 2, so
		# the multiplicity-adjusted critical value strictly exceeds z_{1-a/2}).
		pos <- se > 0
		w_sim <- (sci$ci$simultaneous_ci_high - sci$ci$simultaneous_ci_low)[pos]
		w_pw <- (sci$ci$pointwise_ci_high - sci$ci$pointwise_ci_low)[pos]
		expect_true(all(w_sim >= w_pw - 1e-9), info = st)
	}
})
