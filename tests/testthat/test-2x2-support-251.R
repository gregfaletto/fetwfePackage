# ==============================================================================
# 2x2 and late-adopting single-treatment-effect support (#251, Phase 2 of #112).
#
# #112 added a `num_treats >= 2` guard in prep_for_etwfe_core() that rejected
# panels with a single treatment effect (one post-treatment period). #251 removes
# that guard so num_treats = 1 panels now FIT, down to the true 2x2 case
# (T = 2, G = 1) and the late-adopting single cohort (G = 1, T >= 3, adopting in
# the final period). Motivation: continuity for simulation sweeps over T.
#
# The crux is the warning ASYMMETRY. T == 2 is the unique no-fusion
# configuration (it forces G = 1, num_treats = 1, and a degenerate time-FE block,
# so every fused block is trivial), so fetwfe() warns there that no fusion
# applies and the bridge penalty reduces to individual shrinkage. It is NOT
# equivalent to DiD/ETWFE -- the bridge still shrinks the single effect (spike:
# fetwfe != etwfe on identical data). A late-adopting num_treats = 1 cohort with
# T >= 3 still fuses the TIME fixed effects, so it does NOT warn. betwfe applies
# the bridge to the raw coefficients (never fuses) and etwfe/twfeCovs are
# unpenalized, so none of them warn at T = 2. The warning is therefore
# fetwfe-only, and fires at T == 2 only.
#
# Coverage:
#   1. 2x2 (T = 2, G = 1) and late-adopting (G = 1, T = 5, adopt final period)
#      fixtures: all four estimators fit, accessors return well-formed objects,
#      SEs finite, att_hat is numerically sane (a true effect is injected).
#   2. Warning asymmetry: fetwfe@T=2 warns "no fusion"; fetwfe@late-adopt
#      (num_treats = 1, T >= 3) does NOT warn; fetwfe@num_treats>=2 does NOT
#      warn; betwfe/etwfe@T=2 do NOT warn.
#   3. T-sweep continuity (G = 1, T = 2..6): fetwfe runs at every T, only T = 2
#      warns.
#   4. genCoefs / simulateData round-trip at T = 2.
#   5. Mutation-check (recorded, not executed in CI): deleting the `T == 2`
#      warning in fetwfe() makes the "warns at T = 2" test FAIL.
#
# All assertions are on the canonical `$G` / `$T` (never the deprecated `$R`).
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared fixtures.
# ------------------------------------------------------------------------------

# A hand-built single-cohort panel: `frac` of the units adopt at calendar time
# `adopt` (offset = adopt), the rest are never-treated. When `tau` is non-zero
# the treated post-periods get an additive effect of size `tau`, so ATT recovery
# is non-tautological. The seed stream is fixed for reproducibility.
.s251_make_panel <- function(
	N = 240,
	T = 6,
	adopt = 2,
	tau = 2,
	frac = 0.5,
	seed = 251
) {
	set.seed(seed)
	units <- sprintf("u%03d", seq_len(N))
	treated_units <- units[seq_len(round(N * frac))]
	rows <- do.call(
		rbind,
		lapply(units, function(u) {
			is_tr <- u %in% treated_units
			unit_effect <- rnorm(1)
			time_shocks <- rnorm(T)
			cov1 <- rnorm(T)
			cov2 <- stats::runif(T)
			treat <- as.integer(is_tr & (seq_len(T) >= adopt))
			y <- unit_effect +
				time_shocks +
				0.5 * cov1 +
				tau * treat +
				rnorm(T, sd = 0.3)
			data.frame(
				time = as.integer(seq_len(T)),
				unit = as.character(u),
				treatment = treat,
				cov1 = cov1,
				cov2 = cov2,
				y = y
			)
		})
	)
	rows[order(rows$unit, rows$time), , drop = FALSE]
}

.s251_fit <- function(pdata, estimator = fetwfe, ...) {
	estimator(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		...
	)
}

# ------------------------------------------------------------------------------
# 1a. 2x2 (T = 2, G = 1): all four estimators fit and recover the injected ATT.
#
# tau = 2 is injected, so att_hat recovery is numerical (tolerance), not just
# finiteness. twfeCovs is the naive benchmark -- biased under staggering, but at
# a clean 2x2 with one cohort it should still be in the right ballpark; we assert
# only finiteness for it to stay robust.
# ------------------------------------------------------------------------------

test_that("all four estimators fit a 2x2 (T = 2, G = 1) panel", {
	skip_on_cran()
	p22 <- .s251_make_panel(N = 240, T = 2, adopt = 2, tau = 2, seed = 2202)

	res_f <- suppressWarnings(.s251_fit(p22, fetwfe))
	expect_s3_class(res_f, "fetwfe")
	expect_equal(res_f$G, 1)
	expect_equal(res_f$T, 2)
	expect_equal(nrow(res_f$catt_df), 1L)
	expect_true(is.finite(res_f$att_hat))
	expect_true(is.finite(res_f$att_se))
	expect_gt(res_f$att_se, 0)
	# True effect tau = 2; bridge shrinks slightly toward 0 but stays close.
	expect_equal(res_f$att_hat, 2, tolerance = 0.3)

	res_e <- .s251_fit(p22, etwfe)
	expect_s3_class(res_e, "etwfe")
	expect_equal(res_e$T, 2)
	expect_equal(res_e$att_hat, 2, tolerance = 0.3)

	res_b <- .s251_fit(p22, betwfe)
	expect_s3_class(res_b, "betwfe")
	expect_equal(res_b$T, 2)
	expect_equal(res_b$att_hat, 2, tolerance = 0.3)

	res_t <- .s251_fit(p22, twfeCovs)
	expect_s3_class(res_t, "twfeCovs")
	expect_equal(res_t$T, 2)
	expect_true(is.finite(res_t$att_hat))
	expect_true(is.finite(res_t$att_se))
})

test_that("fetwfe add_ridge fits a 2x2 panel (the relaxed T >= 2L fusion guard)", {
	skip_on_cran()
	# The add_ridge path builds genFullInvFusionTransformMat, whose `T >= 3L`
	# stopifnot used to reject T = 2; #251 relaxed it to `T >= 2L`.
	p22 <- .s251_make_panel(N = 240, T = 2, adopt = 2, tau = 2, seed = 2203)
	res <- suppressWarnings(.s251_fit(p22, fetwfe, add_ridge = TRUE))
	expect_s3_class(res, "fetwfe")
	expect_equal(res$T, 2)
	expect_true(is.finite(res$att_hat))
	expect_equal(res$att_hat, 2, tolerance = 0.3)
})

test_that("accessors and SE paths work at T = 2", {
	skip_on_cran()
	p22 <- .s251_make_panel(N = 240, T = 2, adopt = 2, tau = 2, seed = 2204)
	res <- suppressWarnings(.s251_fit(p22, fetwfe))

	es <- eventStudy(res)
	expect_s3_class(es, "eventStudy")
	expect_equal(nrow(es), 1L)
	expect_true(all(es$n_cohorts == 1L))
	expect_true(all(is.finite(es$estimate)))

	cs <- cohortStudy(res)
	expect_s3_class(cs, "cohortStudy")
	expect_equal(nrow(cs), 1L)

	# cohort family (K = 1, the critical-value bypass) and the K-row families all
	# produce valid simultaneous CIs with a single treatment effect.
	sci_cohort <- simultaneousCIs(res, family = "cohort")
	expect_s3_class(sci_cohort, "simultaneous_cis")
	expect_equal(sci_cohort$K, 1L)

	sci_es <- simultaneousCIs(res, family = "event_study")
	expect_s3_class(sci_es, "simultaneous_cis")
	expect_equal(nrow(sci_es$ci), 1L)
	expect_true(all(is.finite(sci_es$ci$simultaneous_ci_low)))
	expect_true(all(is.finite(sci_es$ci$simultaneous_ci_high)))

	sci_all <- simultaneousCIs(res, family = "all_post_treatment")
	expect_s3_class(sci_all, "simultaneous_cis")

	res_cons <- suppressWarnings(.s251_fit(
		p22,
		fetwfe,
		se_type = "conservative"
	))
	expect_true(is.finite(res_cons$att_se))
	expect_gt(res_cons$att_se, 0)

	res_clus <- suppressWarnings(.s251_fit(p22, fetwfe, se_type = "cluster"))
	expect_true(is.finite(res_clus$att_se))
	expect_gt(res_clus$att_se, 0)
})

# ------------------------------------------------------------------------------
# 1b. Late-adopting single cohort (G = 1, T = 5, adopt final period t = 5):
# num_treats = 1 but T >= 3, so the time FEs still fuse. All four estimators fit
# and recover the injected ATT.
# ------------------------------------------------------------------------------

test_that("all four estimators fit a late-adopting (G = 1, T = 5) num_treats = 1 panel", {
	skip_on_cran()
	plate <- .s251_make_panel(N = 240, T = 5, adopt = 5, tau = 2, seed = 2205)

	res_f <- .s251_fit(plate, fetwfe)
	expect_s3_class(res_f, "fetwfe")
	expect_equal(res_f$G, 1)
	expect_equal(res_f$T, 5)
	# Adopting in the final period -> exactly one post-treatment period.
	expect_equal(nrow(res_f$catt_df), 1L)
	expect_true(is.finite(res_f$att_hat))
	expect_true(is.finite(res_f$att_se))
	expect_gt(res_f$att_se, 0)
	expect_equal(res_f$att_hat, 2, tolerance = 0.4)

	res_e <- .s251_fit(plate, etwfe)
	expect_s3_class(res_e, "etwfe")
	expect_equal(res_e$att_hat, 2, tolerance = 0.4)

	res_b <- .s251_fit(plate, betwfe)
	expect_s3_class(res_b, "betwfe")
	expect_equal(res_b$att_hat, 2, tolerance = 0.4)

	res_t <- .s251_fit(plate, twfeCovs)
	expect_s3_class(res_t, "twfeCovs")
	expect_true(is.finite(res_t$att_hat))
})

# ------------------------------------------------------------------------------
# 2. The warning asymmetry (the whole point of #251's warning design).
# ------------------------------------------------------------------------------

test_that("fetwfe() warns about no fusion at T = 2", {
	# MUTATION TARGET: deleting the `if (T == 2) warning(...)` block in fetwfe()
	# makes this test FAIL (the fit succeeds silently). Verified by the
	# implementer (recorded in the #251 PR).
	p22 <- .s251_make_panel(N = 120, T = 2, adopt = 2, tau = 0, seed = 2210)
	expect_warning(
		.s251_fit(p22, fetwfe),
		"no fusion",
		ignore.case = TRUE
	)
})

test_that("fetwfe() does NOT warn at a late-adopting num_treats = 1 case (T >= 3)", {
	# T = 5, adopt at t = 5 -> num_treats = 1, but the time FEs (T - 1 = 4 of
	# them) still fuse, so there IS fusion -> no warning.
	plate <- .s251_make_panel(N = 120, T = 5, adopt = 5, tau = 0, seed = 2211)
	expect_warning(.s251_fit(plate, fetwfe), NA)
})

test_that("fetwfe() does NOT warn at num_treats >= 2", {
	# Standard multi-period single cohort (T = 6, adopt 2 -> num_treats = 5).
	pmulti <- .s251_make_panel(N = 120, T = 6, adopt = 2, tau = 0, seed = 2212)
	expect_warning(.s251_fit(pmulti, fetwfe), NA)
})

test_that("betwfe() and etwfe() do NOT warn at T = 2", {
	# betwfe applies the bridge to the RAW coefficients (never fuses), and etwfe
	# is unpenalized, so neither emits the no-fusion warning -- it is fetwfe-only.
	p22 <- .s251_make_panel(N = 120, T = 2, adopt = 2, tau = 0, seed = 2213)
	expect_warning(.s251_fit(p22, betwfe), NA)
	expect_warning(.s251_fit(p22, etwfe), NA)
	expect_warning(.s251_fit(p22, twfeCovs), NA)
})

# ------------------------------------------------------------------------------
# 3. T-sweep continuity: fetwfe() runs at every T in 2..6 (G = 1, same DGP
# shape), and ONLY the T = 2 step warns. This is the continuity motivation for
# #251 -- a simulation sweep over T should not break (or spuriously warn) at the
# T = 2 endpoint.
# ------------------------------------------------------------------------------

test_that("fetwfe() runs across a T = 2..6 sweep, warning only at T = 2", {
	skip_on_cran()
	for (Tval in 2:6) {
		pdata <- .s251_make_panel(
			N = 160,
			T = Tval,
			adopt = 2,
			tau = 1,
			seed = 2220 + Tval
		)
		warned <- FALSE
		res <- withCallingHandlers(
			.s251_fit(pdata, fetwfe),
			warning = function(cnd) {
				if (
					grepl(
						"no fusion",
						conditionMessage(cnd),
						ignore.case = TRUE
					)
				) {
					warned <<- TRUE
				}
				invokeRestart("muffleWarning")
			}
		)
		expect_s3_class(res, "fetwfe")
		expect_equal(res$T, Tval)
		expect_true(is.finite(res$att_hat))
		# The no-fusion warning fires iff T == 2.
		expect_equal(warned, Tval == 2L)
	}
})

# ------------------------------------------------------------------------------
# 4. genCoefs / simulateData round-trip at T = 2 (the DGP side of #251).
# ------------------------------------------------------------------------------

test_that("genCoefs and simulateData round-trip at T = 2", {
	# seed = 251 yields a non-zero single treatment effect (att_true = 2), so the
	# att_true / att_hat assertions below exercise a real effect rather than a
	# tautological 0 ~ 0. At density 0.5 the lone num_treats = 1 effect is
	# sparsified out for some seeds (e.g. 2230 / 99 -> att_true = 0; cf. the
	# all-zero-block risk in the project memory); density must be strictly in
	# (0, 1), so density = 1 is not an option.
	coefs <- genCoefs(
		G = 1,
		T = 2,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 251
	)
	expect_s3_class(coefs, "FETWFE_coefs")
	expect_equal(coefs$G, 1)
	expect_equal(coefs$T, 2)

	tes <- getTes(coefs)
	expect_length(tes$actual_cohort_tes, 1L)
	# Overall ATT equals the lone cohort's single-period effect.
	expect_equal(tes$att_true, tes$actual_cohort_tes[1])

	sim <- simulateData(coefs, N = 400, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	expect_equal(sim$G, 1)
	expect_equal(sim$T, 2)
	expect_s3_class(sim$pdata, "data.frame")
	# One treated cohort -> exactly one unique post-adoption start time.
	treated <- sim$pdata[sim$pdata$treatment == 1L, , drop = FALSE]
	first_treat <- stats::aggregate(time ~ unit, data = treated, FUN = min)
	expect_equal(length(unique(first_treat$time)), 1L)

	# The simulated 2x2 fits and recovers the injected ATT (numerical sanity).
	fit <- suppressWarnings(fetwfe(
		pdata = sim$pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = sim$covs,
		response = "y",
		verbose = FALSE
	))
	expect_equal(fit$T, 2)
	expect_true(is.finite(fit$att_hat))
	expect_equal(fit$att_hat, tes$att_true, tolerance = 0.4)
})
