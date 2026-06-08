# Tests for cohortTimeATTs() and tidy.cohortTimeATTs() (#280).
#
# The pivotal correctness test cross-checks the per-cell standard error against
# simultaneousCIs(family = "all_post_treatment") --- the *independent* per-cell
# pointwise SE already in the package --- to floating-point precision. (Verified
# bite-y: scaling the var_1 formula in cohort_time_atts.R by a constant, e.g.
# (N * T) -> (N * T * 2), or dropping the 1/(N*T) scaling, blows this comparison
# far past the 1e-8 tolerance and fails the cross-check for all three
# estimators. Note that dropping the sig_eps_sq factor is a no-op here because
# the fixtures fit with sig_eps_sq = 1, so it is not a valid bite check.)

# ---- Shared fixtures (fit once at collection time) --------------------------
.cta_coefs <- genCoefs(
	G = 4,
	T = 6,
	d = 2,
	density = 0.5,
	eff_size = 2,
	seed = 101
)
.cta_dat <- simulateData(
	.cta_coefs,
	N = 160,
	sig_eps_sq = 1,
	sig_eps_c_sq = 0.5,
	seed = 101
)
.cta_fits <- list(
	fetwfe = fetwfeWithSimulatedData(.cta_dat),
	etwfe = etwfeWithSimulatedData(.cta_dat),
	betwfe = betwfeWithSimulatedData(.cta_dat)
)

# A low-density FETWFE fit that reliably fuses cells to exactly zero, for the
# fused-away convention test (7 of 14 cells fused at this seed/density).
.cta_fused_fit <- fetwfeWithSimulatedData(
	simulateData(
		genCoefs(G = 4, T = 6, d = 2, density = 0.2, eff_size = 2, seed = 101),
		N = 160,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 101
	)
)

test_that("cohortTimeATTs shape, class, and columns per estimator", {
	for (nm in names(.cta_fits)) {
		res <- .cta_fits[[nm]]
		cta <- cohortTimeATTs(res)
		num_treats <- length(res$treat_inds)

		expect_s3_class(cta, "cohortTimeATTs")
		expect_identical(class(cta), c("cohortTimeATTs", "data.frame"))
		expect_equal(nrow(cta), num_treats)

		base_cols <- c(
			"cohort",
			"time",
			"estimate",
			"se",
			"ci_low",
			"ci_high",
			"p_value"
		)
		has_sel <- nm %in% c("fetwfe", "betwfe")
		expect_identical(
			names(cta),
			if (has_sel) c(base_cols, "selected") else base_cols
		)
		expect_type(cta$cohort, "character")
		expect_true(is.numeric(cta$time))
	}
})

test_that("cohortTimeATTs estimates equal beta_hat[treat_inds]", {
	for (nm in names(.cta_fits)) {
		res <- .cta_fits[[nm]]
		cta <- cohortTimeATTs(res)
		expect_equal(cta$estimate, unname(res$beta_hat[res$treat_inds]))
	}
})

test_that("time labels match cohort adoption time + event time", {
	for (nm in names(.cta_fits)) {
		res <- .cta_fits[[nm]]
		cta <- cohortTimeATTs(res)
		for (lab in unique(cta$cohort)) {
			tt <- cta$time[cta$cohort == lab]
			# Synthetic fixtures: cohort label is the (numeric) adoption time,
			# and the cohort's cells run at consecutive calendar times from
			# there.
			expect_equal(tt[1], as.numeric(lab))
			if (length(tt) > 1L) {
				expect_equal(diff(tt), rep(1, length(tt) - 1L))
			}
		}
	}
})

test_that("per-cell SE matches simultaneousCIs(all_post_treatment)", {
	for (nm in names(.cta_fits)) {
		res <- .cta_fits[[nm]]
		cta <- cohortTimeATTs(res)
		sci <- simultaneousCIs(
			res,
			family = "all_post_treatment",
			alpha = res$alpha
		)
		# Same cell ordering (precondition for element-wise SE comparison).
		expect_equal(sci$ci$estimate, cta$estimate)
		pc <- sci$pointwise_critical_value
		se_ref <- (sci$ci$pointwise_ci_high - sci$ci$pointwise_ci_low) /
			(2 * pc)
		# Non-vacuity: these q = 0.5 fits have finite SEs.
		expect_true(all(is.finite(cta$se)))
		expect_equal(length(se_ref), length(cta$se))
		expect_lt(max(abs(cta$se - se_ref)), 1e-8)
	}
})

test_that("cluster se_type per-cell SE matches simultaneousCIs", {
	# Exercises the Liang-Zeger sandwich path (se_type = "cluster") for both the
	# theta-space mask (fetwfe) and the beta-space mask (etwfe). simultaneousCIs()
	# recomputes from the same stored se_type, so it remains the ground truth.
	cluster_fits <- list(
		fetwfe = fetwfeWithSimulatedData(.cta_dat, se_type = "cluster"),
		etwfe = etwfeWithSimulatedData(.cta_dat, se_type = "cluster")
	)
	for (nm in names(cluster_fits)) {
		res <- cluster_fits[[nm]]
		expect_identical(res$se_type, "cluster")
		cta <- cohortTimeATTs(res)
		sci <- simultaneousCIs(
			res,
			family = "all_post_treatment",
			alpha = res$alpha
		)
		expect_equal(sci$ci$estimate, cta$estimate)
		pc <- sci$pointwise_critical_value
		se_ref <- (sci$ci$pointwise_ci_high - sci$ci$pointwise_ci_low) /
			(2 * pc)
		expect_true(all(is.finite(cta$se)))
		expect_lt(max(abs(cta$se - se_ref)), 1e-8)
	}
})

test_that("cell estimates aggregate to the per-cohort cohortStudy estimates", {
	for (nm in names(.cta_fits)) {
		res <- .cta_fits[[nm]]
		cta <- cohortTimeATTs(res)
		cs <- cohortStudy(res)
		for (g in seq_len(res$G)) {
			cells_g <- cta$estimate[cta$cohort == cs$cohort[g]]
			expect_equal(mean(cells_g), cs$estimate[g])
		}
	}
})

test_that("CI bounds are pointwise Wald and respond to alpha", {
	res <- .cta_fits$etwfe
	a <- 0.05
	cta <- cohortTimeATTs(res, alpha = a)
	z <- stats::qnorm(1 - a / 2)
	expect_equal(cta$ci_low, cta$estimate - z * cta$se)
	expect_equal(cta$ci_high, cta$estimate + z * cta$se)
	# Wider alpha => narrower interval.
	cta_wide <- cohortTimeATTs(res, alpha = 0.20)
	w05 <- cta$ci_high - cta$ci_low
	w20 <- cta_wide$ci_high - cta_wide$ci_low
	expect_true(all(w20 < w05))
})

test_that("fused-away cells follow the zero/NA convention", {
	cta <- cohortTimeATTs(.cta_fused_fit)
	fused <- cta$se == 0
	expect_gt(sum(fused), 0) # fixture genuinely fuses cells
	# SE is psi-driven, never gated on the estimate: a zero SE is exactly the
	# set of cells the penalty zeroed out.
	expect_true(all(cta$estimate[fused] == 0))
	expect_true(all(cta$ci_low[fused] == 0))
	expect_true(all(cta$ci_high[fused] == 0))
	expect_true(all(is.na(cta$p_value[fused])))
	expect_true(all(!cta$selected[fused]))
	# Surviving cells keep a positive SE and a finite p-value.
	expect_true(all(cta$se[!fused] > 0))
	expect_true(all(is.finite(cta$p_value[!fused])))
	expect_true(all(cta$selected[!fused]))
})

test_that("selected column equals (estimate != 0) for bridge estimators", {
	for (nm in c("fetwfe", "betwfe")) {
		cta <- cohortTimeATTs(.cta_fits[[nm]])
		expect_identical(cta$selected, cta$estimate != 0)
	}
})

test_that("no-SE fit (q >= 1) yields all-NA SE/CI/p but keeps estimates", {
	# Bridge SEs are undefined for q >= 1 (the fit sets calc_ses = FALSE), the
	# contract-valid way the NA-SE path actually arises.
	res_q2 <- fetwfeWithSimulatedData(.cta_dat, q = 2)
	expect_false(isTRUE(res_q2$internal$calc_ses))
	cta <- cohortTimeATTs(res_q2)
	# Estimates are unaffected by the SE gate.
	expect_equal(cta$estimate, unname(res_q2$beta_hat[res_q2$treat_inds]))
	expect_true(all(is.na(cta$se)))
	expect_true(all(is.na(cta$ci_low)))
	expect_true(all(is.na(cta$ci_high)))
	expect_true(all(is.na(cta$p_value)))
})

test_that("fully fused fit (no cells selected) yields all-NA SEs", {
	# When the penalty zeroes the ENTIRE treatment block there is no selected
	# support to recompute the Gram from, so calc_ses is FALSE and every cell
	# gets se/ci/p = NA (not se = 0) --- distinct from the partial-fusion case,
	# and matching eventStudy()'s identical guard. (Documented in @return `se`.)
	all_fused <- fetwfeWithSimulatedData(
		simulateData(
			genCoefs(
				G = 4,
				T = 6,
				d = 2,
				density = 0.05,
				eff_size = 2,
				seed = 1
			),
			N = 160,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 1
		)
	)
	expect_false(isTRUE(all_fused$att_selected)) # fixture genuinely fuses all
	cta <- cohortTimeATTs(all_fused)
	expect_true(all(cta$estimate == 0))
	expect_true(all(is.na(cta$se)))
	expect_true(all(is.na(cta$ci_low)))
	expect_true(all(is.na(cta$ci_high)))
	expect_true(all(is.na(cta$p_value)))
	expect_true(all(!cta$selected))
})

test_that("tidy.cohortTimeATTs renames columns and computes statistic", {
	cta <- cohortTimeATTs(.cta_fits$fetwfe)
	td <- broom::tidy(cta)

	expect_equal(nrow(td), nrow(cta))
	expect_true(all(
		c(
			"term",
			"time",
			"estimate",
			"std.error",
			"statistic",
			"p.value",
			"conf.low",
			"conf.high",
			"selected"
		) %in%
			names(td)
	))
	# term encodes both axes.
	expect_match(td$term[1], "^cohort_.*_time_.*$")
	expect_equal(td$term, paste0("cohort_", cta$cohort, "_time_", cta$time))
	# Pass-through columns.
	expect_equal(td$std.error, cta$se)
	expect_equal(td$conf.low, cta$ci_low)
	expect_equal(td$conf.high, cta$ci_high)
	expect_equal(td$time, cta$time)
	# statistic = estimate / se where se > 0, else NA (independent recompute).
	expect_equal(
		td$statistic,
		ifelse(cta$se > 0, cta$estimate / cta$se, NA_real_)
	)

	# etwfe (no selection) tidies without a `selected` column.
	td_e <- broom::tidy(cohortTimeATTs(.cta_fits$etwfe))
	expect_false("selected" %in% names(td_e))
})

test_that("tidy.cohortTimeATTs guards against a mutated frame", {
	cta <- cohortTimeATTs(.cta_fits$etwfe)
	bad <- cta
	bad$time <- NULL
	expect_error(broom::tidy(bad), "missing required columns")
})

test_that("twfeCovs and non-estimator inputs are rejected helpfully", {
	tc <- twfeCovsWithSimulatedData(.cta_dat)
	expect_error(cohortTimeATTs(tc), "twfeCovs")
	expect_error(cohortTimeATTs(tc), "cohortStudy")

	expect_error(cohortTimeATTs(list(a = 1)), "requires a fitted object")
	expect_error(
		cohortTimeATTs(.cta_fits$etwfe, alpha = 1.5),
		"alpha"
	)
})
