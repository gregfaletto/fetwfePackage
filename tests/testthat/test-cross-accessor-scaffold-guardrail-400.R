# Cross-accessor byte-identity guardrail for the shared access-time inference
# scaffold (#400, Phase 1 -- the blocking prerequisite for the refactor).
#
# The "resolve selected support -> recompute Gram inverse -> cluster sandwich"
# sequence is duplicated across eventStudy(), cohortTimeATTs(), and
# simultaneousCIs() (R/event_study.R, R/cohort_time_atts.R, R/simultaneous_cis.R).
# If any copy drifts, the accessors silently report INCONSISTENT SEs for the same
# fit. This file pins the three cross-accessor SE/estimate anchors AND literal
# pre-refactor values on a shared fit so the Phase 2 single-sourcing refactor must
# reproduce byte-identical results. There was no unified cross-accessor fixture
# before this (only pairwise checks in test-cohort_time_atts.R / test-simultaneous-cis.R).
#
# simultaneousCIs() exposes no `se` column; the pointwise SE is recovered from the
# pointwise CI half-width. That recovery is algebraically alpha-invariant, so the
# `alpha = fit$alpha` passed below is defensive (matches the fit), not load-bearing.

.g400_dat <- local({
	cf <- genCoefs(G = 4, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 101)
	simulateData(cf, N = 160, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 101)
})

# fetwfe (partial selection, 7/14 cells) + etwfe (no selection) + betwfe (13/14),
# each at se_type default and cluster -> exercises the NA-sentinel, partial-support,
# and cluster-sandwich branches of the scaffold. Built once at collection time.
.g400_fits <- list(
	"fetwfe/default" = fetwfeWithSimulatedData(
		.g400_dat,
		q = 0.5,
		se_type = "default",
		verbose = FALSE
	),
	"etwfe/default" = etwfeWithSimulatedData(
		.g400_dat,
		se_type = "default",
		verbose = FALSE
	),
	"betwfe/default" = betwfeWithSimulatedData(
		.g400_dat,
		se_type = "default",
		verbose = FALSE
	),
	"fetwfe/cluster" = fetwfeWithSimulatedData(
		.g400_dat,
		q = 0.5,
		se_type = "cluster",
		verbose = FALSE
	),
	"etwfe/cluster" = etwfeWithSimulatedData(
		.g400_dat,
		se_type = "cluster",
		verbose = FALSE
	),
	"betwfe/cluster" = betwfeWithSimulatedData(
		.g400_dat,
		se_type = "cluster",
		verbose = FALSE
	)
)

# q = 1 -> has_valid_ses = FALSE: the reachable SE-unavailable path.
.g400_fit_q1 <- fetwfeWithSimulatedData(.g400_dat, q = 1, verbose = FALSE)

.g400_se_of_sci <- function(sci) {
	(sci$ci$pointwise_ci_high - sci$ci$estimate) / sci$pointwise_critical_value
}

# --- Part 1: cross-accessor equality (the first-class regression guard) --------

.g400_expect_anchors <- function(fit, label) {
	a <- fit$alpha

	# Anchor A -- per-(g, t) cell: cohortTimeATTs row <-> all_post_treatment effect.
	# Estimates are the identical contrast (both `tes[idx]`), so exact.
	cta <- cohortTimeATTs(fit)
	ap <- simultaneousCIs(fit, family = "all_post_treatment", alpha = a)
	expect_equal(
		cta$estimate,
		ap$ci$estimate,
		tolerance = 1e-12,
		info = paste(label, "Anchor A estimate")
	)
	expect_lt(max(abs(cta$se - .g400_se_of_sci(ap))), 1e-12)

	# Anchor B -- per-event-time: eventStudy <-> event_study. Compare finite rows
	# only (degenerate event times are NA in eventStudy).
	es <- eventStudy(fit)
	esr <- simultaneousCIs(fit, family = "event_study", alpha = a)
	fin <- is.finite(es$se)
	expect_true(any(fin), info = paste(label, "has a finite event-time SE"))
	expect_equal(
		es$estimate[fin],
		esr$ci$estimate[fin],
		tolerance = 1e-12,
		info = paste(label, "Anchor B estimate")
	)
	expect_lt(max(abs(es$se[fin] - .g400_se_of_sci(esr)[fin])), 1e-12)

	# Anchor C -- per-cohort: cohortStudy (the fit-time scaffold) <-> cohort family.
	# The etwfe estimate side carries ~1e-15 mean-of-block rounding -> tolerance.
	cs <- cohortStudy(fit)
	co <- simultaneousCIs(fit, family = "cohort", alpha = a)
	expect_equal(
		cs$estimate,
		co$ci$estimate,
		tolerance = 1e-9,
		info = paste(label, "Anchor C estimate")
	)
	expect_lt(max(abs(cs$se - .g400_se_of_sci(co))), 1e-12)
}

test_that("the shared scaffold yields cross-accessor-consistent SEs across classes and se_types (#400 guardrail)", {
	for (label in names(.g400_fits)) {
		.g400_expect_anchors(.g400_fits[[label]], label)
	}
})

# --- Part 2: pre-refactor snapshot (the "AND to pre-refactor values" half) ------

test_that("scaffold SEs match pinned pre-refactor values on the fetwfe/default fixture (#400 guardrail)", {
	fit <- .g400_fits[["fetwfe/default"]]
	expect_equal(
		eventStudy(fit)$se,
		c(0.1432551, 0.2040083, 0.2363452, 0.2494955, 0.2219078),
		tolerance = 1e-6
	)
	expect_equal(
		cohortStudy(fit)$se,
		c(0.1738941, 0.1479012, 0.1596468, 0.1438447),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.1817540,
			0.1817540,
			0.1817540,
			0.2219078,
			0.2219078,
			0.1352874,
			0.1840652,
			0.1840652,
			0.2448657,
			0.1352874,
			0.2071185,
			0.2071185,
			0.1352874,
			0.2214705
		),
		tolerance = 1e-6
	)
})

# --- Part 3: SE-unavailable fits degrade consistently ---------------------------

test_that("an SE-unavailable fit degrades consistently: accessors -> NA, simultaneousCIs -> stop (#400 guardrail)", {
	fq <- .g400_fit_q1
	expect_false(isTRUE(fq$internal$calc_ses)) # guard: genuinely has_valid_ses = FALSE

	es_se <- eventStudy(fq)$se
	expect_true(length(es_se) > 0 && all(is.na(es_se)))
	cta_se <- cohortTimeATTs(fq)$se
	expect_true(length(cta_se) > 0 && all(is.na(cta_se)))

	expect_error(simultaneousCIs(fq, family = "event_study"), "calc_ses")
	expect_error(simultaneousCIs(fq, family = "all_post_treatment"), "calc_ses")

	# NB: the *shared-scaffold* singular-Gram branch (recomputed selected-support
	# Gram singular -> res_gram$calc_ses = FALSE; the two accessors degrade to NA
	# while simultaneousCIs STOPs "not invertible") is defensive-only and
	# UNREACHABLE through any public fit -- has_valid_ses = TRUE already implies the
	# fit's own getGramInv() on that support succeeded. Phase 2 guards that
	# on_singular = c("degrade", "stop") policy with a direct unit test of the
	# extracted helper (feeding it a rank-deficient support), which is the only
	# place the asymmetry can be exercised.
})
