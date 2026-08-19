# Cross-accessor byte-identity guardrail for the shared access-time inference
# scaffold (#400, Phase 1 -- the blocking prerequisite for the refactor).
#
# NB the next two sentences are STALE (tracked in #430): the sequence they
# describe has been single-sourced in .recompute_gram_and_sandwich() since
# #413/#414, so there are no longer several copies to drift apart. They are left
# as written -- rewriting them is #430's job, not this file's -- but do not read
# them as a description of the tree. The two sentences after those are still
# accurate: this file does pin the anchors and the literal values, and there
# genuinely was no unified cross-accessor fixture before it.
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

# --- Part 1: cross-accessor equality (no longer the first-class guard) ---------
#
# Measured at a45407b: BOTH mutations inside .recompute_gram_and_sandwich()
# (scaling sandwich_full, scaling gram_inv) leave this entire block green,
# because after the #400 fold both sides of every anchor call that helper -- so
# each anchor now compares the helper against itself. The absolute pins in
# Part 2 are the first-class regression guard now (#429). These anchors are
# kept because they still catch a divergence introduced OUTSIDE the shared
# helper, in an accessor's own assembly of the result.
#
# That residual used to be the whole story, and it is why these anchors were not
# deleted: until #429 PR B they were the ONLY assertions in this file covering
# simultaneousCIs()'s SE path, because the original eighteen Part 2 pins read
# eventStudy(), cohortStudy() and cohortTimeATTs() and none of them went through
# simultaneousCIs() at all. Measured on the #451 review: scaling Sigma one level
# downstream of the scaffold (`Sigma <- 1.5 * (Sigma_1 + Sigma_2)` in
# .simultaneous_cis_impl()) fired 18 failures, every one of them here and none in
# Part 2.
#
# Which set up the same trap one layer up: .assemble_joint_cov_var1/2 are ALREADY
# shared helpers, so the next consolidation routing both sides of these anchors
# through one of them hollows this block exactly as the #400 fold hollowed it
# against the scaffold.
#
# THAT IS NOW PRE-EMPTED (#429 PR B). Each of the six Part 2 blocks carries a
# FOURTH pinned vector, .g400_se_of_sci(simultaneousCIs(fit, family = "cohort")),
# so simultaneousCIs()'s SE assembly has an absolute pin on every fixture --
# twenty-four pins in Part 2, not eighteen. An absolute pin cannot be hollowed by
# consolidating the two things being compared, because there is only one thing.
# Re-measured 2026-08-17 on the finished file: the same Sigma mutation now fires
# 24 failures, 18 in Part 1 and exactly one in each of the six Part 2 blocks.
#
# The fourth vector's VALUE equals the cohortStudy() pin sitting a few lines
# above it in the same block: in ALL SIX blocks the seven printed digits are
# character-identical, which you can confirm without running anything by
# comparing the literal pairs in each block. (They are not bit-identical --
# identical() is FALSE on all six, 8e-17 to 2e-16 apart -- which is why the
# tolerance is 1e-6 and not 0. An earlier draft of this comment said only
# fetwfe/default and fetwfe/cluster matched, implying the other four differed.)
# That is not a defect and not avoidable by changing family (event_study
# reproduces eventStudy()$se, all_post_treatment reproduces
# cohortTimeATTs()$se): the agreement IS what Anchor C asserts. The point of
# pinning it separately is that Anchor C can no longer tell the two paths
# apart, and one absolute value should survive on the simultaneousCIs() side
# when it goes fully hollow.

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

test_that("cross-accessor routing anchors: accessors agree, but are blind to changes inside the shared scaffold (#400/#429)", {
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
	# #429: the fourth pin -- simultaneousCIs()'s OWN SE assembly, recovered by
	# .g400_se_of_sci(). Its value equals the cohortStudy() pin above (here,
	# character-identical at seven digits); that agreement is Anchor C's
	# invariant, and Part 1's header says why pinning it separately is the point.
	# DO NOT "simplify" this to
	# expect_equal(.g400_se_of_sci(...), cohortStudy(fit)$se): that converts an
	# absolute pin back into precisely the parity assertion consolidation
	# hollows. Verified byte-identical on origin/main as well as this branch.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.1738941, 0.1479012, 0.1596468, 0.1438447),
		tolerance = 1e-6
	)
})

# The cluster fixture needs its own ABSOLUTE pins, and the reason is structural.
# Anchor C (Part 1) cross-checks fit-time cohortStudy() against access-time
# simultaneousCIs(). Before the #400 fold those were two independent
# implementations, so that comparison detected drift in either. After the fold
# BOTH sides call .recompute_gram_and_sandwich(), so Anchor C compares the helper
# against itself and a change inside the helper moves both sides together. Proven
# by mutation: perturbing sandwich_full by 5% inside the helper fails Anchor C on
# pre-fold code (3 failures) and passes it entirely on post-fold code.
#
# The general lesson, worth remembering beyond this file: a cross-implementation
# parity test loses ALL of its power the moment the two implementations it
# compares are consolidated. Absolute pins do not have that failure mode.
#
# Part 2 above pins only fetwfe/default, leaving no absolute cluster-sandwich
# value anywhere -- and the one other test that catches such a mutation
# (test-getCohortATTsFinal-unification.R) covers only the OLS branch
# (sel_feat_inds = NULL), so the BRIDGE cluster path had no check at all. Values
# below were computed on this branch and verified byte-identical on origin/main.
test_that("scaffold SEs match pinned pre-refactor values on the fetwfe/cluster fixture (#400 guardrail)", {
	fit <- .g400_fits[["fetwfe/cluster"]]
	expect_equal(
		cohortStudy(fit)$se,
		c(0.1853167, 0.1275883, 0.1469086, 0.1456402),
		tolerance = 1e-6
	)
	expect_equal(
		eventStudy(fit)$se,
		c(0.1466590, 0.2003153, 0.2334116, 0.2379342, 0.2090416),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.1985337,
			0.1985337,
			0.1985337,
			0.2090416,
			0.2090416,
			0.1382049,
			0.1764362,
			0.1764362,
			0.2351856,
			0.1382049,
			0.1930629,
			0.1930629,
			0.1382049,
			0.2308755
		),
		tolerance = 1e-6
	)
	# #429: the fourth pin -- see the fetwfe/default block for the rationale and
	# for the "do not simplify this into a cohortStudy() comparison" prohibition.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.1853167, 0.1275883, 0.1469086, 0.1456402),
		tolerance = 1e-6
	)
})

# The remaining four fixtures, pinned for the same structural reason (#429).
# Before this, Part 2 covered fetwfe/default and fetwfe/cluster only, so a
# mutation inside .recompute_gram_and_sandwich() reddened 3 of what should be 18
# assertions -- Part 1's anchors having gone hollow, nothing else in the file
# could see it. With all six pinned, scaling sandwich_full reddens the three
# `cluster` blocks (the `default` SEs never see it -- the helper returns
# sandwich_full = NULL unless se_type == "cluster") and scaling gram_inv reddens
# the three `default` ones (the `cluster` SEs never see it --
# .compute_cluster_robust_sandwich() solves its own crossprod(X_S_centered) and
# never touches the helper's gram_inv). The two mutations are exact
# complements; together they reach 18 of those 18 assertions. (That count is the
# three eventStudy/cohortStudy/cohortTimeATTs vectors per fixture. #429 PR B
# added a fourth, simultaneousCIs()-based vector to each block, taking Part 2 to
# twenty-four; the two mutations above were measured against the original
# eighteen and have not been re-measured against the fourth.)
#
# Two notes for the reader:
#
#   * Six inline blocks here, while Part 1 loops one helper over all six
#     fixtures. The reason is consistency: the fetwfe/default and fetwfe/cluster
#     blocks were already written this way, and matching them beats mixing two
#     shapes in one section. What matters is that each fixture gets its own name
#     in the reporter, so "how many fixtures did this mutation reach?" is
#     answered by counting failing test names rather than by parsing info=
#     strings off one aggregated test -- but a loop delivers that too, by
#     building the name with paste0(), which is exactly what
#     test-singular-gram-call-site-policy-429.R does. Either shape is fine here;
#     inline is not load-bearing.
#   * betwfe's cohortTimeATTs()$se[4] is EXACTLY zero on both betwfe fixtures.
#     That is a real pinned value, not a placeholder: the bridge selects that
#     cell out (13 of 14 cells selected on this fixture -- see the comment at
#     the top of this file).

test_that("scaffold SEs match pinned pre-refactor values on the etwfe/default fixture (#400 guardrail, #429)", {
	fit <- .g400_fits[["etwfe/default"]]
	expect_equal(
		eventStudy(fit)$se,
		c(0.1639891, 0.2085174, 0.2833031, 0.2903975, 0.3223087),
		tolerance = 1e-6
	)
	expect_equal(
		cohortStudy(fit)$se,
		c(0.2266507, 0.1954749, 0.1990295, 0.2088314),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.2820758,
			0.2887492,
			0.2977382,
			0.3223087,
			0.3223087,
			0.2606832,
			0.2726312,
			0.3009939,
			0.3009939,
			0.2648504,
			0.2952462,
			0.2952462,
			0.2696002,
			0.2696002
		),
		tolerance = 1e-6
	)
	# #429: the fourth pin -- see the fetwfe/default block for the rationale and
	# for the "do not simplify this into a cohortStudy() comparison" prohibition.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.2266507, 0.1954749, 0.1990295, 0.2088314),
		tolerance = 1e-6
	)
})

test_that("scaffold SEs match pinned pre-refactor values on the betwfe/default fixture (#400 guardrail, #429)", {
	fit <- .g400_fits[["betwfe/default"]]
	expect_equal(
		eventStudy(fit)$se,
		c(0.1591995, 0.2037659, 0.2572252, 0.3070540, 0.2826046),
		tolerance = 1e-6
	)
	expect_equal(
		cohortStudy(fit)$se,
		c(0.1490514, 0.1604112, 0.1738290, 0.1978270),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.2515420,
			0.2568852,
			0.2661261,
			0.0000000,
			0.2826046,
			0.2407671,
			0.2529209,
			0.2566070,
			0.2715172,
			0.2542477,
			0.2581047,
			0.2749440,
			0.2484576,
			0.2666360
		),
		tolerance = 1e-6
	)
	# #429: the fourth pin -- see the fetwfe/default block for the rationale and
	# for the "do not simplify this into a cohortStudy() comparison" prohibition.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.1490514, 0.1604112, 0.1738290, 0.1978270),
		tolerance = 1e-6
	)
})

test_that("scaffold SEs match pinned pre-refactor values on the etwfe/cluster fixture (#400 guardrail, #429)", {
	fit <- .g400_fits[["etwfe/cluster"]]
	expect_equal(
		eventStudy(fit)$se,
		c(0.1594857, 0.2023723, 0.2713419, 0.2610870, 0.3509609),
		tolerance = 1e-6
	)
	expect_equal(
		cohortStudy(fit)$se,
		c(0.2396085, 0.1776220, 0.1807051, 0.1994378),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.2793386,
			0.3315840,
			0.2900008,
			0.2623502,
			0.3509609,
			0.2679904,
			0.2633195,
			0.2752361,
			0.2808612,
			0.2747036,
			0.2777359,
			0.2615971,
			0.2568748,
			0.2671456
		),
		tolerance = 1e-6
	)
	# #429: the fourth pin -- see the fetwfe/default block for the rationale and
	# for the "do not simplify this into a cohortStudy() comparison" prohibition.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.2396085, 0.1776220, 0.1807051, 0.1994378),
		tolerance = 1e-6
	)
})

test_that("scaffold SEs match pinned pre-refactor values on the betwfe/cluster fixture (#400 guardrail, #429)", {
	fit <- .g400_fits[["betwfe/cluster"]]
	expect_equal(
		eventStudy(fit)$se,
		c(0.1570818, 0.1998555, 0.2497265, 0.3059113, 0.3083088),
		tolerance = 1e-6
	)
	expect_equal(
		cohortStudy(fit)$se,
		c(0.1648825, 0.1514423, 0.1534372, 0.2044426),
		tolerance = 1e-6
	)
	expect_equal(
		cohortTimeATTs(fit)$se,
		c(
			0.2422136,
			0.2893479,
			0.2629059,
			0.0000000,
			0.3083088,
			0.2450568,
			0.2540328,
			0.2599244,
			0.2661369,
			0.2542525,
			0.2623674,
			0.2549319,
			0.2640755,
			0.2708483
		),
		tolerance = 1e-6
	)
	# #429: the fourth pin -- see the fetwfe/default block for the rationale and
	# for the "do not simplify this into a cohortStudy() comparison" prohibition.
	expect_equal(
		.g400_se_of_sci(
			simultaneousCIs(fit, family = "cohort", alpha = fit$alpha)
		),
		c(0.1648825, 0.1514423, 0.1534372, 0.2044426),
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

	# NB: the fits exercised here degrade because has_valid_ses = FALSE, which is
	# NOT the *shared-scaffold* singular-Gram branch (recomputed selected-support
	# Gram singular -> res_gram$calc_ses = FALSE; the two accessors degrade to NA
	# while simultaneousCIs STOPs "not invertible"). An ordinary fit cannot reach
	# that branch from an ACCESS-TIME caller -- has_valid_ses = TRUE already
	# implies the fit's own getGramInv() on that support succeeded -- but a
	# hand-modified one can, and the fit-time caller getCohortATTsFinal() reaches
	# it outright. For the canonical statement of what is reachable from where,
	# and the tests that pin each route, see the roxygen on
	# .recompute_gram_and_sandwich() in R/variance_machinery.R.
})
