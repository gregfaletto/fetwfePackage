# Tests for the method-entry preconditions introduced in #86.
#
# Two main goals:
#   (1) #73 regression: eventStudy() returns all-NA SEs when the
#       fitted object's calc_ses is FALSE (q in {1, 2} for fetwfe/betwfe).
#   (2) Cross-class method-entry checks: each .check_for_<method>()
#       helper accepts well-formed input, returns the expected shape,
#       and is wired into the corresponding method's entry point.

# ------------------------------------------------------------------------------
# Shared fixture
# ------------------------------------------------------------------------------

.precond_fixture <- function(q = 0.5, seed = 1) {
	sim <- simulateData(
		genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = seed),
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	list(
		sim = sim,
		fetwfe = suppressWarnings(fetwfeWithSimulatedData(sim, q = q)),
		etwfe = suppressWarnings(etwfeWithSimulatedData(sim)),
		betwfe = suppressWarnings(betwfeWithSimulatedData(sim, q = q))
	)
}

# ------------------------------------------------------------------------------
# #73 regression: eventStudy() respects the fit's calc_ses contract
# ------------------------------------------------------------------------------

test_that("eventStudy(fetwfe(q=1)) returns all-NA SEs and p-values (fixes #73)", {
	fx <- .precond_fixture(q = 1)
	# Sanity: the fit's own contract says no SEs.
	expect_true(is.na(fx$fetwfe$att_se))
	expect_false(isTRUE(fx$fetwfe$internal$calc_ses))
	# eventStudy must respect that.
	es <- eventStudy(fx$fetwfe)
	expect_true(all(is.na(es$se)))
	expect_true(all(is.na(es$p_value)))
	# Point estimates ARE valid even without SEs.
	expect_true(any(is.finite(es$estimate)))
})

test_that("eventStudy(fetwfe(q=2)) returns all-NA SEs (ridge, fixes #73)", {
	fx <- .precond_fixture(q = 2)
	expect_true(is.na(fx$fetwfe$att_se))
	expect_false(isTRUE(fx$fetwfe$internal$calc_ses))
	es <- eventStudy(fx$fetwfe)
	expect_true(all(is.na(es$se)))
	expect_true(all(is.na(es$p_value)))
})

test_that("eventStudy(betwfe(q=1)) returns all-NA SEs (fixes #73 for BETWFE)", {
	fx <- .precond_fixture(q = 1)
	expect_true(is.na(fx$betwfe$att_se))
	expect_false(isTRUE(fx$betwfe$calc_ses))
	es <- eventStudy(fx$betwfe)
	expect_true(all(is.na(es$se)))
	expect_true(all(is.na(es$p_value)))
})

test_that("eventStudy(betwfe(q=2)) returns all-NA SEs (ridge, fixes #73)", {
	fx <- .precond_fixture(q = 2)
	expect_true(is.na(fx$betwfe$att_se))
	es <- eventStudy(fx$betwfe)
	expect_true(all(is.na(es$se)))
})

test_that("eventStudy still produces finite SEs on the valid paths (regression check)", {
	fx <- .precond_fixture(q = 0.5)
	# fetwfe q=0.5 has valid SEs.
	expect_true(isTRUE(fx$fetwfe$internal$calc_ses))
	es_f <- eventStudy(fx$fetwfe)
	expect_true(any(is.finite(es_f$se)))
	# etwfe always has valid SEs (pure OLS).
	expect_true(isTRUE(fx$etwfe$calc_ses))
	es_e <- eventStudy(fx$etwfe)
	expect_true(any(is.finite(es_e$se)))
	# betwfe q=0.5 has valid SEs.
	expect_true(isTRUE(fx$betwfe$calc_ses))
	es_b <- eventStudy(fx$betwfe)
	expect_true(any(is.finite(es_b$se)))
})

# ------------------------------------------------------------------------------
# Each .check_for_<method> accepts well-formed input
# ------------------------------------------------------------------------------

test_that(".check_for_event_study returns has_valid_ses derived from calc_ses path", {
	fx05 <- .precond_fixture(q = 0.5)
	contract_f <- fetwfe:::.check_for_event_study(fx05$fetwfe)
	expect_true(contract_f$has_valid_ses)
	contract_e <- fetwfe:::.check_for_event_study(fx05$etwfe)
	expect_true(contract_e$has_valid_ses)
	contract_b <- fetwfe:::.check_for_event_study(fx05$betwfe)
	expect_true(contract_b$has_valid_ses)

	fx1 <- .precond_fixture(q = 1)
	expect_false(fetwfe:::.check_for_event_study(fx1$fetwfe)$has_valid_ses)
	expect_false(fetwfe:::.check_for_event_study(fx1$betwfe)$has_valid_ses)
})

test_that(".check_for_augment / tidy / glance / plot / coef accept all 3 estimators", {
	fx <- .precond_fixture(q = 0.5)
	for (res in list(fx$fetwfe, fx$etwfe, fx$betwfe)) {
		expect_silent(fetwfe:::.check_for_augment(res))
		expect_silent(fetwfe:::.check_for_tidy(res))
		expect_silent(fetwfe:::.check_for_glance(res))
		expect_silent(fetwfe:::.check_for_plot(res))
		expect_silent(fetwfe:::.check_for_coef(res))
	}
})

# ------------------------------------------------------------------------------
# Each helper rejects non-estimator input with a clear message
# ------------------------------------------------------------------------------

test_that("preconditions reject non-estimator objects with a clear message", {
	not_estimator <- list(att_hat = 0.5, foo = "bar")
	expect_error(
		fetwfe:::.assert_estimator_object(not_estimator),
		"Expected a `fetwfe`, `etwfe`, `betwfe`, or `twfeCovs` object",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.check_for_event_study(not_estimator),
		"Expected a `fetwfe`, `etwfe`, `betwfe`, or `twfeCovs` object",
		fixed = TRUE
	)
	expect_error(
		fetwfe:::.check_for_augment(not_estimator),
		"Expected a `fetwfe`, `etwfe`, `betwfe`, or `twfeCovs` object",
		fixed = TRUE
	)
})

# ------------------------------------------------------------------------------
# Mutation: hand-modify a class object to violate the C1 contract; the
# precondition (via .validate_<class>) catches it
# ------------------------------------------------------------------------------

test_that("precondition catches hand-modified objects that violate C1 (att_se NA but calc_ses FALSE)", {
	fx <- .precond_fixture(q = 0.5)
	res <- fx$fetwfe
	# Mutate: set calc_ses=FALSE (now C1 requires att_se to be NA, but
	# the original fit had a finite att_se).
	res$internal$calc_ses <- FALSE
	expect_true(is.finite(res$att_se)) # mutation set up the inconsistency
	expect_error(
		fetwfe:::.check_for_event_study(res),
		"C1 SE consistency",
		fixed = TRUE
	)
})

# ------------------------------------------------------------------------------
# Cluster-SE × q=1: contract gate cascades through the cluster path too
# ------------------------------------------------------------------------------

test_that("eventStudy respects calc_ses under se_type='cluster' × q=1", {
	# Use a panel large enough for the cluster path to be well-conditioned.
	sim <- simulateData(
		genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 1),
		N = 200,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- suppressWarnings(
		fetwfeWithSimulatedData(sim, q = 1, se_type = "cluster")
	)
	# Sanity: q=1 still means calc_ses=FALSE regardless of se_type.
	expect_true(is.na(res$att_se))
	expect_false(isTRUE(res$internal$calc_ses))
	es <- eventStudy(res)
	expect_true(all(is.na(es$se)))
})
