# Tests for getTes()'s `actual_event_time_tes` accessor and the internal
# `.true_event_time_effects()` helper (#388). The load-bearing test binds the
# helper to the estimator's `eventStudy()` estimand. As of #389 the two estimator
# event-study loops (`.event_study_fetwfe()` / `.event_study_etwfe_betwfe()`) are
# folded onto this helper (single source of truth); the byte-exact test below
# locks that `eventStudy()`'s estimand IS the helper's output for all three
# estimators.

# Helper: the event-time effects the estimator targets, computed from a fit's
# estimated cells + resolved cohort structure.
.helper_est_from_fit <- function(fit) {
	offs <- fetwfe:::.resolve_cohort_offsets_and_first_inds(fit, fit$G, fit$T)
	fetwfe:::.true_event_time_effects(
		cell_effects = fit$beta_hat[fit$treat_inds],
		first_inds = offs$first_inds,
		cohort_offsets_int = offs$cohort_offsets_int,
		cohort_probs = fit$cohort_probs_overall,
		T = fit$T
	)
}

test_that(".true_event_time_effects reproduces eventStudy()'s estimand (marginal)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 3
	)
	sim <- simulateData(
		coefs,
		N = 400,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 3
	)
	fit <- fetwfeWithSimulatedData(sim)
	expect_equal(
		unname(.helper_est_from_fit(fit)),
		eventStudy(fit)$estimate,
		tolerance = 1e-10
	)
})

test_that("getTes()'s assumed indices equal a fit's resolved indices (genCoefs)", {
	for (cfg in list(c(3L, 5L), c(2L, 6L), c(4L, 7L))) {
		G <- cfg[1]
		T <- cfg[2]
		coefs <- genCoefs(
			G = G,
			T = T,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = 5
		)
		sim <- simulateData(
			coefs,
			N = 300,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 5
		)
		fit <- fetwfeWithSimulatedData(sim)
		offs <- fetwfe:::.resolve_cohort_offsets_and_first_inds(fit, G, T)
		expect_equal(offs$first_inds, fetwfe:::getFirstInds(G = G, T = T))
		expect_equal(offs$cohort_offsets_int, seq.int(2L, G + 1L))
	}
})

test_that("getTes()$actual_event_time_tes: shape, names, and helper consistency", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 9
	)
	te <- getTes(coefs)
	expect_length(te$actual_event_time_tes, coefs$T - 1L)
	expect_named(te$actual_event_time_tes, as.character(0:(coefs$T - 2L)))
	num_treats <- fetwfe:::getNumTreats(G = coefs$G, T = coefs$T)
	treat_inds <- fetwfe:::getTreatInds(
		G = coefs$G,
		T = coefs$T,
		d = coefs$d,
		num_treats = num_treats
	)
	expect_equal(
		te$actual_event_time_tes,
		fetwfe:::.true_event_time_effects(
			cell_effects = coefs$beta[treat_inds],
			first_inds = fetwfe:::getFirstInds(G = coefs$G, T = coefs$T),
			cohort_offsets_int = seq.int(2L, coefs$G + 1L),
			cohort_probs = te$cohort_weights,
			T = coefs$T
		)
	)
})

test_that("event-time truth matches the estimand under a multinomial DGP", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		seed = 11
	)
	sim <- simulateData(
		coefs,
		N = 400,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 11
	)
	fit <- fetwfeWithSimulatedData(sim)
	expect_equal(
		unname(.helper_est_from_fit(fit)),
		eventStudy(fit)$estimate,
		tolerance = 1e-10
	)
})

test_that("event-time truth matches the estimand under event_study fusion", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.25,
		eff_size = 2,
		fusion_structure = "event_study",
		seed = 4
	)
	sim <- simulateData(
		coefs,
		N = 400,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 4
	)
	fit <- fetwfeWithSimulatedData(sim)
	expect_equal(
		unname(.helper_est_from_fit(fit)),
		eventStudy(fit)$estimate,
		tolerance = 1e-10
	)
})

test_that(".true_event_time_effects returns NA at an event time no cohort reaches", {
	# Synthetic offsets 5,6 with T=6: V_e = which(offsets <= 6 - e).
	# e=0 -> {1,2}; e=1 -> {1}; e=2 -> {} (offsets <= 4) -> NA; likewise e=3,4.
	res <- fetwfe:::.true_event_time_effects(
		cell_effects = rep(1, 20),
		first_inds = c(1L, 6L),
		cohort_offsets_int = c(5L, 6L),
		cohort_probs = c(0.5, 0.5),
		T = 6L
	)
	expect_length(res, 5L)
	expect_named(res, as.character(0:4))
	expect_false(is.na(res[["0"]]))
	expect_true(is.na(res[["2"]]))
	expect_true(is.na(res[["4"]]))
})

test_that(".true_event_time_effects() single-sources eventStudy()'s estimand byte-exactly (fetwfe/etwfe/betwfe, #389)", {
	# After #389 both estimator loops compute their point estimates by calling
	# `.true_event_time_effects()` (with `empty_value = 0`), using the SAME
	# sparse-weight crossprod form, so `eventStudy(fit)$estimate` equals the
	# helper's output BIT-FOR-BIT. Loop seeds 1-6: whether the old inline `sum`
	# form diverged from `crossprod` is BLAS/DGP-dependent (e.g. betwfe seed 3
	# does not diverge on some platforms), so a single pinned seed could make an
	# estimator's assertion vacuous; looping guarantees the mutation-check bites.
	for (seed in 1:6) {
		coefs <- genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = seed
		)
		sim <- simulateData(
			coefs,
			N = 400,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = seed
		)
		fits <- list(
			fetwfeWithSimulatedData(sim),
			etwfeWithSimulatedData(sim),
			betwfeWithSimulatedData(sim)
		)
		for (fit in fits) {
			offs <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
				fit,
				fit$G,
				fit$T
			)
			est_via_helper <- unname(fetwfe:::.true_event_time_effects(
				cell_effects = fit$beta_hat[fit$treat_inds],
				first_inds = offs$first_inds,
				cohort_offsets_int = offs$cohort_offsets_int,
				cohort_probs = fit$cohort_probs_overall,
				T = fit$T,
				empty_value = 0
			))
			# expect_identical (not expect_equal): byte-for-byte is the point.
			# unname() is load-bearing -- a named vector would flip the data
			# frame's row.names, failing df-level identical().
			expect_identical(est_via_helper, eventStudy(fit)$estimate)
		}
	}
})

test_that(".true_event_time_effects() `empty_value` controls the unreached-event-time fill (#389)", {
	# Synthetic offsets 5,6 with T=6: event times 2,3,4 are reached by no cohort
	# (V_e empty). The default fills NA (truth convention); the estimator loops
	# pass empty_value = 0 (a zero pooled estimate). Reached rows are unaffected.
	args <- list(
		cell_effects = rep(1, 20),
		first_inds = c(1L, 6L),
		cohort_offsets_int = c(5L, 6L),
		cohort_probs = c(0.5, 0.5),
		T = 6L
	)
	res_na <- do.call(fetwfe:::.true_event_time_effects, args)
	res_0 <- do.call(
		fetwfe:::.true_event_time_effects,
		c(args, list(empty_value = 0))
	)
	expect_true(is.na(res_na[["2"]]))
	expect_true(is.na(res_na[["4"]]))
	expect_identical(res_0[["2"]], 0)
	expect_identical(res_0[["4"]], 0)
	expect_false(is.na(res_0[["0"]]))
	expect_identical(res_na[["0"]], res_0[["0"]])
})
