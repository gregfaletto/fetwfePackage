# Tests for getTes()'s `actual_event_time_tes` accessor and the internal
# `.true_event_time_effects()` helper (#388). The load-bearing test binds the
# helper to the estimator's `eventStudy()` estimand. As of #389 the two estimator
# event-study loops (`.event_study_fetwfe()` / `.event_study_etwfe_betwfe()`) are
# folded onto this helper (single source of truth).

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

test_that(".true_event_time_effects() single-sources eventStudy()'s estimand (fetwfe/etwfe/betwfe, #389)", {
	# After #389 both estimator loops compute their point estimates by calling
	# `.true_event_time_effects()` (with `empty_value = 0`), so
	# `eventStudy(fit)$estimate` equals the helper's output.
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
			expect_equal(est_via_helper, eventStudy(fit)$estimate)
		}
	}
})

# The functions that call .finish_event_study(), walked by the #389 blocks below.
event_study_loops <- c(".event_study_etwfe_betwfe", ".event_study_fetwfe")

test_that("the #389 blocks walk every function that calls .finish_event_study() (#389)", {
	calls_finish <- function(fn) {
		found <- FALSE
		for (code in .ns_code_exprs(fn)) {
			.ns_walk_ast(code, function(node, ancestors) {
				if (.ns_is_call_to(node, ".finish_event_study")) {
					found <<- TRUE
				}
			})
		}
		found
	}
	finishers <- names(Filter(calls_finish, .ns_functions("fetwfe")))
	expect_setequal(finishers, event_study_loops)
})

test_that("eventStudy()'s estimator loops compute their estimates by calling .true_event_time_effects() (#389)", {
	for (fn_name in event_study_loops) {
		fn <- get(fn_name, envir = asNamespace("fetwfe"))
		assignments <- list()
		.ns_walk_ast(body(fn), function(node, ancestors) {
			if (
				.ns_is_call_to(node, c("<-", "=", "<<-")) &&
					.ns_subtree_has_symbol(node[[2]], "estimates")
			) {
				assignments[[length(assignments) + 1L]] <<- node
			}
		})
		expect_equal(length(assignments), 1L, info = fn_name)
		calls_helper <- FALSE
		if (length(assignments) == 1L) {
			.ns_walk_ast(assignments[[1]][[3]], function(node, ancestors) {
				if (.ns_is_call_to(node, ".true_event_time_effects")) {
					calls_helper <<- TRUE
				}
			})
		}
		expect_true(calls_helper, info = fn_name)
	}
})

test_that("eventStudy()'s loops hand .finish_event_study() the helper's `estimates` (#389)", {
	for (fn_name in event_study_loops) {
		fn <- get(fn_name, envir = asNamespace("fetwfe"))
		finish_calls <- list()
		.ns_walk_ast(body(fn), function(node, ancestors) {
			if (.ns_is_call_to(node, ".finish_event_study")) {
				finish_calls[[length(finish_calls) + 1L]] <<- node
			}
		})
		expect_equal(length(finish_calls), 1L, info = fn_name)
		passes_estimates <- FALSE
		if (length(finish_calls) == 1L) {
			mc <- match.call(
				get(".finish_event_study", envir = asNamespace("fetwfe")),
				finish_calls[[1]]
			)
			passes_estimates <- identical(mc$estimates, as.name("estimates"))
		}
		expect_true(passes_estimates, info = fn_name)
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
