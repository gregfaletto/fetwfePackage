library(testthat)
library(fetwfe)

# Issue #238: NEW tests exercising the estimator and its accessors against
# NATIVELY generated event-study-sparse truth via
# genCoefs(fusion_structure = "event_study") (added in #237), complementing the
# hand-built-truth tests in test-event-study-fusion-40.R (left untouched as a
# regression anchor). The shared forward-differencing helper .es_setup() lives
# in tests/testthat/helper-es-fusion.R.
#
# Empirical note (density- AND shape-dependence): on native event-study-sparse
# truth the event-study penalty's MSE advantage over the cohort penalty is
# reliable only at favorable shapes (large treatment block + covariates) and low
# density. At the vignette shape (G=6,T=9,d=1,N=400,density=0.25) it wins on the
# majority of seeds (median ratio ~0.47); at d=0 / smaller shapes the MSE race is
# a coin flip or worse (e.g. (5,8,0) median ~1.5). Per-cell RECOVERY (correlation
# with the truth) is excellent (~1.0) regardless. So the "beats cohort" claim is
# asserted only at the favorable shape; the d=0 transform path is covered by
# recovery.

# Per-seed MSE ratio (event_study / cohort) for native event-study-sparse truth;
# NA for a degenerate all-zero treatment block.
.es_recovery_ratios <- function(G, T, d, N, seeds, density = 0.25) {
	num_treats <- fetwfe:::getNumTreats(G = G, T = T)
	treat_inds <- fetwfe:::getTreatInds(
		G = G,
		T = T,
		d = d,
		num_treats = num_treats
	)
	vapply(
		seeds,
		function(sd) {
			co <- genCoefs(
				G = G,
				T = T,
				d = d,
				density = density,
				eff_size = 2,
				seed = sd,
				fusion_structure = "event_study"
			)
			tb <- co$beta[treat_inds]
			if (all(tb == 0)) {
				return(NA_real_)
			}
			sim <- simulateData(
				co,
				N = N,
				sig_eps_sq = 1,
				sig_eps_c_sq = 0.5,
				seed = sd
			)
			fd <- suppressWarnings(fetwfeWithSimulatedData(
				sim,
				lambda_selection = "bic"
			))
			fe <- suppressWarnings(fetwfeWithSimulatedData(
				sim,
				lambda_selection = "bic",
				fusion_structure = "event_study"
			))
			mean((fe$beta_hat[treat_inds] - tb)^2) /
				mean((fd$beta_hat[treat_inds] - tb)^2)
		},
		numeric(1)
	)
}

test_that("event_study beats cohort on native ES truth at a favorable shape (#238)", {
	skip_on_cran()
	# Headline recovery claim, asserted where it is robust: large treatment block
	# + covariates, low density. median(ratio) < 1 means event_study wins on the
	# majority of a fixed (non-cherry-picked) seed set.
	ratios <- .es_recovery_ratios(G = 6, T = 9, d = 1, N = 400, seeds = 1:12)
	expect_gte(sum(is.finite(ratios)), 10L) # not mass-degenerate
	expect_lt(stats::median(ratios, na.rm = TRUE), 1)
})

test_that("event_study fit recovers native ES truth on both transform paths (#238)", {
	skip_on_cran()
	# Per-cell recovery (correlation of beta_hat with the native truth) is robust
	# (~1.0) even where the MSE race is close. Covers BOTH inverse-transform sites
	# #237 routes through .gen_inv_treat_block(): d = 0 (base treatment effects
	# only) and d > 0 (also the treatment x covariate interaction site).
	recover_cor <- function(G, T, d, N, sd) {
		num_treats <- fetwfe:::getNumTreats(G = G, T = T)
		treat_inds <- fetwfe:::getTreatInds(
			G = G,
			T = T,
			d = d,
			num_treats = num_treats
		)
		co <- genCoefs(
			G = G,
			T = T,
			d = d,
			density = 0.25,
			eff_size = 2,
			seed = sd,
			fusion_structure = "event_study"
		)
		sim <- simulateData(
			co,
			N = N,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = sd
		)
		fe <- suppressWarnings(fetwfeWithSimulatedData(
			sim,
			lambda_selection = "bic",
			fusion_structure = "event_study"
		))
		stats::cor(fe$beta_hat[treat_inds], co$beta[treat_inds])
	}
	expect_gt(recover_cor(6, 9, 1, 400, 10), 0.9) # d > 0 path
	expect_gt(recover_cor(4, 7, 0, 300, 7), 0.9) # d = 0 path
})

test_that("accessors recover and cover native ES truth (#238)", {
	skip_on_cran()
	# eventStudy / cohortStudy / simultaneousCIs on a native event-study-sparse
	# fit (shape (6,9,1), seed 10 -- the #237 vignette seed). Behavioral
	# assertions (recovery / coverage / structural ordering), not is.finite-only.
	G <- 6L
	T <- 9L
	d <- 1L
	N <- 400L
	# Treatment-block layout from the shared helper (tests/testthat/helper-es-fusion.R).
	s <- .es_setup(G, T)
	num_treats <- s$num_treats
	first_inds <- s$first_inds
	n_k <- s$n_k
	treat_inds <- fetwfe:::getTreatInds(
		G = G,
		T = T,
		d = d,
		num_treats = num_treats
	)
	event_time <- unlist(lapply(n_k, function(nk) 0:(nk - 1)))

	co <- genCoefs(
		G = G,
		T = T,
		d = d,
		density = 0.25,
		eff_size = 2,
		seed = 10,
		fusion_structure = "event_study"
	)
	true_te <- co$beta[treat_inds]
	true_prof <- as.numeric(tapply(true_te, event_time, mean))
	true_cohort <- vapply(
		seq_len(G),
		function(g) {
			mean(true_te[fetwfe:::.cohort_block_inds(
				g,
				G,
				first_inds,
				num_treats
			)])
		},
		numeric(1)
	)

	sim <- simulateData(
		co,
		N = N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 10
	)
	fit <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic",
		fusion_structure = "event_study"
	))

	# eventStudy: aggregated per-event-time profile tracks the true profile.
	es <- suppressMessages(eventStudy(fit))
	expect_equal(nrow(es), length(true_prof))
	expect_gt(stats::cor(es$estimate, true_prof), 0.9)

	# cohortStudy: per-cohort ATTs track the true per-cohort means.
	cs <- suppressMessages(cohortStudy(fit))
	expect_equal(nrow(cs), G)
	expect_gt(stats::cor(cs$estimate, true_cohort), 0.9)

	# simultaneousCIs (event_study family): structural ordering + covers truth.
	sci <- suppressMessages(simultaneousCIs(
		fit,
		family = "event_study",
		alpha = 0.05
	))
	# Simultaneous band is wider than pointwise but no wider than Bonferroni.
	expect_lt(sci$pointwise_critical_value, sci$critical_value)
	expect_lte(sci$critical_value, sci$bonferroni_critical_value)
	ci <- sci$ci
	sim_w <- ci$simultaneous_ci_high - ci$simultaneous_ci_low
	pt_w <- ci$pointwise_ci_high - ci$pointwise_ci_low
	expect_true(all(sim_w >= pt_w - 1e-9))
	# The simultaneous band covers the (native) true event-time profile.
	covered <- true_prof >= ci$simultaneous_ci_low &
		true_prof <= ci$simultaneous_ci_high
	expect_gte(sum(covered), length(true_prof) - 1L)
})
