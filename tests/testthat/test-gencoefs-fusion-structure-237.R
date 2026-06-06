library(testthat)
library(fetwfe)

# Issue #237: the fusion_structure argument on genCoefs() / genCoefsCore(), the
# simulation-side companion to #234 (which added the same argument to the
# estimator fetwfe()). These tests pin:
#   (1) the "cohort" default is byte-identical to the pre-#237 output;
#   (2) "event_study" truth is genuinely sparse in the event-study basis;
#   (3) the slot is stored/printed and bad values are rejected;
#   (4) on a representative event-study-sparse draw, event_study fitting
#       recovers the truth better than cohort fitting.
#
# The independent forward differencing .es_difference() / .es_setup() live in
# tests/testthat/helper-es-fusion.R.

test_that("cohort default is byte-identical to the no-arg result (#237)", {
	# Golden vectors captured from the pre-#237 code (main @ 6797a41). Byte-
	# identity of the "cohort" path is also structural: the random draws in
	# genCoefsCore() precede the deterministic inverse transforms, and
	# match.arg() does not consume the RNG, so the default output is unchanged.
	golden_d0 <- c(2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 0, 2, 4, -2, 0)
	golden_d2 <- c(
		2,
		2,
		2,
		2,
		2,
		2,
		0,
		0,
		2,
		2,
		4,
		2,
		4,
		2,
		2,
		0,
		0,
		-2,
		0,
		-2,
		0,
		-2,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		2,
		2,
		0,
		0,
		2,
		2,
		2,
		2,
		2,
		2,
		0,
		0,
		0,
		-2,
		0,
		0,
		2,
		0,
		0,
		2
	)

	expect_equal(
		genCoefs(
			G = 3,
			T = 5,
			d = 0,
			density = 0.5,
			eff_size = 2,
			seed = 11
		)$beta,
		golden_d0,
		tolerance = 1e-10
	)
	expect_equal(
		genCoefs(
			G = 3,
			T = 5,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = 11
		)$beta,
		golden_d2,
		tolerance = 1e-10
	)

	# Default == explicit "cohort", same seed (same code path), across shapes.
	for (cfg in list(c(3, 5, 0), c(4, 6, 1), c(3, 5, 2))) {
		def <- genCoefs(
			G = cfg[1],
			T = cfg[2],
			d = cfg[3],
			density = 0.5,
			eff_size = 2,
			seed = 7
		)$beta
		coh <- genCoefs(
			G = cfg[1],
			T = cfg[2],
			d = cfg[3],
			density = 0.5,
			eff_size = 2,
			seed = 7,
			fusion_structure = "cohort"
		)$beta
		expect_identical(def, coh)
	}

	# Same at the core level.
	expect_identical(
		genCoefsCore(
			G = 4,
			T = 6,
			d = 1,
			density = 0.5,
			eff_size = 2,
			seed = 7
		)$beta,
		genCoefsCore(
			G = 4,
			T = 6,
			d = 1,
			density = 0.5,
			eff_size = 2,
			seed = 7,
			fusion_structure = "cohort"
		)$beta
	)
})

test_that("event_study truth is sparse in the event-study basis (#237)", {
	# The treatment-effect block of an "event_study"-generated beta, pushed
	# through the *independent* forward event-study differencing, must return
	# the stored sparse theta. This fails if the cohort transform were used.
	for (cfg in list(c(3, 6, 0), c(4, 6, 1), c(5, 7, 2))) {
		G <- cfg[1]
		T <- cfg[2]
		d <- cfg[3]
		s <- .es_setup(G, T)
		treat_inds <- fetwfe:::getTreatInds(
			G = G,
			T = T,
			d = d,
			num_treats = s$num_treats
		)
		co <- genCoefs(
			G = G,
			T = T,
			d = d,
			density = 0.5,
			eff_size = 2,
			seed = 3,
			fusion_structure = "event_study"
		)

		theta_es <- .es_difference(co$beta[treat_inds], s$first_inds, s$n_k, G)
		expect_equal(
			theta_es,
			co$theta[treat_inds],
			tolerance = 1e-9,
			info = sprintf("G=%d T=%d d=%d", G, T, d)
		)

		# The recovered representation is genuinely sparse (strictly fewer
		# nonzeros than its length under the density-0.5 draw).
		expect_lt(sum(abs(theta_es) > 1e-8), length(theta_es))
	}
})

test_that("event_study differs from cohort at d=0 and d>0; slot stored/printed; bad value rejected (#237)", {
	for (d in c(0L, 2L)) {
		a <- genCoefs(
			G = 4,
			T = 6,
			d = d,
			density = 0.5,
			eff_size = 2,
			seed = 5
		)
		b <- genCoefs(
			G = 4,
			T = 6,
			d = d,
			density = 0.5,
			eff_size = 2,
			seed = 5,
			fusion_structure = "event_study"
		)
		expect_false(identical(a$beta, b$beta))
		expect_identical(a$fusion_structure, "cohort")
		expect_identical(b$fusion_structure, "event_study")
	}

	# Printed provenance.
	expect_output(
		print(genCoefs(
			G = 3,
			T = 5,
			d = 1,
			density = 0.5,
			eff_size = 2,
			seed = 1,
			fusion_structure = "event_study"
		)),
		"event_study"
	)

	# match.arg rejects an unknown structure at both layers.
	expect_error(genCoefs(
		G = 3,
		T = 5,
		d = 1,
		density = 0.5,
		eff_size = 2,
		fusion_structure = "bogus"
	))
	expect_error(genCoefsCore(
		G = 3,
		T = 5,
		d = 1,
		density = 0.5,
		eff_size = 2,
		fusion_structure = "bogus"
	))

	# Composes with a covariate-dependent assignment DGP (orthogonal options).
	cm <- genCoefs(
		G = 4,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 5,
		fusion_structure = "event_study",
		assignment_type = "multinomial"
	)
	expect_identical(cm$fusion_structure, "event_study")
	expect_identical(cm$assignment_type, "multinomial")
})

test_that("event_study recovers native event-study truth better than cohort (#237)", {
	skip_on_cran()
	# Deterministic, single-seed check pinned to the SAME seed/shape the
	# fusion_structure vignette discloses (seed 10, density 0.25, G=6/T=9/d=1,
	# N=400). On this representative event-study-sparse draw, fitting with the
	# event-study penalty recovers the per-cell treatment effects with lower
	# MSE than the default cohort penalty. The truth is generated NATIVELY via
	# fusion_structure = "event_study" (no hand-built profile). Across-seed
	# robustness (~80% win at this configuration) is documented in the #237
	# plan; this asserts the pinned showcase seed reproducibly wins.
	G <- 6L
	T <- 9L
	d <- 1L
	N <- 400L
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
		seed = 10,
		fusion_structure = "event_study"
	)
	true_te <- co$beta[treat_inds]
	expect_true(any(true_te != 0)) # non-degenerate treatment block

	sim <- simulateData(
		co,
		N = N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 10
	)
	fit_cohort <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic"
	))
	fit_event <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic",
		fusion_structure = "event_study"
	))

	mse_cohort <- mean((fit_cohort$beta_hat[treat_inds] - true_te)^2)
	mse_event <- mean((fit_event$beta_hat[treat_inds] - true_te)^2)
	expect_lt(mse_event, mse_cohort)
})
