library(testthat)
library(fetwfe)

# Issue #40 Phase 1: the event-study fusion transform. These pin the new
# constructor genInvEventStudyFusionTransformMat() against the paper's closed
# form (Faletto 2025, Lemma event.study.sing.val.lem): it must invert the
# event-study differencing exactly, and the implied D must meet the Lemma's
# singular-value bounds sigma_max <= sqrt(6), sigma_min >= 1/(T sqrt(2T)).

# Direct event-study differencing theta = D^{(2)}_ES beta, implemented from the
# paper's *prose* (independent of the inverse's matrix form, so a common error
# is unlikely): cohort 1 contributes its base term tau_{g1,g1} plus
# adjacent-event-time differences; each later cohort k >= 2 contributes
# same-event-time differences against the preceding cohort.
.es_difference <- function(beta, first_inds, n_k, G) {
	theta <- numeric(length(beta))
	c1 <- first_inds[1]:(first_inds[1] + n_k[1] - 1L)
	theta[c1[1]] <- beta[c1[1]]
	if (n_k[1] >= 2L) {
		for (i in 2:n_k[1]) {
			theta[c1[i]] <- beta[c1[i]] - beta[c1[i] - 1L]
		}
	}
	if (G >= 2L) {
		for (k in 2:G) {
			for (e in 0:(n_k[k] - 1L)) {
				theta[first_inds[k] + e] <-
					beta[first_inds[k] + e] - beta[first_inds[k - 1L] + e]
			}
		}
	}
	theta
}

.es_setup <- function(G, T) {
	num_treats <- fetwfe:::getNumTreats(G = G, T = T)
	first_inds <- fetwfe:::getFirstInds(G = G, T = T)
	n_k <- c(diff(first_inds), num_treats - first_inds[G] + 1L)
	list(num_treats = num_treats, first_inds = first_inds, n_k = n_k)
}

test_that("genInvEventStudyFusionTransformMat exactly inverts the event-study differencing (#40)", {
	for (cfg in list(c(3, 6), c(4, 6), c(2, 5), c(5, 7), c(2, 3))) {
		G <- cfg[1]
		T <- cfg[2]
		s <- .es_setup(G, T)
		Dinv <- fetwfe:::genInvEventStudyFusionTransformMat(
			s$num_treats,
			s$first_inds,
			G
		)
		expect_equal(dim(Dinv), c(s$num_treats, s$num_treats))

		# (a) Round-trip: D^{-1} applied to the event-study differences of a
		# random effect vector recovers that effect vector exactly.
		set.seed(40)
		beta <- rnorm(s$num_treats)
		theta <- .es_difference(beta, s$first_inds, s$n_k, G)
		expect_equal(
			as.numeric(Dinv %*% theta),
			beta,
			tolerance = 1e-10,
			info = sprintf("G=%d T=%d", G, T)
		)

		# (b) Gold standard: build the forward operator D independently (the
		# differencing applied to each basis column) and check D^{-1} D = I.
		D <- apply(
			diag(s$num_treats),
			2,
			.es_difference,
			first_inds = s$first_inds,
			n_k = s$n_k,
			G = G
		)
		expect_equal(
			Dinv %*% D,
			diag(s$num_treats),
			tolerance = 1e-10,
			info = sprintf("G=%d T=%d", G, T)
		)
	}
})

test_that("event-study D meets the Lemma's singular-value bounds (#40)", {
	for (cfg in list(c(3, 6), c(4, 6), c(5, 7), c(6, 8))) {
		G <- cfg[1]
		T <- cfg[2]
		s <- .es_setup(G, T)
		Dinv <- fetwfe:::genInvEventStudyFusionTransformMat(
			s$num_treats,
			s$first_inds,
			G
		)
		# D's singular values are the reciprocals of D^{-1}'s.
		sv_inv <- svd(Dinv)$d
		sigma_max_D <- 1 / min(sv_inv)
		sigma_min_D <- 1 / max(sv_inv)
		expect_lte(sigma_max_D, sqrt(6) + 1e-8)
		expect_gte(sigma_min_D, 1 / (T * sqrt(2 * T)) - 1e-8)
	}
})

test_that("event-study fusion reduces to the all-ones spine when G = 1 (#40)", {
	s <- .es_setup(1, 5)
	Dinv <- fetwfe:::genInvEventStudyFusionTransformMat(
		s$num_treats,
		s$first_inds,
		1
	)
	# With a single cohort there is no cross-cohort fusion, so D^{-1} is the
	# n1 x n1 all-ones lower-triangular matrix -- identical to the default
	# two-way fusion inverse.
	default <- fetwfe:::genInvTwoWayFusionTransformMat(
		s$num_treats,
		s$first_inds,
		1
	)
	expect_equal(Dinv, default)
	expect_true(all(Dinv[lower.tri(Dinv, diag = TRUE)] == 1))
	expect_true(all(Dinv[upper.tri(Dinv)] == 0))
})

# ---------------------------------------------------------------------------
# Recovery (#40): when the true treatment effects vary ONLY by event time
# (constant across cohorts at each e), the event-study penalty recovers them
# with smaller error than the default within-cohort fusion. The truth is built
# directly in (g,t) space and `simulateData()` generates y = X %*% beta, so the
# comparison is not tautological.
# ---------------------------------------------------------------------------
test_that("event_study recovers event-time-structured truth better than default (#40)", {
	G <- 3L
	T <- 6L
	d <- 1L
	num_treats <- fetwfe:::getNumTreats(G = G, T = T)
	first_inds <- fetwfe:::getFirstInds(G = G, T = T)
	treat_inds <- fetwfe:::getTreatInds(
		G = G,
		T = T,
		d = d,
		num_treats = num_treats
	)
	n_k <- c(diff(first_inds), num_treats - first_inds[G] + 1L)

	skip_on_cran()
	# Pure event-study truth: tau_{g, g+e} = 1 + 0.5 * e, identical across
	# cohorts at each event time e.
	es_tes <- numeric(num_treats)
	for (g in seq_len(G)) {
		for (e in 0:(n_k[g] - 1L)) {
			es_tes[first_inds[g] + e] <- 1 + 0.5 * e
		}
	}

	# Averaged over seeds to be flake-resistant: the per-seed improvement is
	# real but modest (mse_e / mse_d ~= 0.88 at seed 40), so assert on the mean
	# ratio and that event_study wins on a clear majority of seeds.
	ratios <- vapply(
		c(40L, 41L, 42L, 7L, 99L),
		function(sd) {
			coefs <- genCoefs(
				G = G,
				T = T,
				d = d,
				density = 0.5,
				eff_size = 2,
				seed = sd
			)
			coefs$beta[treat_inds] <- es_tes
			sim <- simulateData(
				coefs,
				N = 300,
				sig_eps_sq = 1,
				sig_eps_c_sq = 0.5
			)
			fit_d <- suppressWarnings(fetwfeWithSimulatedData(
				sim,
				lambda_selection = "bic"
			))
			fit_e <- suppressWarnings(fetwfeWithSimulatedData(
				sim,
				lambda_selection = "bic",
				fusion_structure = "event_study"
			))
			mse_d <- mean((fit_d$beta_hat[treat_inds] - es_tes)^2)
			mse_e <- mean((fit_e$beta_hat[treat_inds] - es_tes)^2)
			mse_e / mse_d
		},
		numeric(1)
	)

	expect_lt(mean(ratios), 1)
	expect_gte(sum(ratios < 1), 4L)
})

# ---------------------------------------------------------------------------
# Coverage (#40): event_study composes with the configurations the threading
# touches (d = 0, se_type = "cluster", add_ridge = TRUE), and an in-repo
# byte-identical guard that fusion_structure = "cohort" is the no-arg path.
# ---------------------------------------------------------------------------
test_that("event_study composes with d=0 / cluster / add_ridge; default is byte-identical (#40)", {
	skip_on_cran()
	sim <- simulateData(
		genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 11),
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	# Byte-identity guard: fusion_structure = "cohort" == omitting it.
	f_omit <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic"
	))
	f_def <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic",
		fusion_structure = "cohort"
	))
	expect_identical(f_omit$beta_hat, f_def$beta_hat)
	expect_identical(f_omit$att_hat, f_def$att_hat)

	# (a) d = 0.
	sim0 <- simulateData(
		genCoefs(G = 3, T = 5, d = 0, density = 0.5, eff_size = 2, seed = 12),
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	f0 <- suppressWarnings(fetwfeWithSimulatedData(
		sim0,
		lambda_selection = "bic",
		fusion_structure = "event_study"
	))
	expect_identical(f0$fusion_structure, "event_study")
	expect_true(is.finite(f0$att_hat) && is.finite(f0$att_se))

	# (c) add_ridge = TRUE: ridge rows are built in the event-study D^{-1} basis.
	fr <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		lambda_selection = "bic",
		fusion_structure = "event_study",
		add_ridge = TRUE
	))
	expect_true(is.finite(fr$att_hat) && is.finite(fr$att_se))

	# (b) cluster SE: exercises the variance-machinery threading.
	fc <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		fusion_structure = "event_study",
		se_type = "cluster"
	))
	expect_true(is.finite(fc$att_hat))
	expect_true(any(is.finite(fc$catt_df$se) & fc$catt_df$se > 0))
})

# ---------------------------------------------------------------------------
# Plumbing sanity (#40): a small event_study fit yields finite, well-formed
# accessor output across catt_df / eventStudy() / simultaneousCIs().
# ---------------------------------------------------------------------------
test_that("event_study fit produces well-formed accessor output (#40)", {
	skip_if_not_installed("mvtnorm")
	sim <- simulateData(
		genCoefs(G = 2L, T = 4L, d = 1L, density = 0.5, eff_size = 2, seed = 7),
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	fit <- suppressWarnings(fetwfeWithSimulatedData(
		sim,
		fusion_structure = "event_study"
	))
	expect_identical(fit$fusion_structure, "event_study")
	expect_true(is.finite(fit$att_hat))
	expect_true(all(is.finite(fit$catt_df$estimate)))
	fin <- is.finite(fit$catt_df$se) & fit$catt_df$se > 0
	expect_true(any(fin))
	es <- suppressMessages(eventStudy(fit))
	expect_true(all(is.finite(es$estimate)))
	sci <- suppressMessages(simultaneousCIs(fit, family = "cohort"))
	expect_true(all(is.finite(sci$ci$estimate)))
})
