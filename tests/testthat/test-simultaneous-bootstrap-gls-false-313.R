# Tests for #313: Omega-free simultaneousCIs() bootstrap bands on gls = FALSE
# high-dimensional (p >= NT) fits.
#
# #307 made debiasedATT()'s overall-ATT SE Omega-free (a gls = FALSE fit -- un-
# whitened, calc_ses = FALSE, no GLS/REML, no sig_eps_sq -- gets a cluster-robust
# SE). The high-dim simultaneousCIs() bootstrap band is built the same way (the
# empirical per-unit influence function, cluster-robust by construction; the
# center re-fits its own q=1 nuisance on the un-whitened design), so it composes
# un-whitened too. Previously the `has_valid_ses` gate (R/simultaneous_cis.R)
# errored on every calc_ses = FALSE fit; now it lets the high-dim fetwfe bootstrap
# path through while STILL erroring for method = "analytic" (needs sig_eps_sq) and
# for fixed-p (#312) / non-fetwfe (#308).

# ---- gls = FALSE high-dim fixture (raw panel, no supplied variances) ----------
# N=20, T=8, d=12 => p = 376 >> NT = 160 (mirrors test-debiased-att-highdim.R).
.make_hd_panel_313 <- function() {
	set.seed(3)
	N <- 20L
	T <- 8L
	d <- 12L
	cohort_of_unit <- c(rep(0L, 8), rep(2L, 4), rep(3L, 4), rep(4L, 4))
	eff <- c(`2` = 0.5, `3` = 1.75, `4` = 3.0)
	covs <- matrix(stats::rnorm(N * d), N, d)
	do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cohort_of_unit[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T,
				treat = as.integer(g > 0 & (1:T) >= g)
			)
			for (j in 1:d) {
				df[[paste0("x", j)]] <- covs[i, j]
			}
			te <- if (g > 0) eff[[as.character(g)]] else 0
			df$y <- 0.3 *
				(1:T) /
				T +
				te * df$treat +
				0.2 * covs[i, 1] +
				stats::rnorm(T, 0, 0.4)
			df
		})
	)
}

hd_nogls <- fetwfe(
	pdata = .make_hd_panel_313(),
	time_var = "year",
	unit_var = "unit",
	treatment = "treat",
	covs = paste0("x", 1:12),
	response = "y",
	q = 0.5,
	verbose = FALSE,
	gls = FALSE
)

test_that("the #313 fixture is high-dim, gls = FALSE (calc_ses = FALSE, sig_eps_sq = NA)", {
	expect_gte(
		ncol(hd_nogls$internal$X_final),
		nrow(hd_nogls$internal$X_final)
	)
	expect_false(isTRUE(hd_nogls$internal$calc_ses))
	expect_true(is.na(hd_nogls$sig_eps_sq))
})

test_that("gls = FALSE high-dim bootstrap bands run for every family (finite, high-dim)", {
	for (fam in c("event_study", "cohort", "all_post_treatment")) {
		bo <- simultaneousCIs(
			hd_nogls,
			family = fam,
			method = "bootstrap",
			B = 300,
			seed = 1
		)
		expect_s3_class(bo, "simultaneous_cis")
		expect_identical(bo$regime, "high-dimensional")
		expect_true(all(is.finite(bo$ci$simultaneous_ci_low)))
		expect_true(all(is.finite(bo$ci$simultaneous_ci_high)))
		expect_true(all(
			c("feasibility", "converged", "lambda_node") %in% names(bo)
		))
		# studentized sup-t critical value in [pointwise, Bonferroni].
		expect_gte(bo$critical_value, bo$pointwise_critical_value - 1e-8)
		expect_lte(bo$critical_value, bo$bonferroni_critical_value + 1e-8)
	}
})

test_that("gls = FALSE high-dim bootstrap band is deterministic at a fixed seed", {
	a <- simultaneousCIs(
		hd_nogls,
		family = "event_study",
		method = "bootstrap",
		B = 300,
		seed = 7
	)
	b <- simultaneousCIs(
		hd_nogls,
		family = "event_study",
		method = "bootstrap",
		B = 300,
		seed = 7
	)
	expect_identical(a$critical_value, b$critical_value)
	expect_identical(a$ci$estimate, b$ci$estimate)
})

test_that("gls = FALSE high-dim band center matches debiasedATT() on the overall ATT (#303 cross-check)", {
	# The overall-ATT custom contrast: cohort g loads its treatment cells with
	# weight cohort_probs[g] / (#cells). The desparsified band centers on the
	# Theorem 6.6 debiased estimate, equal to debiasedATT()$att.
	G <- length(hd_nogls$cohort_probs)
	Tt <- hd_nogls$T
	nt <- length(hd_nogls$treat_inds)
	sizes <- (Tt - 1):(Tt - G)
	cohort_of_treat <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(nt)
	for (g in seq_len(G)) {
		idx <- which(cohort_of_treat == g)
		a_beta[idx] <- hd_nogls$cohort_probs[g] / length(idx)
	}
	sc <- simultaneousCIs(
		hd_nogls,
		family = "custom",
		contrasts = matrix(a_beta, nrow = 1),
		method = "bootstrap",
		B = 100,
		seed = 1
	)
	expect_identical(sc$regime, "high-dimensional")
	expect_equal(sc$ci$estimate, debiasedATT(hd_nogls)$att, tolerance = 1e-9)
})

test_that("method = 'analytic' on a gls = FALSE fit still errors cleanly (not a cryptic NA crash)", {
	# The analytic band needs sig_eps_sq (NA here); the gate must catch it with the
	# calc_ses = FALSE message, NOT crash downstream in Sigma_1.
	for (fam in c("cohort", "event_study")) {
		expect_error(
			simultaneousCIs(hd_nogls, family = fam, method = "analytic"),
			"calc_ses = FALSE"
		)
	}
})

test_that("a fixed-p gls = FALSE fit still errors for the bootstrap band (the #312 scope boundary)", {
	# Small d => p < NT (fixed-p). gls = FALSE => calc_ses = FALSE. The Omega-free
	# band is high-dim-only here; the fixed-p cluster-robust band is #312.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	fp_nogls <- fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
	expect_lt(ncol(fp_nogls$internal$X_final), nrow(fp_nogls$internal$X_final))
	expect_false(isTRUE(fp_nogls$internal$calc_ses))
	expect_error(
		simultaneousCIs(fp_nogls, family = "cohort", method = "bootstrap"),
		"calc_ses = FALSE"
	)
})
