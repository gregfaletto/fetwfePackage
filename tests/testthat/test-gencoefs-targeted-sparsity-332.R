# Tests for the genCoefs() / genCoefsCore() targeted-sparsity mode (#332):
# n_signal_cohorts (count) and treat_base_levels (explicit per-cohort fused base
# levels). The mode places a small, heterogeneous treatment signal so the DGP is
# simultaneously sparse (||theta||_0 = signal count) and non-degenerate
# (att != 0, cohort-weight variance V2 > 0) at p >> num_treats -- which the
# uniform-`density` path cannot do in the high-dim regime. Both new args NULL is
# the unchanged uniform path, and must be BYTE-IDENTICAL.
#
# V2 has no field on the FETWFE_tes object, so "V2 > 0" is operationalized as
# sum(cohort_weights * (actual_cohort_tes - att_true)^2) > 0 (the same sign as
# the fit-side .plugin_v2): positive iff the per-cohort effects are heterogeneous.

# Golden vectors captured from the pre-#332 code (origin/main) and verified
# byte-identical there; they pin the default uniform path against drift.
GOLD0_BETA <- c(0, 0, 0, 2, 2, 2, 0, 0, 2, 2, 2, 0, 2, 4, 2, 2)
GOLD0_THETA <- c(0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 2, 2, 2, 0)
GOLD2_BETA <- c(
	0,
	0,
	0,
	-2,
	-2,
	-2,
	0,
	0,
	2,
	0,
	0,
	0,
	2,
	-2,
	-2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	-2,
	0,
	0,
	0,
	-2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	-2,
	-2,
	0,
	0,
	-2,
	-2,
	0,
	0,
	-2,
	-2,
	2,
	0,
	-2,
	-2,
	0,
	0,
	-2,
	-2,
	0,
	0
)
GOLD2_THETA <- c(
	0,
	0,
	0,
	0,
	0,
	-2,
	0,
	0,
	2,
	0,
	0,
	0,
	2,
	-2,
	-2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	-2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	2,
	0,
	0,
	0,
	2,
	0,
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
	0,
	2,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
)

v2_truth <- function(te) {
	sum(te$cohort_weights * (te$actual_cohort_tes - te$att_true)^2)
}

test_that("default (both new args NULL) is byte-identical to the uniform path (#332)", {
	# Golden pin: the uniform-`density` output is unchanged by the new arguments.
	set.seed(7)
	g0 <- genCoefs(G = 3, T = 5, d = 0, density = 0.3, eff_size = 2, seed = 11L)
	expect_identical(g0$beta, GOLD0_BETA)
	expect_identical(g0$theta, GOLD0_THETA)

	set.seed(7)
	g2 <- genCoefs(G = 3, T = 5, d = 4, density = 0.3, eff_size = 2, seed = 11L)
	expect_identical(g2$beta, GOLD2_BETA)
	expect_identical(g2$theta, GOLD2_THETA)

	# Explicitly passing NULL == omitting the args (Convention B: re-seed before
	# each compared call since data-gen advances the ambient stream).
	set.seed(7)
	a <- genCoefs(G = 3, T = 5, d = 4, density = 0.3, eff_size = 2, seed = 11L)
	set.seed(7)
	b <- genCoefs(
		G = 3,
		T = 5,
		d = 4,
		density = 0.3,
		eff_size = 2,
		n_signal_cohorts = NULL,
		treat_base_levels = NULL,
		seed = 11L
	)
	expect_identical(a$beta, b$beta)
	expect_identical(a$theta, b$theta)
})

test_that("count mode yields exact sparsity + att != 0 + V2 > 0 at p >> num_treats (#332)", {
	# High-dim fixture: G=3,T=5,d=40 -> num_treats=9, p=696, NT=200 (p > NT and
	# p >> num_treats), sub-second.
	num_treats <- getNumTreats(G = 3, T = 5)
	p <- getP(G = 3, T = 5, d = 40, num_treats = num_treats)
	expect_gt(p, 3 * 5 * 40 / 3) # sanity: p >> num_treats
	expect_gt(p, 40 * 5) # p > N*T at N = 40

	for (fs in c("cohort", "event_study")) {
		for (k in 1:2) {
			coefs <- genCoefs(
				G = 3,
				T = 5,
				d = 40,
				density = 0.2,
				eff_size = 2,
				fusion_structure = fs,
				n_signal_cohorts = k,
				seed = 1L
			)
			info <- paste0(fs, " k=", k)
			# exact sparsity s = k
			expect_identical(sum(coefs$theta != 0), k, info = info)
			te <- getTes(coefs)
			expect_true(te$att_true != 0, info = info) # non-degenerate ATT
			expect_gt(v2_truth(te), 0) # heterogeneity == V2 > 0
			expect_gt(length(unique(round(te$actual_cohort_tes, 8))), 1L)
			# cohort 1 stays at the baseline (the count mode skips first_inds[1])
			expect_equal(te$actual_cohort_tes[1], 0, info = info)
		}
	}
})

test_that("explicit treat_base_levels places exact bases and stays heterogeneous (#332)", {
	# c(0, m, 0): only cohort 2's base is set -> s = 1, cohorts {2,3} lifted.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 40,
		density = 0.2,
		eff_size = 2,
		treat_base_levels = c(0, 1.5, 0),
		seed = 1L
	)
	expect_identical(sum(coefs$theta != 0), 1L)
	te <- getTes(coefs)
	expect_true(te$att_true != 0)
	expect_gt(v2_truth(te), 0)
	# under cohort fusion the bases map cumulatively: (0, 1.5, 0) -> (0, 1.5, 1.5)
	expect_equal(te$actual_cohort_tes, c(0, 1.5, 1.5))
})

test_that("explicit treat_base_levels under event_study fusion is non-degenerate (#343)", {
	# Same explicit base vector as the cohort-fusion test above, but mapped
	# through the event-study inverse fusion transform. Pins the cumulative ->
	# fused mapping for the explicit-vector x event_study combination (previously
	# only the count mode exercised event_study fusion; the explicit-vector mode
	# tested only the default "cohort" fusion).
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 40,
		density = 0.2,
		eff_size = 2,
		treat_base_levels = c(0, 1.5, 0),
		fusion_structure = "event_study",
		seed = 1L
	)
	expect_identical(sum(coefs$theta != 0), 1L)
	te <- getTes(coefs)
	expect_true(te$att_true != 0)
	expect_gt(v2_truth(te), 0)
	# event_study fusion maps the bases differently from cohort fusion (which
	# gives c(0, 1.5, 1.5) above).
	expect_equal(te$actual_cohort_tes, c(0, 0.5, 0.75))
})

test_that("simulateData() / getTes() consume a targeted coefs object unchanged (#332)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 40,
		density = 0.2,
		eff_size = 2,
		n_signal_cohorts = 2,
		seed = 1L
	)
	sim <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1,
		seed = 7L
	)
	# simulateData reads only beta + metadata; the truth round-trips.
	expect_identical(sim$coefs, coefs$beta)
})

test_that("targeted-sparsity validation + non-degeneracy guards fire (#332)", {
	mk <- function(...) {
		genCoefs(
			G = 3,
			T = 5,
			d = 40,
			density = 0.2,
			eff_size = 2,
			seed = 1L,
			...
		)
	}
	# mutual exclusion
	expect_error(
		mk(n_signal_cohorts = 1, treat_base_levels = c(0, 1, 0)),
		"at most one of"
	)
	# count out of range / non-integer
	expect_error(mk(n_signal_cohorts = 3), "<= G - 1")
	expect_error(mk(n_signal_cohorts = 1.5), "positive ")
	# explicit vector wrong length
	expect_error(mk(treat_base_levels = c(1, 2)), "length G")
	# G < 2
	expect_error(
		genCoefs(
			G = 1,
			T = 5,
			d = 40,
			density = 0.2,
			eff_size = 2,
			n_signal_cohorts = 1,
			seed = 1L
		),
		"requires G >= 2"
	)
	# build-time att == 0 guard (eff_size 0 zeroes the staircase)
	expect_error(
		genCoefs(
			G = 3,
			T = 5,
			d = 40,
			density = 0.2,
			eff_size = 0,
			n_signal_cohorts = 2,
			seed = 1L
		),
		"all-zero cohort effects"
	)
	# build-time V2 == 0 guard: c(m, 0, 0) is the all-ones baseline column ->
	# homogeneous (m, m, m).
	expect_error(
		mk(treat_base_levels = c(2, 0, 0)),
		"homogeneous cohort effects"
	)
})
