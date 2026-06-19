# ------------------------------------------------------------------------------
# #308: betwfe (non-fetwfe) high-dimensional (p >= NT) bootstrap band coverage.
#
# A non-fetwfe high-dim fit has no desparsified (Theorem 6.6) band, so
# simultaneousCIs(method = "bootstrap") falls back to the post-selection fixed-p
# selected-support band (preserving #305: it returns a band, not an error). This
# test pins the empirical finding that motivated the #308 warning -- that
# fallback band substantially UNDER-COVERS the true cohort effects in the
# p >= NT regime, *even in a favorable sparse regime* where the selected support
# is low-dimensional. The fixture's truth is (4.5, 0, 0): a single non-zero
# cohort effect. The bridge shrinks the band center on that cohort toward zero
# (biasing it away from 4.5) and the post-selection SE understates sampling
# variability, so simultaneous coverage collapses far below the nominal 0.95.
#
# Measured (#308 fixed-p path): simultaneous coverage ~ 0.09. The assertion uses
# a loose < 0.6 threshold -- hugely Monte-Carlo-robust at this rep count (true
# coverage ~0.09 would need >= ~20/33 covered reps to breach it, probability
# ~1e-17). The data, bridge fit, and bootstrap are all seeded, so coverage is
# deterministic across runs (not flaky). Mutation-checkable: removing the
# warning() in the non-fetwfe high-dim bootstrap branch of .simultaneous_cis_impl()
# fails the expect_warning in part (a); the coverage assertion in part (b) is
# independent of the warning.
# ------------------------------------------------------------------------------

test_that("non-fetwfe (betwfe) high-dim bootstrap band under-covers + warns (#308)", {
	skip_on_cran()

	G <- 3L
	T_ <- 5L
	# density 0.08 + large eff_size 6 => sparse truth (one non-zero cohort): the
	# regime most favorable to the fixed-p fallback, yet it still under-covers.
	coefs <- genCoefs(
		G = G,
		T = T_,
		d = 20,
		density = 0.08,
		eff_size = 6,
		seed = 11
	)
	truth <- getTes(coefs)$actual_cohort_tes # length-G cohort ATTs; ~ (4.5, 0, 0)

	fit_one <- function(rep_seed) {
		dat <- simulateData(
			coefs,
			N = 60,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = rep_seed
		)
		dat$indep_counts <- NA
		betwfe(
			pdata = dat$pdata,
			time_var = dat$time_var,
			unit_var = dat$unit_var,
			treatment = dat$treatment,
			response = dat$response,
			covs = dat$covs,
			q = 0.5,
			verbose = FALSE,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5
		)
	}

	# (a) Deterministic warning + the p >= NT and fixed-p-fallback invariants on a
	# single fit (this is the warning's mutation anchor).
	bfit1 <- fit_one(1001L)
	expect_gte(ncol(bfit1$internal$X_final), nrow(bfit1$internal$X_final)) # p >= NT
	expect_warning(
		bo1 <- simultaneousCIs(
			bfit1,
			family = "cohort",
			method = "bootstrap",
			B = 200,
			seed = 1L
		),
		"under-?cover|unreliable"
	)
	expect_identical(bo1$regime, "fixed-p") # NOT the fetwfe-only desparsified path
	expect_identical(nrow(bo1$ci), G) # G cohort rows, in cohort order

	# (b) Coverage Monte Carlo over the #308 fixed-p fallback band. A minority of
	# reps' bridge zeroes the ENTIRE selected support; those hit the separate #304
	# degenerate early-exit (a zero-width band, `regime` NULL) rather than the
	# fixed-p fallback band, so to measure exactly the band this issue is about we
	# restrict the coverage estimate to the genuine #308-path reps (regime
	# "fixed-p"). (Coverage is far below nominal on the blended set too, but the
	# attribution -- bridge-shrunk center + understated post-selection SE -- is the
	# #308 path's.)
	nreps <- 40L
	covered <- logical(nreps)
	on_308_path <- logical(nreps)
	for (r in seq_len(nreps)) {
		bfit <- fit_one(1000L + r)
		bo <- suppressWarnings(simultaneousCIs(
			bfit,
			family = "cohort",
			method = "bootstrap",
			B = 200,
			seed = r
		))
		on_308_path[r] <- identical(bo$regime, "fixed-p")
		covered[r] <- all(
			truth >= bo$ci$simultaneous_ci_low &
				truth <= bo$ci$simultaneous_ci_high
		)
	}
	# the MC must actually exercise the fixed-p fallback path in most reps
	expect_gt(sum(on_308_path), 20L)
	coverage <- mean(covered[on_308_path]) # nominal 0.95; ~0.09 on the #308 path
	expect_lt(coverage, 0.6)
})
