# Behavior-visible items from the 2026-07 review sweep (#399).
#
# Item 10: an invalid-argument error on the data-generation entry points must
#          NOT advance the caller's ambient RNG (the seed is applied only after
#          all argument validation succeeds).
# Item 11: processFactors() names factor dummies as `factorVar_levelName`
#          (e.g. `group_B`), not the leaked `model.matrix` term label
#          (`group_pdata[[v]]B`).
# Item 12: the targeted-sparsity degeneracy guards compare against a
#          scale-relative tolerance, so a floating-point-cancelling DGP is
#          rejected while a legitimately small (but nonzero-variance) one is not.

# ---- Item 11: processFactors clean dummy names -------------------------------

test_that("processFactors() names factor dummies as factorVar_levelName (#399)", {
	df <- data.frame(
		y = 1:6,
		group = factor(c("A", "B", "A", "C", "B", "C")),
		x1 = 11:16
	)
	res <- processFactors(df, c("group", "x1"))

	# Baseline level "A" is dropped; "B"/"C" are retained as group_B / group_C.
	# The names must equal paste(factorVar, level, sep = "_") EXACTLY, not the
	# leaked model.matrix term label (which produced e.g. `group_pdata[[v]]B`).
	dummy_covs <- setdiff(res$covs, "x1")
	expect_identical(sort(dummy_covs), c("group_B", "group_C"))
	expect_true("x1" %in% res$covs)
	expect_false("group" %in% res$covs)

	# No leak of the internal `pdata[[v]]` term label into any covariate name.
	expect_false(any(grepl("pdata", res$covs, fixed = TRUE)))
	expect_false(any(grepl("[[", res$covs, fixed = TRUE)))

	# The same clean names appear as columns of the returned data frame.
	expect_true(all(c("group_B", "group_C") %in% colnames(res$pdata)))
	expect_false("group" %in% colnames(res$pdata))
})

# ---- Item 10: invalid args do not advance the ambient RNG --------------------

test_that("genCoefs() invalid-arg error leaves the ambient RNG untouched (#399)", {
	# A non-NULL `seed` makes the pre-#399 bug observable: the old order applied
	# set.seed(seed) BEFORE the scalar validations, so an invalid `T` reseeded
	# (and thus perturbed) the caller's ambient .Random.seed on the error path.
	set.seed(1)
	before <- .Random.seed
	expect_error(
		genCoefs(G = 3, T = -1, d = 2, density = 0.5, eff_size = 1, seed = 123)
	)
	expect_identical(.Random.seed, before)

	# Draw-based cross-check: the failed call must not shift the stream, so the
	# next draw matches a fresh set.seed(1) draw.
	set.seed(1)
	expect_error(
		genCoefs(G = 3, T = -1, d = 2, density = 0.5, eff_size = 1, seed = 123)
	)
	got <- runif(1)
	set.seed(1)
	expect_identical(got, runif(1))
})

test_that("simulateDataCore() invalid-arg error leaves the ambient RNG untouched (#399)", {
	num_treats <- getNumTreats(G = 3, T = 4)
	beta_ok <- rep(0.1, getP(G = 3, T = 4, d = 1, num_treats = num_treats))

	set.seed(2)
	before <- .Random.seed
	expect_error(
		simulateDataCore(
			N = -5,
			T = 4,
			d = 1,
			G = 3,
			sig_eps_sq = 1,
			sig_eps_c_sq = 1,
			beta = beta_ok,
			seed = 123,
			distribution = "gaussian"
		)
	)
	expect_identical(.Random.seed, before)
})

# ---- Item 12: scale-relative targeted-sparsity degeneracy guards -------------

test_that("genCoefs() targeted degeneracy guards use a scale-relative tolerance (#399)", {
	G <- 3
	T <- 5
	d <- 0
	mk <- function(tbl) {
		genCoefs(
			G = G,
			T = T,
			d = d,
			density = 0.5,
			eff_size = 1,
			treat_base_levels = tbl,
			seed = 1
		)
	}

	# Truly degenerate cases are still rejected.
	expect_error(mk(rep(0, G)), "all-zero cohort effects") # att == 0 (guard 1)
	expect_error(mk(c(1, 0, 0)), "homogeneous cohort effects") # V2 == 0 (guard 2)

	# Near-homogeneous below tolerance: the cohort effects differ only at ~1e-9.
	# The OLD exact `length(unique(catt_built)) == 1L` guard would have PASSED
	# this (three distinct doubles); the scale-relative tolerance rejects it as
	# numerically homogeneous (V2 ~ 0).
	expect_error(mk(c(1, 1e-9, 0)))

	# Legitimately small-but-nonzero heterogeneity must NOT be false-rejected.
	expect_s3_class(mk(c(0, 0.01, 0.02)), "FETWFE_coefs")
	expect_s3_class(mk(c(1, 1e-7, 0)), "FETWFE_coefs")
})
