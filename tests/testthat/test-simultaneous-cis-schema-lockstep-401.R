library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #401 batch 5 -- forward lockstep guards for the Group B consolidations that
# single-source drift-prone STRUCTURE (as opposed to behavior, which the
# pre/post identical() battery + the existing eventStudy/simultaneousCIs suite
# already cover). These assert the invariants the consolidations exist to
# protect, so a future edit that re-introduces the duplication -- or drifts one
# copy -- breaks loudly here.
#   - item 3: the four simultaneousCIs.<class> S3 methods are ONE implementation.
#   - item 6: the analytic and bootstrap `simultaneous_cis` objects (assembled by
#             `.assemble_simultaneous_cis_result()`) share the same core schema.
# ------------------------------------------------------------------------------

test_that("simultaneousCIs.<class> methods are aliases of one implementation (#401 item 3)", {
	# The #351 (riesz_tol), #402, and #400 edits each had to touch all four
	# method signatures in lockstep. Aliasing three to `.fetwfe` makes that
	# impossible to drift; if a future change re-expands one method into its own
	# body, these expectations fail.
	fe <- getS3method("simultaneousCIs", "fetwfe")
	expect_identical(getS3method("simultaneousCIs", "etwfe"), fe)
	expect_identical(getS3method("simultaneousCIs", "betwfe"), fe)
	expect_identical(getS3method("simultaneousCIs", "twfeCovs"), fe)
})

test_that("analytic and bootstrap simultaneous_cis share the same core schema (#401 item 6)", {
	# print()/tidy()/plot() for `simultaneous_cis` read these fields positionally
	# and by name; the analytic path, the analytic degenerate early-return, and
	# the bootstrap path must stay aligned. `.assemble_simultaneous_cis_result()`
	# single-sources them -- this pins the contract so a field added to one path
	# but not another fails here.
	core <- c(
		"ci",
		"adjusted_p_values",
		"critical_value",
		"pointwise_critical_value",
		"bonferroni_critical_value",
		"family",
		"alpha",
		"K"
	)
	cf <- genCoefs(G = 4, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 270)
	sim <- simulateData(
		cf,
		N = 160,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 270
	)
	fit <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		verbose = FALSE
	)
	an <- simultaneousCIs(fit, family = "cohort", method = "analytic")
	bo <- simultaneousCIs(
		fit,
		family = "cohort",
		method = "bootstrap",
		B = 100L,
		seed = 1L
	)
	expect_s3_class(an, "simultaneous_cis")
	expect_s3_class(bo, "simultaneous_cis")
	# Analytic carries EXACTLY the core fields, in order.
	expect_identical(names(an), core)
	# Bootstrap carries the SAME core fields first, in the same order, then its
	# method/B/seed/multiplier/regime extras.
	expect_identical(names(bo)[seq_along(core)], core)
	expect_true(all(
		c("method", "B", "seed", "multiplier", "regime") %in% names(bo)
	))
	# The per-effect `ci` data frame column schema is shared (what tidy/plot read).
	expect_identical(names(an$ci), names(bo$ci))
})
