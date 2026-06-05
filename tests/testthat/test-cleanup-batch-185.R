# Regression tests for the actionable items of the #185 cleanup batch.
#
# Covered here:
#   SB2 -- the estimator input validator rejects sig_eps_sq = 0 (a zero noise
#          variance makes the GLS transform 1 / sqrt(sig_eps_sq) non-finite).
#          Previously only `< 0` was rejected.
#   SB6 -- getNumTreats() / getFirstInds() / getTreatInds() return integer
#          vectors, honoring their @return contracts (previously double via the
#          `/ 2` divisions and `seq()`).
#
# Other #185 items, covered elsewhere or not actionable:
#   SB5 (ridge no-op) -- guarded by an invariant stopifnot in
#       R/gls_machinery.R; exercised by the existing add_ridge tests.
#   A2  (validation-timing asymmetry) -- documentation-only.
#   G4  (synthetic scattered-cohort fixture) -- in
#       test-event-study-present-in-print-summary-174.R.
#   IC1 (CV ignores lambda.* args) -- already warns; tested in
#       test-lambda-selection-164.R.
#   SB4 (1e-6 cohort-probs tolerance) -- won't-fix: the 1e-6 band is the
#       deliberate, tested choice of #56 (test-soft-bugs-56.R), and firing at
#       <= 1e-6 never-treated mass is protective, not spurious.

# ------------------------------------------------------------------------------
# SB6: the base-treatment-effect index helpers return integer.
# ------------------------------------------------------------------------------

test_that("getNumTreats / getFirstInds / getTreatInds return integer vectors (#185 SB6)", {
	nt <- fetwfe:::getNumTreats(G = 3L, T = 6L)
	expect_type(nt, "integer")
	expect_identical(nt, 12L) # 6 * 3 - 3 * 4 / 2 = 12

	fi <- fetwfe:::getFirstInds(G = 3L, T = 6L)
	expect_type(fi, "integer")
	expect_identical(fi, c(1L, 6L, 10L)) # hand-checked for G = 3, T = 6

	ti <- fetwfe:::getTreatInds(G = 3L, T = 6L, d = 2L, num_treats = nt)
	expect_type(ti, "integer")
	expect_length(ti, 12L)
})

# ------------------------------------------------------------------------------
# SB2: the estimator input validator rejects sig_eps_sq = 0 (not just < 0).
# ------------------------------------------------------------------------------

.fm185_panel <- function() {
	set.seed(185)
	mk <- function(ids, adopt) {
		do.call(
			rbind,
			lapply(ids, function(u) {
				tr <- if (is.na(adopt)) {
					integer(4L)
				} else {
					as.integer(1:4 >= adopt)
				}
				data.frame(
					unit = as.character(u),
					year = 1:4,
					x1 = rnorm(1),
					treated = tr,
					y = rnorm(4),
					stringsAsFactors = FALSE
				)
			})
		)
	}
	rbind(mk(1:8, 3L), mk(9:16, 4L), mk(17:24, NA))
}

test_that("input validator rejects sig_eps_sq = 0 with a positive-only message (#185 SB2)", {
	pd <- .fm185_panel()
	v0 <- fetwfe:::.collect_etwfe_input_violations(
		pdata = pd,
		time_var = "year",
		unit_var = "unit",
		treatment = "treated",
		response = "y",
		covs = "x1",
		sig_eps_sq = 0
	)
	expect_true(any(grepl("sig_eps_sq must be positive", v0)))

	# A valid positive variance produces no sig_eps_sq violation.
	v1 <- fetwfe:::.collect_etwfe_input_violations(
		pdata = pd,
		time_var = "year",
		unit_var = "unit",
		treatment = "treated",
		response = "y",
		covs = "x1",
		sig_eps_sq = 1
	)
	expect_false(any(grepl("sig_eps_sq", v1)))
})
