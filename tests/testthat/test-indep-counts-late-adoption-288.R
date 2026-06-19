# Regression tests for #288: cohortTimeATTs() / eventStudy() / simultaneousCIs()
# errored ("first_inds[G] <= n_vars is not TRUE", or "subscript out of bounds"
# for simultaneousCIs) on a fit given `indep_counts` when the treated cohorts
# adopt LATER than the earliest panel periods.
#
# Root cause: the indep_counts path stored cohort_probs/cohort_probs_overall
# UNNAMED, so .derive_cohort_offsets_from_fit() could not recover the adoption
# years and the accessors fell back to earliest-timing getFirstInds(), producing
# a first_inds inconsistent with the actual num_treats. Two layered fixes:
#   (1) .select_att_branch() now carries the cohort-year names onto the indep
#       vectors (forward fix for new fits) -- bitten by the name/first_inds
#       parity assertions below.
#   (2) .derive_cohort_offsets_from_fit() falls back to catt_df$cohort when
#       cohort_probs names are absent (rescues legacy/persisted pre-fix objects)
#       -- bitten by the "legacy object" test below.
#
# The fixture is pure base R (no bacondecomp / dplyr) and uses
# `indep_counts == in-sample counts`, so the cohort probabilities -- and hence
# all POINT estimates -- are identical with and without indep_counts (only the
# SE *formula* differs by design, so SEs are not asserted equal). Synthetic
# genCoefs fixtures do NOT expose the bug: their cohorts adopt consecutively
# from period 2, where the earliest-timing fallback is already correct.

# A balanced panel whose cohorts adopt at t = 6, 7, 8, 9 in a 1..11 panel
# (earliest offset 6 > 2 = the trigger), plus never-treated units. d = 0.
.make_late_adoption_panel <- function() {
	set.seed(11)
	cohort_of_unit <- c(
		rep(0L, 22),
		rep(6L, 8),
		rep(7L, 7),
		rep(8L, 6),
		rep(9L, 9)
	) # N = 52
	N <- length(cohort_of_unit)
	T <- 11L
	pdata <- do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cohort_of_unit[i]
			treated <- as.integer(g > 0 & (1:T) >= g)
			data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T,
				treat = treated,
				y = 0.5 *
					i /
					N +
					0.3 * (1:T) / T +
					0.8 * treated +
					rnorm(T, 0, 0.5),
				stringsAsFactors = FALSE
			)
		})
	)
	# indep_counts in cohort order c(never, 6, 7, 8, 9); == in-sample counts.
	ic <- as.integer(table(factor(cohort_of_unit, levels = c(0, 6, 7, 8, 9))))
	list(pdata = pdata, indep_counts = ic)
}

.fit_estimator <- function(estimator, pdata, indep_counts = NULL) {
	fn <- switch(estimator, fetwfe = fetwfe, etwfe = etwfe, betwfe = betwfe)
	args <- list(
		pdata = pdata,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		verbose = FALSE
	)
	if (!is.null(indep_counts)) {
		args$indep_counts <- indep_counts
	}
	do.call(fn, args)
}

for (est in c("fetwfe", "etwfe", "betwfe")) {
	test_that(
		paste0(
			"indep_counts + late adoption: accessors succeed and labels match (",
			est,
			")"
		),
		{
			panel <- .make_late_adoption_panel()
			fit_no <- .fit_estimator(est, panel$pdata)
			fit_ic <- .fit_estimator(est, panel$pdata, panel$indep_counts)

			# (a) Forward-fix smoke test: the three accessors that crashed must now
			# succeed on the indep_counts fit. These error ("first_inds[G] <=
			# n_vars" / "subscript out of bounds") only if BOTH fixes are reverted
			# --- either the names-copy (fix #1) OR the catt_df fallback (fix #2)
			# rescues them. Fix #1 is isolated by the name/first_inds parity
			# assertions in (b); fix #2 is isolated by the legacy test below.
			expect_no_error(eventStudy(fit_ic))
			expect_no_error(cohortTimeATTs(fit_ic))
			expect_no_error(simultaneousCIs(fit_ic))

			# (b) The indep path now carries the same cohort-year names + resolves
			# the same first_inds as the no-indep path (bites the source fix).
			expect_identical(
				names(fit_ic$cohort_probs),
				names(fit_no$cohort_probs)
			)
			expect_identical(names(fit_ic$cohort_probs), c("6", "7", "8", "9"))
			fi_no <- .resolve_cohort_offsets_and_first_inds(
				fit_no,
				fit_no$G,
				fit_no$T
			)$first_inds
			fi_ic <- .resolve_cohort_offsets_and_first_inds(
				fit_ic,
				fit_ic$G,
				fit_ic$T
			)$first_inds
			expect_identical(as.integer(fi_ic), as.integer(fi_no))
			expect_identical(as.integer(fi_ic), c(1L, 7L, 12L, 16L))

			# (c) indep_counts == in-sample counts, so all POINT estimates are
			# identical across paths (indep changes only the SE form, not asserted).
			expect_equal(fit_ic$att_hat, fit_no$att_hat)
			expect_equal(
				eventStudy(fit_ic)$estimate,
				eventStudy(fit_no)$estimate
			)
			expect_equal(
				cohortTimeATTs(fit_ic)$estimate,
				cohortTimeATTs(fit_no)$estimate
			)
		}
	)
}

test_that("legacy indep fit with unnamed cohort_probs is rescued via catt_df", {
	# Simulate a pre-#288 saved object: the indep_counts fix did not yet attach
	# names, so cohort_probs / cohort_probs_overall are unnamed. The catt_df
	# fallback in .derive_cohort_offsets_from_fit() must still recover the
	# correct adoption timing so the accessors do not crash. (Reverting that
	# fallback makes this test error with "first_inds[G] <= n_vars".)
	panel <- .make_late_adoption_panel()
	legacy <- .fit_estimator("fetwfe", panel$pdata, panel$indep_counts)
	names(legacy$cohort_probs) <- NULL
	names(legacy$cohort_probs_overall) <- NULL

	# catt_df$cohort still carries the adoption years -> correct offsets.
	off <- .derive_cohort_offsets_from_fit(legacy)
	expect_identical(as.integer(off), c(6L, 7L, 8L, 9L))
	expect_identical(
		as.integer(
			.resolve_cohort_offsets_and_first_inds(
				legacy,
				legacy$G,
				legacy$T
			)$first_inds
		),
		c(1L, 7L, 12L, 16L)
	)
	expect_no_error(eventStudy(legacy))
	expect_no_error(cohortTimeATTs(legacy))
	expect_no_error(simultaneousCIs(legacy))
})

test_that("indep_counts re-weights eventStudy / att_hat but not cohortTimeATTs cells", {
	# Executable guard for the invariance nuance in NEWS: with indep_counts that
	# DIFFER from the in-sample cohort counts, the overall att_hat and the
	# eventStudy() event-time estimates are re-weighted by the independent-sample
	# cohort probabilities, while the per-(cohort, time) cohortTimeATTs() cells are
	# invariant (each cell is a single coefficient with no cross-cohort weighting).
	panel <- .make_late_adoption_panel()
	fit_no <- .fit_estimator("fetwfe", panel$pdata)
	# Skewed indep_counts: same total (N) and same never-treated count, but the
	# treated-cohort weights differ from the in-sample c(8, 7, 6, 9).
	ic_skew <- c(22L, 12L, 3L, 6L, 9L)
	stopifnot(sum(ic_skew) == sum(panel$indep_counts))
	fit_sk <- .fit_estimator("fetwfe", panel$pdata, ic_skew)

	# The #288 fix still applies (accessors succeed), and the per-cell estimates
	# are invariant to the cohort weighting.
	expect_no_error(eventStudy(fit_sk))
	expect_equal(
		cohortTimeATTs(fit_sk)$estimate,
		cohortTimeATTs(fit_no)$estimate
	)
	# ...but the overall ATT and the aggregated event-time estimates re-weight.
	expect_false(isTRUE(all.equal(fit_sk$att_hat, fit_no$att_hat)))
	expect_false(isTRUE(all.equal(
		eventStudy(fit_sk)$estimate,
		eventStudy(fit_no)$estimate
	)))
})
