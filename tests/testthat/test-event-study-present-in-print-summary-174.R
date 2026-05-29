library(testthat)
library(fetwfe)

# Issue #174, event-study-present guard.
#
# After #174 fixes `eventStudy()` to handle non-consecutive cohort
# offsets (the `bacondecomp::divorce` reproducer in the issue body),
# the print/summary "Event Study" section must render on BOTH
# consecutive-cohort fixtures (the synthetic genCoefs-based panels)
# and scattered-cohort fixtures (divorce or an equivalent synthetic
# panel where cohort years are non-consecutive). Pre-#174 the divorce
# case produced no Event Study section at all because of the
# silent-swallow tryCatch.
#
# Plus a helper-parity test asserting `getFirstIndsFromOffsets()`
# (the new offset-aware variant added in `R/utility.R`) reduces to
# the existing `getFirstInds(R, T)` on consecutive offsets `2..R+1`.
# This guards byte-identical results on synthetic fixtures that use
# `genCoefs()`-based panels (their `cohort_probs` carry no integer-
# coercible names, so the offset-aware path falls back to
# `getFirstInds`; the parity test guarantees the math is the same
# regardless of which path is taken).

# ----------------------------------------------------------------------
# 1) Event-study-present on the consecutive-cohort synthetic fixture.
#    `genCoefs(R, T)` panels have cohorts at consecutive offsets
#    2, 3, ..., R + 1; this is what the synthetic test suite has been
#    exercising since PR #156.
# ----------------------------------------------------------------------

# The print method renders "Event-Study Average Treatment Effects"
# (hyphenated); the summary method renders "Event Study (preview):"
# (space). Both forms count as "the section rendered"; the regex
# matches either spelling.
.es_header_pattern <- "Event[ -]Study"

test_that("print(res) and summary(res) render an Event Study section on consecutive-cohort synthetic fits", {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	sim <- simulateData(coefs, N = 60, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
	for (res in list(
		fetwfeWithSimulatedData(sim, q = 0.5),
		etwfeWithSimulatedData(sim),
		betwfeWithSimulatedData(sim, q = 0.5)
	)) {
		out_print <- paste(capture.output(print(res)), collapse = "\n")
		expect_match(out_print, .es_header_pattern)

		out_summary <- paste(
			capture.output(print(summary(res))),
			collapse = "\n"
		)
		expect_match(out_summary, .es_header_pattern)
	}
})

# ----------------------------------------------------------------------
# 2) Event-study-present on the scattered-cohort case. Gated on
#    `bacondecomp` (a Suggests:-style dependency under the test gate
#    `skip_if_not_installed`). Pre-#174 this fit would produce no Event
#    Study section because `eventStudy()` raised
#    `first_inds[R] <= n_vars is not TRUE` and the silent-swallow
#    tryCatch hid it.
# ----------------------------------------------------------------------

test_that("print(res) and summary(res) render an Event Study section on divorce (scattered cohorts)", {
	skip_if_not_installed("bacondecomp")
	bacondecomp::divorce
	data(divorce, package = "bacondecomp")
	res <- suppressWarnings(fetwfe(
		pdata = divorce[divorce$sex == 2, ],
		time_var = "year",
		unit_var = "st",
		treatment = "changed",
		covs = c("murderrate", "lnpersinc", "afdcrolls"),
		response = "suiciderate_elast_jag",
		sig_eps_sq = 0.0344,
		sig_eps_c_sq = 0.1507,
		add_ridge = TRUE,
		verbose = FALSE
	))
	# Confirm we're exercising the scattered-cohort path (the bug's
	# reproducer numbers): R = 12, T = 33, non-consecutive cohort offsets.
	expect_identical(as.integer(res$R), 12L)
	expect_identical(as.integer(res$T), 33L)
	offsets <- fetwfe:::.derive_cohort_offsets_from_fit(res)
	expect_false(identical(offsets, seq.int(2L, res$R + 1L)))

	# Sanity: eventStudy() returns a non-empty data frame
	# (pre-#174 it errored before reaching this line).
	es <- eventStudy(res)
	expect_s3_class(es, "eventStudy")
	expect_true(nrow(es) > 0L)

	# Print / summary render the section.
	out_print <- paste(capture.output(print(res)), collapse = "\n")
	expect_match(out_print, .es_header_pattern)
	out_summary <- paste(
		capture.output(print(summary(res))),
		collapse = "\n"
	)
	expect_match(out_summary, .es_header_pattern)
})

# ----------------------------------------------------------------------
# 3) Helper-parity. `getFirstIndsFromOffsets(2:(R+1), T)` must match
#    `getFirstInds(R, T)` byte-identically across a battery of
#    (R, T) pairs covering the existing synthetic fixtures. Without
#    this guard, a future refactor of `getFirstInds()` could drift
#    from the offset-aware variant and break the consecutive-cohort
#    path.
# ----------------------------------------------------------------------

test_that("getFirstIndsFromOffsets reduces to getFirstInds on consecutive offsets", {
	cases <- list(
		c(R = 1, T = 3),
		c(R = 2, T = 4),
		c(R = 2, T = 5),
		c(R = 3, T = 5),
		c(R = 3, T = 6),
		c(R = 5, T = 7),
		c(R = 10, T = 20),
		c(R = 12, T = 33) # divorce-scale, with consecutive offsets
	)
	for (rc in cases) {
		R_v <- as.integer(rc[["R"]])
		T_v <- as.integer(rc[["T"]])
		old <- fetwfe:::getFirstInds(R = R_v, T = T_v)
		new <- fetwfe:::getFirstIndsFromOffsets(
			cohort_offsets_int = seq.int(2L, R_v + 1L),
			T = T_v
		)
		expect_identical(
			as.integer(new),
			as.integer(old),
			info = paste0("R = ", R_v, ", T = ", T_v)
		)
	}
})

# ----------------------------------------------------------------------
# 4) Helper-parity for the scattered-offset case: hand-verified
#    expected first_inds for a small explicit panel.
# ----------------------------------------------------------------------

test_that("getFirstIndsFromOffsets matches a hand-verified scattered-cohort layout", {
	# Three cohorts at offsets (2, 4, 5) in a T = 6 panel.
	#   cohort 1 (offset 2): T - 2 + 1 = 5 effects, indices  1..5
	#   cohort 2 (offset 4): T - 4 + 1 = 3 effects, indices  6..8
	#   cohort 3 (offset 5): T - 5 + 1 = 2 effects, indices  9..10
	# first_inds = (1, 6, 9); num_treats = 10.
	got <- fetwfe:::getFirstIndsFromOffsets(
		cohort_offsets_int = c(2L, 4L, 5L),
		T = 6L
	)
	expect_identical(as.integer(got), c(1L, 6L, 9L))

	# Another: divorce's actual offsets (1969-1985 with gaps) on T = 33.
	divorce_offsets <- c(6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 17L, 21L, 22L)
	got_div <- fetwfe:::getFirstIndsFromOffsets(
		cohort_offsets_int = divorce_offsets,
		T = 33L
	)
	expect_identical(
		as.integer(got_div),
		c(1L, 29L, 56L, 82L, 107L, 131L, 154L, 176L, 197L, 217L, 234L, 247L)
	)
	# Sanity: last first_ind + last cohort effects - 1 == num_treats.
	last_first <- got_div[length(got_div)]
	last_count <- 33L - divorce_offsets[length(divorce_offsets)] + 1L
	expect_identical(
		as.integer(last_first + last_count - 1L),
		258L
	)
})
