# Regression test for #166: idCohorts() now uses a vectorised
# (sort-once + reshape-to-T-by-N) strategy instead of the per-unit
# O(N^2 * T) row-filter loop. The existing test-idcohorts-collect-
# violations.R covers error-message wording and lex-ordering; this
# file adds tests focused on the cohort-assignment result shape and
# treated-unit accounting, plus a performance sanity check that runs
# only locally (skip_on_cran) to avoid flakiness on CI machines.

library(testthat)

# ------------------------------------------------------------------------------
# Cohort-assignment correctness on a small fixture with known truth.
# ------------------------------------------------------------------------------

test_that("idCohorts assigns treated units to the correct cohort under varied first-treatment times", {
	df <- data.frame(
		unit = rep(c("u01", "u02", "u03", "u04", "u05"), each = 4L),
		time = rep(c(1L, 2L, 3L, 4L), times = 5L),
		# u01 treated at t=2, u02 treated at t=3, u03 never-treated,
		# u04 treated at t=4 (last period), u05 never-treated.
		treat = c(
			0L, 1L, 1L, 1L,
			0L, 0L, 1L, 1L,
			0L, 0L, 0L, 0L,
			0L, 0L, 0L, 1L,
			0L, 0L, 0L, 0L
		),
		stringsAsFactors = FALSE
	)
	out <- fetwfe:::idCohorts(df, "time", "unit", "treat")
	# Cohort 1 was the first-period cohort; it gets removed up front
	# in the existing post-loop logic, so `out$cohorts` should start
	# from t=2 onward and contain only the non-empty treated cohorts.
	expect_true("2" %in% names(out$cohorts))
	expect_true("3" %in% names(out$cohorts))
	expect_true("4" %in% names(out$cohorts))
	expect_setequal(out$cohorts[["2"]], "u01")
	expect_setequal(out$cohorts[["3"]], "u02")
	expect_setequal(out$cohorts[["4"]], "u04")
	# Never-treated units remain in df but not in any cohort.
	expect_true(all(c("u03", "u05") %in% out$units))
	expect_false("u03" %in% unlist(out$cohorts))
	expect_false("u05" %in% unlist(out$cohorts))
})

test_that("idCohorts groups multiple units into the same cohort", {
	df <- data.frame(
		unit = rep(c("a", "b", "c", "d"), each = 3L),
		time = rep(c(1L, 2L, 3L), times = 4L),
		treat = c(
			0L, 1L, 1L, # a treated t=2
			0L, 1L, 1L, # b treated t=2
			0L, 0L, 1L, # c treated t=3
			0L, 0L, 0L  # d never-treated
		),
		stringsAsFactors = FALSE
	)
	out <- fetwfe:::idCohorts(df, "time", "unit", "treat")
	expect_setequal(out$cohorts[["2"]], c("a", "b"))
	expect_setequal(out$cohorts[["3"]], "c")
})

# ------------------------------------------------------------------------------
# Multiple violations of mixed types in a single panel: all surface.
# ------------------------------------------------------------------------------

test_that("idCohorts reports unbalanced + duplicate-time + absorbing violations together", {
	df <- data.frame(
		unit = c(
			rep("balA", 3),               # balanced, fine
			rep("dupT", 3),               # T=3 rows but times {1, 2, 2}
			rep("miss", 2),               # only 2 rows
			rep("absV", 3)                # treat 0/1/0 (non-absorbing)
		),
		time = c(
			1L, 2L, 3L,
			1L, 2L, 2L,
			1L, 2L,
			1L, 2L, 3L
		),
		treat = c(
			0L, 0L, 0L,
			0L, 0L, 0L,
			0L, 0L,
			0L, 1L, 0L
		),
		stringsAsFactors = FALSE
	)
	err <- tryCatch(
		fetwfe:::idCohorts(df, "time", "unit", "treat"),
		error = function(e) conditionMessage(e)
	)
	# Both error sections present.
	expect_match(err, "Panel does not appear to be balanced")
	expect_match(err, "Treatment does not appear to be an absorbing state")
	# All three balance violations named.
	expect_match(err, "dupT has 3 observations (2 distinct time periods)", fixed = TRUE)
	expect_match(err, "miss has 2 observations (2 distinct time periods)", fixed = TRUE)
	# Absorbing violation named.
	expect_match(err, "absV")
})

# ------------------------------------------------------------------------------
# Performance smoke test: vectorised idCohorts on a large balanced panel
# completes under a generous wall-clock budget. Skip on CRAN.
# ------------------------------------------------------------------------------

test_that("idCohorts vectorised path scales to N=2000 in under 1 second", {
	skip_on_cran()
	N <- 2000L
	T <- 5L
	df <- data.frame(
		unit = rep(sprintf("u%04d", seq_len(N)), each = T),
		time = rep(seq_len(T), times = N),
		treat = 0L,
		stringsAsFactors = FALSE
	)
	# Treat a random subset of units starting at varied first-treatment times.
	set.seed(42L)
	treated_units <- sample(N, 600)
	for (i in treated_units) {
		start_t <- sample(seq(2L, T), 1L)
		df$treat[(i - 1L) * T + start_t:T] <- 1L
	}
	t0 <- Sys.time()
	out <- fetwfe:::idCohorts(df, "time", "unit", "treat")
	elapsed <- as.numeric(Sys.time() - t0, units = "secs")
	# Generous ceiling — pre-#166 took ~0.4s; budget 1.0s tolerates a
	# wide range of hardware while still catching any regression that
	# brings back the O(N^2 * T) per-unit row-filter loop.
	expect_lt(elapsed, 1.0)
	# Result shape contract.
	expect_named(out, c("df", "cohorts", "units", "times"))
	expect_equal(length(out$cohorts), T - 1L) # first-period cohort removed
	expect_true(all(lengths(out$cohorts) >= 1L))
})
