# Tests for the idCohorts() multi-unit-violation collection (v1.9.6, #64).
#
# Before this PR, idCohorts() stop()-ed on the FIRST offending unit for
# both the balance check (unit doesn't appear in T periods) and the
# absorbing-state check (treatment, once 1, must stay 1). A user with
# multiple malformed units had to fix one, re-run, hit the next stop(),
# fix again, re-run — one round-trip per bad unit.
#
# After this PR, idCohorts() collects violations across the whole loop
# into two buckets (balance + absorbing-state) and stop()-s once at the
# end with a grouped, lex-sorted listing, truncated to 20 entries per
# bucket. The error message wording starts with the canonical prefix
# ("Panel does not appear to be balanced" / "Treatment does not appear
# to be an absorbing state") so downstream grep-based consumers still
# match; the post-prefix content is the new listing.

# ------------------------------------------------------------------------------
# Test helpers
# ------------------------------------------------------------------------------

# Build a balanced T=3 panel with `n_good` well-behaved units, then mutate
# selected units to introduce balance / absorbing-state violations. Returns
# a data.frame with columns (unit, time, treat).
.make_panel <- function(n_good = 3L) {
	good_units <- sprintf("good%02d", seq_len(n_good))
	data.frame(
		unit = rep(good_units, each = 3L),
		time = rep(c(1L, 2L, 3L), times = n_good),
		treat = rep(c(0L, 0L, 0L), times = n_good),
		stringsAsFactors = FALSE
	)
}

# Capture conditionMessage from a stop() inside idCohorts.
.catch_idcohorts <- function(df) {
	tryCatch(
		fetwfe:::idCohorts(
			df = df,
			time_var = "time",
			unit_var = "unit",
			treat_var = "treat",
			covs = NULL
		),
		error = function(e) conditionMessage(e)
	)
}

# ------------------------------------------------------------------------------
# .truncate_violations() direct unit test
# ------------------------------------------------------------------------------

test_that(".truncate_violations returns xs unchanged when below max_show", {
	xs <- letters[1:5]
	expect_identical(fetwfe:::.truncate_violations(xs, max_show = 20L), xs)
})

test_that(".truncate_violations truncates and appends '... and N more' when above max_show", {
	xs <- as.character(1:25)
	out <- fetwfe:::.truncate_violations(xs, max_show = 20L)
	expect_equal(length(out), 21L)
	expect_equal(out[1:20], as.character(1:20))
	expect_equal(out[21], "... and 5 more")
})

test_that(".truncate_violations boundary: length equal to max_show returns unchanged", {
	xs <- letters[1:20]
	expect_identical(fetwfe:::.truncate_violations(xs, max_show = 20L), xs)
})

# ------------------------------------------------------------------------------
# Balance violations: single unit (sanity / wording-compat with v1.9.5).
# ------------------------------------------------------------------------------

test_that("idCohorts: single balance violation produces canonical-prefix message", {
	df <- .make_panel(n_good = 2L)
	# Drop one row of good01 → only 2 observations.
	df <- df[!(df$unit == "good01" & df$time == 3L), ]
	err <- .catch_idcohorts(df)
	expect_match(err, "Panel does not appear to be balanced")
	expect_match(err, "good01 has 2 observations", fixed = TRUE)
})

# ------------------------------------------------------------------------------
# Balance violations: multiple units — error message names ALL of them.
# ------------------------------------------------------------------------------

test_that("idCohorts: 3 balance violations reported together (lex-ordered)", {
	# Use lex-visible names ("unit10" < "unit2") to make the sort order
	# observable; mirrors the v1.9.3 lex-cohort-sort caveat (#53).
	df <- rbind(
		data.frame(
			unit = rep(c("unit2", "unit3", "unit10"), each = 3L),
			time = rep(c(1L, 2L, 3L), times = 3L),
			treat = 0L
		),
		data.frame(
			unit = "good01",
			time = 1L:3L,
			treat = 0L
		),
		stringsAsFactors = FALSE
	)
	# Truncate one row from each of unit2, unit3, unit10.
	df <- df[!(df$unit == "unit2" & df$time == 3L), ]
	df <- df[!(df$unit == "unit3" & df$time == 3L), ]
	df <- df[!(df$unit == "unit10" & df$time == 3L), ]
	err <- .catch_idcohorts(df)
	# All three names appear.
	expect_match(err, "unit2 has 2 observations", fixed = TRUE)
	expect_match(err, "unit3 has 2 observations", fixed = TRUE)
	expect_match(err, "unit10 has 2 observations", fixed = TRUE)
	# Ordering is lexicographic (unit10 < unit2 < unit3, the v1.9.3 trap).
	pos_10 <- regexpr("unit10 has", err, fixed = TRUE)
	pos_2 <- regexpr("unit2 has", err, fixed = TRUE)
	pos_3 <- regexpr("unit3 has", err, fixed = TRUE)
	expect_true(pos_10 < pos_2)
	expect_true(pos_2 < pos_3)
})

# ------------------------------------------------------------------------------
# Absorbing-state violations: multiple units.
# ------------------------------------------------------------------------------

test_that("idCohorts: multiple absorbing-state violations reported together", {
	# Build a balanced panel where two units have treatment that turns on
	# at time=2 and then off at time=3 (non-absorbing).
	df <- data.frame(
		unit = rep(c("absA", "absB", "good01"), each = 3L),
		time = rep(c(1L, 2L, 3L), times = 3L),
		treat = c(0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
		stringsAsFactors = FALSE
	)
	err <- .catch_idcohorts(df)
	expect_match(err, "Treatment does not appear to be an absorbing state")
	expect_match(err, "absA")
	expect_match(err, "absB")
	expect_false(grepl("good01", err, fixed = TRUE))
})

# ------------------------------------------------------------------------------
# Both kinds of violations simultaneously: error message has both sections.
# ------------------------------------------------------------------------------

test_that("idCohorts: balance + absorbing violations both surfaced in one message", {
	df <- data.frame(
		unit = rep(c("badA", "badB", "good01"), each = 3L),
		time = rep(c(1L, 2L, 3L), times = 3L),
		treat = c(0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
		stringsAsFactors = FALSE
	)
	# Drop one row from badA → balance violation.
	df <- df[!(df$unit == "badA" & df$time == 3L), ]
	err <- .catch_idcohorts(df)
	expect_match(err, "Panel does not appear to be balanced")
	expect_match(err, "badA has 2 observations", fixed = TRUE)
	expect_match(err, "Treatment does not appear to be an absorbing state")
	expect_match(err, "badB")
})

# ------------------------------------------------------------------------------
# Truncation: 25 balance violations → first 20 + "... and 5 more".
# ------------------------------------------------------------------------------

test_that("idCohorts: more than 20 balance violations are truncated to 20 + summary", {
	# Build 26 units, all balanced (3 rows each), then drop the third row
	# from 25 of them so they fail balance.
	n_bad <- 25L
	unit_names <- sprintf("u%03d", seq_len(n_bad + 1L))
	df <- data.frame(
		unit = rep(unit_names, each = 3L),
		time = rep(c(1L, 2L, 3L), times = length(unit_names)),
		treat = 0L,
		stringsAsFactors = FALSE
	)
	bad_units <- unit_names[seq_len(n_bad)]
	df <- df[!(df$unit %in% bad_units & df$time == 3L), ]
	err <- .catch_idcohorts(df)
	expect_match(err, "Panel does not appear to be balanced")
	# Summary tail.
	expect_match(err, "... and 5 more", fixed = TRUE)
	# Spot-check: the first lex-sorted bad unit is "u001"; the 20th is
	# "u020"; "u025" should NOT appear (it's in the truncated tail).
	expect_match(err, "u001 has 2 observations", fixed = TRUE)
	expect_match(err, "u020 has 2 observations", fixed = TRUE)
	expect_false(grepl("u025 has 2 observations", err, fixed = TRUE))
})
