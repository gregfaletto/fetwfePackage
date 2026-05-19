library(testthat)
library(fetwfe)

# Snapshot guardrail for the four estimator entry points' input-validation
# error messages. Locks the canonical `<predicate> is not TRUE` strings
# produced by `checkEtwfeInputs()` (and `checkFetwfeInputs()` for the
# fetwfe-specific extension) against six markdown goldens under
# tests/testthat/_snaps/estimator-entry-snapshot.md so any byte drift
# surfaces as a failing test. Pre-refactor guardrail for issue #79's
# entry-point consolidation: PR B will eliminate the verbatim 12-step
# pipeline duplication, and any accidental error-message rewording will
# trip this file.
#
# Each test mutates one field of a known-good 4-unit x 3-period panel
# to trigger one canonical malformed-input pattern, then runs all 4
# entry points (or 2, for the fetwfe/betwfe-specific `q` test) and
# asserts byte-equality of the captured `conditionMessage()`. The
# canonical message is then snapshotted so it can't drift without
# explicit re-blessing.
#
# Snapshot tests require testthat edition 3 and skip under CRAN by
# default (so `R CMD check --as-cran` shows SKIPs for these tests;
# `devtools::test()` runs them).

testthat::local_edition(3)

# In-file panel fixture (NOT shared with other test files; this PR's
# scope is intentionally narrow). 4 units x 3 time periods = 12 rows.
# U1 never treated, U2/U3 in cohort 2, U4 in cohort 3.
.build_valid_panel <- function() {
	data.frame(
		unit = rep(c("U1", "U2", "U3", "U4"), each = 3),
		time = rep(1L:3L, 4),
		treat = c(
			0L,
			0L,
			0L, # U1: never
			0L,
			1L,
			1L, # U2: cohort 2
			0L,
			1L,
			1L, # U3: cohort 2
			0L,
			0L,
			1L # U4: cohort 3
		),
		y = as.numeric(1:12),
		stringsAsFactors = FALSE
	)
}

# Capture conditionMessage() from each entry point under a given
# malformed call. Returns a list with one entry per estimator name.
.capture_entry_point_errors <- function(
	args,
	estimators = c(
		"fetwfe",
		"etwfe",
		"betwfe",
		"twfeCovs"
	)
) {
	fns <- list(
		fetwfe = fetwfe,
		etwfe = etwfe,
		betwfe = betwfe,
		twfeCovs = twfeCovs
	)
	out <- list()
	for (est in estimators) {
		out[[est]] <- tryCatch(
			do.call(fns[[est]], args),
			error = function(e) conditionMessage(e)
		)
	}
	out
}

# ------------------------------------------------------------------------------
# Test 1: wrong column type. `pdata[[time_var]]` is numeric (not integer).
# Trips `stopifnot(is.integer(pdata[[time_var]]))` in checkEtwfeInputs.
# ------------------------------------------------------------------------------
test_that("entry-point error: pdata[[time_var]] not integer (#79 PR A)", {
	pdata <- .build_valid_panel()
	pdata$time <- as.numeric(pdata$time) # cast away the integer marker
	args <- list(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treat",
		response = "y"
	)
	msgs <- .capture_entry_point_errors(args)

	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 2: missing column. `time_var` references a column that's not in pdata.
# Trips `stopifnot(time_var %in% colnames(pdata))`.
# ------------------------------------------------------------------------------
test_that("entry-point error: time_var not a column (#79 PR A)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = "nonexistent",
		unit_var = "unit",
		treatment = "treat",
		response = "y"
	)
	msgs <- .capture_entry_point_errors(args)

	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 3: non-character flag arg. `time_var` is an integer.
# Trips `stopifnot(is.character(time_var))`.
# ------------------------------------------------------------------------------
test_that("entry-point error: time_var not character (#79 PR A)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = 42L,
		unit_var = "unit",
		treatment = "treat",
		response = "y"
	)
	msgs <- .capture_entry_point_errors(args)

	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 4: non-logical flag arg. `verbose` is a character.
# Trips `stopifnot(is.logical(verbose))`.
# ------------------------------------------------------------------------------
test_that("entry-point error: verbose not logical (#79 PR A)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		verbose = "yes"
	)
	msgs <- .capture_entry_point_errors(args)

	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 5: range violation. `alpha` is negative.
# Trips `stopifnot(alpha > 0)`.
# ------------------------------------------------------------------------------
test_that("entry-point error: alpha out of range (#79 PR A)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		alpha = -0.1
	)
	msgs <- .capture_entry_point_errors(args)

	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 6: fetwfe / betwfe `q` violation. Both call checkFetwfeInputs(), so
# both should produce a byte-identical "q > 0 is not TRUE" message.
# etwfe / twfeCovs don't have a `q` arg and are excluded.
# ------------------------------------------------------------------------------
test_that("entry-point error: fetwfe + betwfe q out of range (#79 PR A)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		q = -1
	)
	msgs <- .capture_entry_point_errors(
		args,
		estimators = c("fetwfe", "betwfe")
	)

	expect_identical(msgs$fetwfe, msgs$betwfe)

	expect_snapshot(cat(msgs$fetwfe))
})

# ------------------------------------------------------------------------------
# Test 7 (#84): NULL-input handling. Plan-review explicitly called out
# `class(NULL)[1] = "NULL"`, `length(NULL) = 0` so the per-arg sprintf
# interpolation should produce a clean, non-crashing message. We
# exercise two NULL args at once (`time_var`, `verbose`) and assert
# that BOTH violations appear (collect-across-args) in a single error
# message — this is the load-bearing behavior change for #84.
# ------------------------------------------------------------------------------
test_that("entry-point error: NULL inputs produce clean multi-arg message (#84)", {
	pdata <- .build_valid_panel()
	args <- list(
		pdata = pdata,
		time_var = NULL,
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		verbose = NULL
	)
	msgs <- .capture_entry_point_errors(args)

	# 4-way byte-equality still holds across the entry points.
	expect_identical(msgs$fetwfe, msgs$etwfe)
	expect_identical(msgs$etwfe, msgs$betwfe)
	expect_identical(msgs$betwfe, msgs$twfeCovs)

	# Both NULL args contribute their own bullet (collect-across-args).
	expect_match(msgs$fetwfe, "time_var must be a single character")
	expect_match(msgs$fetwfe, "verbose must be a single logical")
	# Sanity: message header is the new "Invalid inputs:" header.
	expect_match(msgs$fetwfe, "^Invalid inputs:")
	# Sanity: the message references NULL (the class of a NULL arg)
	# and length 0 — confirms sprintf interpolation didn't crash.
	expect_match(msgs$fetwfe, "NULL")
	expect_match(msgs$fetwfe, "length 0")
})
