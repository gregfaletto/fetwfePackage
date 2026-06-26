# Regression tests for the behavior-touching items of the 2026-06-19 review
# cleanup sweep (#325). The other items (theorem-anchor corrections, the dead
# `# model,` comments, the riesz_lasso feasibility-semantics comment, the
# col_ss/ses degeneracy comment) are comment-only and carry no testable behavior.

test_that("eventStudy() dispatch error echoes the offending class (#325 item 4)", {
	# Matches the cohortStudy() sibling pattern: the message names the class it got.
	expect_error(
		eventStudy(42L),
		"requires an object of class.*Got class: integer",
		fixed = FALSE
	)
	expect_error(
		eventStudy(structure(list(), class = "lm")),
		"Got class: lm",
		fixed = FALSE
	)
	# a genuinely fitted object still dispatches (no error from the guard).
	# (covered elsewhere; here we only pin the failure message.)
})

test_that(".validate_in_sample_counts_named_unique enforces named + unique (#325 item 3)", {
	validate <- fetwfe:::.validate_in_sample_counts_named_unique
	msg <- "must have all unique named entries"

	# valid: fully named, unique -> silent (returns invisibly).
	expect_silent(validate(c(a = 5L, b = 3L, c = 2L)))
	expect_null(validate(c(a = 5L, b = 3L)))

	# unnamed -> error (names() is NULL, length 0 != length(x)). This and the
	# duplicate-names case are exactly what the original two inline blocks caught;
	# a single empty name among otherwise-named entries is NOT caught by either the
	# old code or the helper (byte-faithful refactor), so it is not asserted here.
	expect_error(validate(c(5L, 3L, 2L)), msg)
	# duplicate names -> error.
	expect_error(validate(c(a = 5L, a = 3L)), msg)
})
