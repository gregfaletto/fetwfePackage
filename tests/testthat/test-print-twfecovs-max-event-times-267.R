library(testthat)
library(fetwfe)

# Panel-data fixture helpers (generate_panel_data, ...) live in
# tests/testthat/helper-panel-fixture.R and are sourced by testthat first (#91).

# ------------------------------------------------------------------------------
# #267: print.twfeCovs() must accept `max_event_times` like its three sibling print
# methods, instead of erroring via a formals/body mismatch. The body forwarded a
# literal `max_event_times = 10L` AND `...`, so a user-supplied value double-bound
# into .print_estimator_output() ("matched by multiple actual arguments"). The arg
# is a documented no-op for twfeCovs (no event-study section), so the correct
# behavior is "accepted and ignored", not an error.
# ------------------------------------------------------------------------------
test_that("print.twfeCovs() accepts max_event_times without erroring (#267)", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	tc <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	# The formal now exists, matching the three sibling print methods.
	expect_true(
		"max_event_times" %in% names(formals(getS3method("print", "twfeCovs")))
	)
	# A user-supplied value no longer double-binds into .print_estimator_output().
	expect_no_error(invisible(capture.output(print(tc, max_event_times = 5))))
	# Sibling args and the plain call still work.
	expect_no_error(invisible(capture.output(print(tc, max_cohorts = 2))))
	expect_no_error(invisible(capture.output(print(tc))))
})
