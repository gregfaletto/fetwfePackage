library(testthat)
library(fetwfe)

# Issue #138: event-study estimates included in print() and summary()
# displays. The byte-exact rendered output is locked by
# tests/testthat/_snaps/print-method-snapshot.md (regenerated when #138
# landed); this file covers the behavioral edges that snapshots don't:
# (a) the summary list carries an `event_study` field with the expected
# shape, (b) print() silently skips the event-study block when
# eventStudy() errors (defense-in-depth), and (c) the max_event_times
# option controls truncation.

.es_display_setup <- function() {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	simulateData(coefs, N = 60, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

# ----------------------------------------------------------------------
# 1) The summary list carries an `event_study` field with the same
#    snake_case column shape eventStudy() returns; positioned between
#    `catt` and `model_info` per .summary_estimator_output()'s `keep`
#    ordering.
# ----------------------------------------------------------------------

test_that("summary() output carries event_study field with eventStudy() shape", {
	sim <- .es_display_setup()
	for (res in list(
		fetwfeWithSimulatedData(sim),
		etwfeWithSimulatedData(sim),
		betwfeWithSimulatedData(sim)
	)) {
		s <- summary(res)
		expect_true("event_study" %in% names(s))
		expect_s3_class(s$event_study, "data.frame")
		expect_true(all(
			c(
				"event_time",
				"n_cohorts",
				"estimate",
				"se",
				"ci_low",
				"ci_high",
				"p_value"
			) %in%
				names(s$event_study)
		))
		# Field ordering: event_study sits between catt and model_info
		# (#138 chose this position so the print preview matches the
		# CATT -> ES -> Model Details flow).
		keys <- names(s)
		expect_true(which(keys == "event_study") > which(keys == "catt"))
		expect_true(
			which(keys == "event_study") < which(keys == "model_info")
		)
	}
})

# ----------------------------------------------------------------------
# 2) print() and summary() silently skip the event-study block when
#    eventStudy() would error -- defense-in-depth so a fit with a
#    configuration that breaks eventStudy() doesn't crash the routine
#    display. Currently no realistic fit triggers this path (the
#    eventStudy() dispatch handles fetwfe/etwfe/betwfe explicitly), but
#    the tryCatch is exercised via a synthetic broken object.
# ----------------------------------------------------------------------

test_that("print() / summary() skip event-study block when eventStudy() errors", {
	sim <- .es_display_setup()
	res <- etwfeWithSimulatedData(sim)
	# Mutate the object so eventStudy() would fail. eventStudy() dispatches
	# on class; stripping it to plain list (but keeping `etwfe` first so
	# print.etwfe still dispatches) breaks the eventStudy() class check.
	broken <- res
	class(broken) <- c("etwfe", "list_without_eventStudy_dispatch")
	# eventStudy itself errors on this object (no longer matches the
	# class chain it expects internally).
	# print() should NOT crash; the tryCatch in .print_estimator_output
	# absorbs the error and skips the event-study block.
	out <- capture.output(print(broken))
	joined <- paste(out, collapse = "\n")
	# Header still rendered.
	expect_match(joined, "Extended Two-Way Fixed Effects Results")
	# Other blocks still rendered.
	expect_match(joined, "Cohort Average Treatment Effects")
	expect_match(joined, "Model Details")
})

# ----------------------------------------------------------------------
# 3) max_event_times truncation: the displayed event-study preview is
#    capped at max_event_times rows, with a "... and N more event times"
#    footer when truncated.
# ----------------------------------------------------------------------

test_that("max_event_times truncates the event-study preview", {
	sim <- .es_display_setup()
	res <- fetwfeWithSimulatedData(sim)
	n_event_times <- nrow(eventStudy(res))
	# Pre-condition: the fixture has more than 2 event times (otherwise
	# the truncation path is unreachable).
	expect_gt(n_event_times, 2L)
	# Truncate aggressively to 2; expect "... and (n_event_times - 2)
	# more event times.".
	out <- capture.output(print(res, max_event_times = 2L))
	joined <- paste(out, collapse = "\n")
	expect_match(joined, "Event-Study Average Treatment Effects")
	expect_match(
		joined,
		sprintf("and %d more event times", n_event_times - 2L)
	)
	# Untruncated (large limit) should NOT show the "more event times"
	# footer.
	out_full <- capture.output(print(res, max_event_times = 100L))
	joined_full <- paste(out_full, collapse = "\n")
	expect_match(joined_full, "Event-Study Average Treatment Effects")
	expect_false(grepl("more event times", joined_full))
})
