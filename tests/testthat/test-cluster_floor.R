library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Tests for `.floor_cluster_quad()` and its integration at the five
# cluster-sandwich quadratic-form floor sites (issue #139, version 1.11.2).
#
# `.floor_cluster_quad()` layers a two-tier diagnostic on top of the
# pre-existing `max(q, 0)` floor at each cluster-sandwich quadratic-form
# site. The forms are PSD in exact arithmetic, so any negative value is
# either FP-noise (silently clipped to 0) or a bug signal (warn / stop).
#
# The five sites:
#   1. R/variance_machinery.R    getTeResultsOLS         (att_var_1)
#   2. R/variance_machinery.R    getTeResults2           (att_var_1)
#   3. R/variance_machinery.R    getCohortATTsFinal      (cohort_te_ses)
#   4. R/event_study.R           .event_study_etwfe_betwfe (var_1_e)
#   5. R/event_study.R           .event_study_fetwfe       (var_1_e)
#
# (PR #111 added 6 sites in 3 files; PR #118 consolidated
#  getCohortATTsFinalOLS into getCohortATTsFinal, dropping the count to 5.
#  Surprises & Discoveries in .plans/cluster-sandwich-floor-diagnostics-139
#  documents the 6 -> 5 reduction.)
# ------------------------------------------------------------------------------

# --- Unit tests on the helper -------------------------------------------------

test_that(".floor_cluster_quad returns positive values unchanged", {
	expect_equal(fetwfe:::.floor_cluster_quad(0.5, "test"), 0.5)
	expect_equal(fetwfe:::.floor_cluster_quad(1e-12, "test"), 1e-12)
	expect_equal(fetwfe:::.floor_cluster_quad(100, "test"), 100)
})

test_that(".floor_cluster_quad returns 0 unchanged", {
	expect_equal(fetwfe:::.floor_cluster_quad(0, "test"), 0)
})

test_that(".floor_cluster_quad silently floors FP-noise negatives to 0", {
	# Magnitudes inside the FP-noise band ([-1e-10, 0]) clip without warning.
	expect_silent(out1 <- fetwfe:::.floor_cluster_quad(-1e-15, "test"))
	expect_equal(out1, 0)
	expect_silent(out2 <- fetwfe:::.floor_cluster_quad(-1e-11, "test"))
	expect_equal(out2, 0)
	# Exactly at the warning threshold is still silent (strict `<`).
	expect_silent(out3 <- fetwfe:::.floor_cluster_quad(-1e-10, "test"))
	expect_equal(out3, 0)
})

test_that(".floor_cluster_quad warns on warning-range negatives", {
	expect_warning(
		out <- fetwfe:::.floor_cluster_quad(-1e-8, "test_site"),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out, 0)
	# Site name in the message.
	expect_warning(
		fetwfe:::.floor_cluster_quad(-1e-8, "test_site"),
		"test_site"
	)
	# Boundary: just-inside the error threshold still warns, doesn't stop.
	expect_warning(
		out_boundary <- fetwfe:::.floor_cluster_quad(-0.999, "test"),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out_boundary, 0)
})

test_that(".floor_cluster_quad errors on catastrophic negatives", {
	expect_error(
		fetwfe:::.floor_cluster_quad(-2, "test_site"),
		"catastrophically negative"
	)
	expect_error(
		fetwfe:::.floor_cluster_quad(-1.5, "another_site"),
		"another_site"
	)
	expect_error(
		fetwfe:::.floor_cluster_quad(-1000, "test"),
		"file an issue"
	)
})

test_that(".floor_cluster_quad pass-through on NA / multi-element / non-numeric", {
	# NA passes through.
	expect_silent(out_na <- fetwfe:::.floor_cluster_quad(NA_real_, "test"))
	expect_true(is.na(out_na))
	# Multi-element vector passes through unchanged.
	expect_silent(out_vec <- fetwfe:::.floor_cluster_quad(c(-2, 1), "test"))
	expect_equal(out_vec, c(-2, 1))
	# Non-numeric (character) passes through.
	expect_silent(
		out_chr <- fetwfe:::.floor_cluster_quad("not_numeric", "test")
	)
	expect_equal(out_chr, "not_numeric")
	# Length-0 numeric passes through.
	expect_silent(out_e <- fetwfe:::.floor_cluster_quad(numeric(0), "test"))
	expect_equal(out_e, numeric(0))
})

test_that(".floor_cluster_quad respects custom thresholds", {
	# Tighter warn threshold -> a value that would normally pass silently warns.
	expect_warning(
		out <- fetwfe:::.floor_cluster_quad(
			-1e-12,
			"test",
			warn_threshold = -1e-13
		),
		"Negative cluster-sandwich quadratic form"
	)
	expect_equal(out, 0)
	# Looser error threshold -> a value that would normally warn errors.
	expect_error(
		fetwfe:::.floor_cluster_quad(
			-0.5,
			"test",
			err_threshold = -0.1
		),
		"catastrophically negative"
	)
})

# --- Coverage / regression test ----------------------------------------------
# Assert every cluster-sandwich quadratic-form site in the package code
# routes through `.floor_cluster_quad`, catching a future revert to the bare
# `max(., 0)` floor.

test_that("all 5 cluster-sandwich floor sites use .floor_cluster_quad", {
	pkg_root <- system.file(package = "fetwfe")
	# When the package is loaded via devtools::load_all(), the installed
	# `system.file()` may return the install dir (no R/ source). Prefer the
	# source tree if it is present (the test repo root).
	src_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
	if (
		dir.exists(file.path(src_root, "R")) &&
			file.exists(file.path(src_root, "R", "variance_machinery.R"))
	) {
		root <- src_root
	} else {
		root <- pkg_root
	}

	files <- c("R/variance_machinery.R", "R/event_study.R")
	# Count `.floor_cluster_quad(` call sites across both files. The helper
	# is invoked once per site, so the call count is exactly 5. We skip
	# comment lines (so the matching helper-name mention in each site's
	# explanatory comment does not double-count).
	call_count <- 0L
	for (f in files) {
		path <- file.path(root, f)
		skip_if_not(
			file.exists(path),
			paste("source file not in test bundle:", f)
		)
		src <- readLines(path, warn = FALSE)
		non_comment <- src[!grepl("^\\s*#", src)]
		call_count <- call_count +
			sum(grepl(".floor_cluster_quad(", non_comment, fixed = TRUE))
	}
	expect_equal(call_count, 5L)

	# Per-site label coverage: each site uses its expected per-site label.
	# This catches a future copy-paste regression where two sites end up
	# with the same label.
	expected_labels <- c(
		"getTeResultsOLS/att_var_1",
		"getTeResults2/att_var_1",
		"getCohortATTsFinal/cohort_te_se",
		"event_study_etwfe_betwfe/var_1_e",
		"event_study_fetwfe/var_1_e"
	)
	src_all <- character(0)
	for (f in files) {
		path <- file.path(root, f)
		if (file.exists(path)) {
			src_all <- c(src_all, readLines(path, warn = FALSE))
		}
	}
	for (lbl in expected_labels) {
		expect_true(
			any(grepl(lbl, src_all, fixed = TRUE)),
			info = paste("missing site label:", lbl)
		)
	}

	# Negative-direction guard: no bare `max(..., 0)` floor remains around a
	# `sandwich_full` quadratic form.
	bare_max_on_sandwich <- 0L
	for (f in files) {
		path <- file.path(root, f)
		if (!file.exists(path)) {
			next
		}
		src <- readLines(path, warn = FALSE)
		# Look for the historical pattern: a `max(` line followed shortly by
		# `sandwich_full %*%` and a closing `, 0)`. The wired-up helper does
		# not use `max(`, so this catches a partial revert.
		for (i in seq_along(src)) {
			if (grepl("^\\s*max\\(", src[i])) {
				# Inspect the next 6 lines for the sandwich pattern.
				window <- src[i:min(length(src), i + 6L)]
				if (
					any(grepl("sandwich_full %*%", window, fixed = TRUE)) &&
						any(grepl(", 0\\)", window))
				) {
					bare_max_on_sandwich <- bare_max_on_sandwich + 1L
				}
			}
		}
	}
	expect_equal(bare_max_on_sandwich, 0L)
})

# --- Integration smoke test --------------------------------------------------
# Fit a small cluster-SE model on well-conditioned simulated data and
# verify no warning/error from `.floor_cluster_quad` fires.

test_that("cluster-SE fit on well-conditioned data does not trigger the diagnostic", {
	# Small but well-conditioned: same recipe as several existing tests.
	set.seed(2026)
	coefs <- genCoefs(
		R = 3,
		T = 5,
		density = 0.5,
		eff_size = 1,
		d = 1,
		seed = 2026
	)
	sim <- simulateData(
		coefs,
		N = 80,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	# Run the fit; assert no warning at all surfaces with the message
	# pattern from `.floor_cluster_quad()`. (Other unrelated warnings are
	# not the concern of this test; we filter by substring.)
	withCallingHandlers(
		res <- fetwfeWithSimulatedData(sim, q = 0.5, se_type = "cluster"),
		warning = function(w) {
			msg <- conditionMessage(w)
			if (grepl("cluster-sandwich quadratic form", msg, fixed = TRUE)) {
				stop(
					"Unexpected .floor_cluster_quad warning on well-conditioned data: ",
					msg
				)
			}
			invokeRestart("muffleWarning")
		}
	)
	expect_true(!is.null(res$att_se))
	expect_true(is.finite(res$att_se))
})
