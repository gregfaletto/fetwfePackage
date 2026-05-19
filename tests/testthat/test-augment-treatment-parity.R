# Tests for the augment column-set parity fix between the auto-trim and
# no-trim code paths in `.augment_estimator_output()` (R/broom_methods.R).
#
# Background: prior to v1.9.4, augment had two code paths. When `nrow(data)
# == nrow(X_ints)` (user pre-trimmed), `data` was used as-is and every user
# column flowed through to the output. When `nrow(data) != nrow(X_ints)`
# (auto-trim path), the function called `idCohorts()` and assigned
# `data <- ret$df`. But `idCohorts()` strips the `treat_var` column
# (R/utility.R:149), so the auto-trim output was missing the `treatment`
# column — asymmetric with the no-trim path and surprising for users who
# rely on `(unit, time, treatment)` downstream. See GitHub issue #51.

if (!requireNamespace("broom", quietly = TRUE)) {
	testthat::skip("broom not installed")
}

# ------------------------------------------------------------------------------
# Shared fixture: a panel with one first-period-treated unit, so the
# augment auto-trim path is exercised. We reuse the simulator regime from
# test-broom-methods.R rather than coupling to that file's helper.
# ------------------------------------------------------------------------------

.parity_fixture <- function(seed = 20260516) {
	coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = seed
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	pdata <- sim$pdata

	# Force one never-treated unit into the first-period-treated cohort so
	# augment's auto-trim path will fire. Treatment is absorbing, so set 1
	# for all periods of that unit.
	never_treated <- unique(
		pdata$unit[
			vapply(
				unique(pdata$unit),
				function(u) all(pdata$treatment[pdata$unit == u] == 0),
				logical(1)
			)
		]
	)
	testthat::skip_if(
		length(never_treated) == 0,
		"no never-treated units in fixture"
	)
	u_first <- never_treated[1]
	pdata$treatment[pdata$unit == u_first] <- 1L

	# Fit fetwfe on the full panel — the estimator internally drops the
	# first-period-treated unit via idCohorts(), so X_ints has (N - 1) * T
	# rows.
	res <- suppressWarnings(fetwfe(
		pdata = pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		q = 0.5,
		verbose = FALSE
	))

	# Pre-trimmed version: the same panel with the offending unit removed.
	# nrow(pdata_pre) == res$N * res$T == nrow(X_ints), so augment's
	# auto-trim path is skipped on this input.
	pdata_pre <- pdata[pdata$unit != u_first, , drop = FALSE]

	list(res = res, pdata = pdata, pdata_pre = pdata_pre, dropped = u_first)
}

test_that("augment auto-trim and no-trim paths return identical column sets", {
	fx <- .parity_fixture()

	aug_auto <- suppressWarnings(broom::augment(fx$res, data = fx$pdata))
	aug_pre <- broom::augment(fx$res, data = fx$pdata_pre)

	# Column-set parity: both paths must return exactly the same columns.
	expect_equal(sort(names(aug_auto)), sort(names(aug_pre)))

	# The specific bug: treatment column must survive the auto-trim path.
	expect_true("treatment" %in% names(aug_auto))
	expect_true("treatment" %in% names(aug_pre))

	# Sanity: the auto-trim path still produces correct row counts and
	# .fitted / .resid behaves correctly.
	expect_equal(nrow(aug_auto), fx$res$N * fx$res$T)
	expect_false(fx$dropped %in% aug_auto$unit)
	# Round-trip identity (.fitted + .resid == aug_auto$y) was removed per
	# #76 Item 7b -- it's tautological (holds by construction of the .resid
	# subtraction regardless of the auto-trim's row-alignment correctness).
	# The load-bearing row-order check lives in the
	# "augment column parity holds ..." test_that block below
	# (state_name / panel_id column alignment after sort). See
	# feedback_tautological_roundtrip_tests.md for the broader pattern.
})

test_that("augment column parity holds for etwfe and betwfe on auto-trim path", {
	# `.augment_estimator_output()` is the shared body for augment.fetwfe,
	# augment.etwfe, and augment.betwfe — so the fix in v1.9.4 applies to all
	# three. This block locks in that behavior so a future class-conditional
	# regression that only breaks etwfe/betwfe doesn't slip through.
	fx <- .parity_fixture()

	res_etwfe <- suppressWarnings(etwfe(
		pdata = fx$pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		verbose = FALSE
	))
	aug_etwfe_auto <- suppressWarnings(broom::augment(
		res_etwfe,
		data = fx$pdata
	))
	aug_etwfe_pre <- broom::augment(res_etwfe, data = fx$pdata_pre)
	expect_equal(sort(names(aug_etwfe_auto)), sort(names(aug_etwfe_pre)))
	expect_true("treatment" %in% names(aug_etwfe_auto))
	expect_false(fx$dropped %in% aug_etwfe_auto$unit)

	res_betwfe <- suppressWarnings(betwfe(
		pdata = fx$pdata,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		q = 0.5,
		verbose = FALSE
	))
	aug_betwfe_auto <- suppressWarnings(broom::augment(
		res_betwfe,
		data = fx$pdata
	))
	aug_betwfe_pre <- broom::augment(res_betwfe, data = fx$pdata_pre)
	expect_equal(sort(names(aug_betwfe_auto)), sort(names(aug_betwfe_pre)))
	expect_true("treatment" %in% names(aug_betwfe_auto))
	expect_false(fx$dropped %in% aug_betwfe_auto$unit)
})

test_that("augment auto-trim path preserves extra user columns", {
	fx <- .parity_fixture()

	# Add columns the estimator does not know about. These should flow
	# through both augment paths unchanged. `state_name` is character;
	# `panel_id` is integer — exercising both types. `unit` is character,
	# so we build panel_id via match() against a stable level table rather
	# than coercing unit names to integer directly.
	pdata <- fx$pdata
	pdata$state_name <- paste0("state_", pdata$unit)
	unit_levels <- unique(pdata$unit)
	pdata$panel_id <- match(pdata$unit, unit_levels) * 10L

	aug_auto <- suppressWarnings(broom::augment(fx$res, data = pdata))

	expect_true("state_name" %in% names(aug_auto))
	expect_true("panel_id" %in% names(aug_auto))
	# Values still align with the (unit, time) sort.
	expect_equal(aug_auto$state_name, paste0("state_", aug_auto$unit))
	expect_equal(aug_auto$panel_id, match(aug_auto$unit, unit_levels) * 10L)

	# The dropped first-period-treated unit must not appear.
	expect_false(fx$dropped %in% aug_auto$unit)
	expect_false(paste0("state_", fx$dropped) %in% aug_auto$state_name)
})
