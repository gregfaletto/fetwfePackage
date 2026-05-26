library(testthat)
library(fetwfe)

# Internal-slot parity (#144): all four estimators expose `$internal` with
# the canonical inner slots. For the three OLS-family estimators
# (etwfe / betwfe / twfeCovs), the same five inner slots are also
# duplicated at top level for backward compat ("pure additive" strategy);
# `$internal` is the canonical access path going forward and the only
# location for these slots on `fetwfe` results.

.parity_setup <- function() {
	set.seed(2026)
	coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	dat
}

# ----------------------------------------------------------------------
# 1) All four estimators expose `$internal` with the canonical inner
#    slots. fetwfe additionally carries `theta_hat` (bridge-selection
#    parameter); the OLS-family does not.
# ----------------------------------------------------------------------

test_that("all four estimators expose `$internal` with the canonical inner slots", {
	sim <- .parity_setup()
	res_fe <- fetwfeWithSimulatedData(sim)
	res_et <- etwfeWithSimulatedData(sim)
	res_be <- betwfeWithSimulatedData(sim)
	res_tc <- twfeCovsWithSimulatedData(sim)

	expected_ols <- c("X_ints", "y", "X_final", "y_final", "calc_ses")
	expected_fe <- c(expected_ols, "theta_hat")

	expect_true(!is.null(res_fe$internal))
	expect_true(all(expected_fe %in% names(res_fe$internal)))

	for (res in list(res_et, res_be, res_tc)) {
		expect_true(!is.null(res$internal))
		expect_true(all(expected_ols %in% names(res$internal)))
	}
})

# ----------------------------------------------------------------------
# 2) Pure additive: for the OLS-family, `$internal` slots equal the
#    top-level slots (same value, different access path). Regression
#    guard against silent drift between the two paths.
# ----------------------------------------------------------------------

test_that("OLS-family: $internal slots equal top-level slots (parity)", {
	sim <- .parity_setup()
	for (res in list(
		etwfeWithSimulatedData(sim),
		betwfeWithSimulatedData(sim),
		twfeCovsWithSimulatedData(sim)
	)) {
		expect_identical(res$internal$X_ints, res$X_ints)
		expect_identical(res$internal$y, res$y)
		expect_identical(res$internal$X_final, res$X_final)
		expect_identical(res$internal$y_final, res$y_final)
		expect_identical(res$internal$calc_ses, res$calc_ses)
	}
})

# ----------------------------------------------------------------------
# 3) fetwfe: the inner slots are EXCLUSIVELY under `$internal` (existing
#    behavior, regression guard). Documents the asymmetry that #144
#    consciously preserves: fetwfe gates these slots through
#    `$internal` as the only access path; the OLS-family duplicates for
#    backward compat.
# ----------------------------------------------------------------------

test_that("fetwfe: $internal slots are NOT duplicated at top level", {
	sim <- .parity_setup()
	res <- fetwfeWithSimulatedData(sim)
	expect_true(!is.null(res$internal$X_ints))
	# Use `[[` (exact match) -- `$` does partial matching, so `res$y`
	# would match `res$y_mean` (a real top-level slot) and produce a
	# false positive here.
	expect_null(res[["X_ints"]])
	expect_null(res[["y"]])
	expect_null(res[["X_final"]])
	expect_null(res[["y_final"]])
	expect_null(res[["theta_hat"]])
	expect_null(res[["calc_ses"]])
})

# ----------------------------------------------------------------------
# 4) Validators enforce the new contract: if `$internal` is dropped from
#    the constructed object, the validator stops with a structured
#    error.
# ----------------------------------------------------------------------

test_that("validators enforce `$internal` presence for OLS-family", {
	sim <- .parity_setup()
	for (entry in list(
		list(res = etwfeWithSimulatedData(sim), fn = fetwfe:::.validate_etwfe),
		list(
			res = betwfeWithSimulatedData(sim),
			fn = fetwfe:::.validate_betwfe
		),
		list(
			res = twfeCovsWithSimulatedData(sim),
			fn = fetwfe:::.validate_twfeCovs
		)
	)) {
		bad <- entry$res
		bad$internal <- NULL
		expect_error(
			entry$fn(bad),
			"Missing slot.*internal"
		)
	}
})

test_that("validators enforce inner-slot presence for OLS-family", {
	sim <- .parity_setup()
	for (entry in list(
		list(res = etwfeWithSimulatedData(sim), fn = fetwfe:::.validate_etwfe),
		list(
			res = betwfeWithSimulatedData(sim),
			fn = fetwfe:::.validate_betwfe
		),
		list(
			res = twfeCovsWithSimulatedData(sim),
			fn = fetwfe:::.validate_twfeCovs
		)
	)) {
		bad <- entry$res
		bad$internal$X_ints <- NULL
		expect_error(
			entry$fn(bad),
			"Missing slot.*X_ints"
		)
	}
})
