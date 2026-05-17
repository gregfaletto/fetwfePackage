# Tests for the constructor validators introduced in #85.
#
# Three sections:
#   (a) well-formed objects pass the validator (covers q in {0.5, 1, 2}
#       and the all-zero-theta fallback);
#   (b) each contract C1-C8 catches its target malformation with a clear
#       error message naming the failed contract;
#   (c) cross-class consistency: the validator's expected-slots list
#       matches `names(<live fit>)` via expect_setequal (set membership;
#       order doesn't matter);
#   (d) slot-deletion mutation: setting key slots to NULL is caught.
#
# Mutation tests for (b) deliberately malform an object and assert the
# validator errors with a message containing the contract name.

# ------------------------------------------------------------------------------
# Shared fixture
# ------------------------------------------------------------------------------

.validator_fixture <- function(q = 0.5, seed = 17052026) {
	sim <- simulateData(
		genCoefs(R = 2, T = 3, d = 2, density = 0.5, eff_size = 1, seed = seed),
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	list(
		sim = sim,
		fetwfe = suppressWarnings(fetwfeWithSimulatedData(sim, q = q)),
		etwfe = suppressWarnings(etwfeWithSimulatedData(sim)),
		betwfe = suppressWarnings(betwfeWithSimulatedData(sim, q = q)),
		twfeCovs = suppressWarnings(twfeCovsWithSimulatedData(sim))
	)
}

# ------------------------------------------------------------------------------
# (a) Well-formed objects pass
# ------------------------------------------------------------------------------

test_that("validators accept well-formed objects across q in {0.5, 1, 2}", {
	# q = 0.5: bridge actively selects (default case)
	fx05 <- .validator_fixture(q = 0.5)
	expect_silent(fetwfe:::.validate_fetwfe(fx05$fetwfe))
	expect_silent(fetwfe:::.validate_etwfe(fx05$etwfe))
	expect_silent(fetwfe:::.validate_betwfe(fx05$betwfe))
	expect_silent(fetwfe:::.validate_twfeCovs(fx05$twfeCovs))

	# q = 1: lasso, calc_ses = FALSE (the #73-relevant path).
	# Only exercises fetwfe and betwfe (etwfe/twfeCovs have no q).
	fx1 <- .validator_fixture(q = 1)
	expect_silent(fetwfe:::.validate_fetwfe(fx1$fetwfe))
	expect_silent(fetwfe:::.validate_betwfe(fx1$betwfe))
	# Confirm calc_ses is indeed FALSE on the q=1 path (the validator's
	# C1 contract is meaningfully exercised here).
	expect_false(isTRUE(fx1$fetwfe$internal$calc_ses))
	expect_false(isTRUE(fx1$betwfe$calc_ses))

	# q = 2: ridge, calc_ses = FALSE
	fx2 <- .validator_fixture(q = 2)
	expect_silent(fetwfe:::.validate_fetwfe(fx2$fetwfe))
	expect_silent(fetwfe:::.validate_betwfe(fx2$betwfe))
})

test_that("validators accept well-formed all-zero-theta fallback objects", {
	# A degenerate panel where the bridge selects nothing for fetwfe.
	# Use small N + tiny effects to maximize chance of full selection-out.
	sim_small <- simulateData(
		genCoefs(
			R = 2,
			T = 3,
			d = 0,
			density = 0.001,
			eff_size = 0.001,
			seed = 42
		),
		N = 30,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- suppressWarnings(fetwfeWithSimulatedData(sim_small, q = 0.5))
	# Confirm this IS the all-zero-theta path (skip if seed picked a
	# regime where bridge selected something).
	if (isTRUE(res$att_hat == 0)) {
		expect_false(isTRUE(res$att_selected))
		expect_silent(fetwfe:::.validate_fetwfe(res))
	} else {
		skip("seed didn't trigger all-zero fallback; smoke-test only")
	}
})

# ------------------------------------------------------------------------------
# (b) Contract violations are caught with clear messages
# ------------------------------------------------------------------------------

test_that("C1 SE-consistency: validator catches att_se non-NA when calc_ses is FALSE", {
	fx <- .validator_fixture(q = 1)
	res <- fx$fetwfe
	# Sanity check: this IS a calc_ses=FALSE object.
	stopifnot(!isTRUE(res$internal$calc_ses))
	stopifnot(is.na(res$att_se))
	# Mutate: force att_se to a non-NA value while calc_ses stays FALSE.
	res$att_se <- 1.5
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C1 SE consistency",
		fixed = TRUE
	)
})

test_that("C1 SE-consistency: validator catches catt_ses non-NA when calc_ses is FALSE", {
	fx <- .validator_fixture(q = 1)
	res <- fx$fetwfe
	stopifnot(!isTRUE(res$internal$calc_ses))
	stopifnot(all(is.na(res$catt_ses)))
	# Mutate: force one catt_ses entry non-NA.
	res$catt_ses[1] <- 0.1
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C1 SE consistency",
		fixed = TRUE
	)
})

test_that("C2 selection-consistency: validator catches inconsistent att_selected", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	stopifnot(isTRUE(res$att_selected))
	stopifnot(res$att_hat != 0)
	# Mutate: flip att_selected to FALSE while att_hat is nonzero.
	res$att_selected <- FALSE
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C2 selection consistency",
		fixed = TRUE
	)
})

test_that("C3 p-value NA: validator catches p_value non-NA when att_se is NA", {
	fx <- .validator_fixture(q = 1)
	res <- fx$fetwfe
	stopifnot(is.na(res$att_se))
	stopifnot(is.na(res$att_p_value))
	# Mutate: force att_p_value non-NA while att_se stays NA.
	res$att_p_value <- 0.05
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C3 p-value NA",
		fixed = TRUE
	)
})

test_that("C4 catt_df shape: validator catches wrong nrow", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	stopifnot(nrow(res$catt_df) == res$R)
	# Mutate: trim catt_df to 1 row.
	res$catt_df <- res$catt_df[1, , drop = FALSE]
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C4 catt_df: nrow == R",
		fixed = TRUE
	)
})

test_that("C5 cohort_probs: validator catches wrong length", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	stopifnot(length(res$cohort_probs) == res$R)
	# Mutate: drop one entry.
	res$cohort_probs <- res$cohort_probs[1]
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C5 cohort_probs length == R",
		fixed = TRUE
	)
})

test_that("C6 dimensions: validator catches beta_hat length mismatch", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	stopifnot(length(res$beta_hat) == res$p)
	# Mutate: trim beta_hat by 1.
	res$beta_hat <- res$beta_hat[-1]
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C6 length(beta_hat) == p",
		fixed = TRUE
	)
})

test_that("C7 lambda monotonicity: validator catches inverted model sizes", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	# Mutate: swap max and min model sizes to violate monotonicity.
	res$lambda.max_model_size <- res$lambda.min_model_size + 5L
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C7 lambda.max_model_size <= lambda_star_model_size",
		fixed = TRUE
	)
})

test_that("C8 type sanity: validator catches malformed se_type", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	stopifnot(res$se_type %in% c("default", "cluster"))
	res$se_type <- "standard" # not a valid choice
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"C8 se_type",
		fixed = TRUE
	)
})

# ------------------------------------------------------------------------------
# (c) Cross-class consistency: validator expected slots == live names()
# ------------------------------------------------------------------------------

test_that("validator expected-slots lists match live names() for all 4 classes", {
	fx <- .validator_fixture(q = 0.5)
	expect_setequal(names(fx$fetwfe), fetwfe:::.EXPECTED_SLOTS_FETWFE)
	expect_setequal(names(fx$etwfe), fetwfe:::.EXPECTED_SLOTS_ETWFE)
	expect_setequal(names(fx$betwfe), fetwfe:::.EXPECTED_SLOTS_BETWFE)
	expect_setequal(names(fx$twfeCovs), fetwfe:::.EXPECTED_SLOTS_TWFECOVS)
	# Also: fetwfe internal sub-slots
	expect_setequal(
		names(fx$fetwfe$internal),
		fetwfe:::.EXPECTED_INTERNAL_SLOTS_FETWFE
	)
})

# ------------------------------------------------------------------------------
# (d) Slot-deletion mutation: NULL-ing a key slot is caught
# ------------------------------------------------------------------------------

test_that("validator catches NULL key slots across all 4 classes", {
	fx <- .validator_fixture(q = 0.5)
	# Each class: NULL a different key slot.
	res_f <- fx$fetwfe
	res_f$catt_df <- NULL
	expect_error(
		fetwfe:::.validate_fetwfe(res_f),
		"Missing slot(s): catt_df",
		fixed = TRUE
	)

	res_e <- fx$etwfe
	res_e$cohort_probs_overall <- NULL
	expect_error(
		fetwfe:::.validate_etwfe(res_e),
		"Missing slot(s): cohort_probs_overall",
		fixed = TRUE
	)

	res_b <- fx$betwfe
	res_b$lambda_star <- NULL
	expect_error(
		fetwfe:::.validate_betwfe(res_b),
		"Missing slot(s): lambda_star",
		fixed = TRUE
	)

	res_t <- fx$twfeCovs
	res_t$se_type <- NULL
	expect_error(
		fetwfe:::.validate_twfeCovs(res_t),
		"Missing slot(s): se_type",
		fixed = TRUE
	)
})

test_that("validator catches missing internal sub-slot on fetwfe", {
	fx <- .validator_fixture(q = 0.5)
	res <- fx$fetwfe
	res$internal$theta_hat <- NULL
	expect_error(
		fetwfe:::.validate_fetwfe(res),
		"in `internal`",
		fixed = TRUE
	)
})
