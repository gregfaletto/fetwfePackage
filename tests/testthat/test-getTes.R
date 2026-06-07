library(testthat)
library(fetwfe)

# Helper function to compute the expected number of treatment effects.
compute_num_treats <- function(T, R) {
	T * R - (R * (R + 1)) / 2
}

# ------------------------------------------------------------------------------
# Test 1: Check that getTes returns a list with the expected output components.
# ------------------------------------------------------------------------------
test_that("getTes returns a list with att_true and actual_cohort_tes", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	# Call getTes() with the coefs object
	res <- genCoefs(
		G = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	) |>
		getTes()

	expect_type(res, "list")
	expect_true("att_true" %in% names(res))
	expect_true("actual_cohort_tes" %in% names(res))
})

# ------------------------------------------------------------------------------
# Test 2: Check that att_true equals the mean of actual_cohort_tes.
# ------------------------------------------------------------------------------
test_that("att_true equals mean(actual_cohort_tes)", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res <- genCoefs(
		G = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	) |>
		getTes()

	expect_equal(
		res$att_true,
		as.numeric(mean(res$actual_cohort_tes)),
		tolerance = 1e-6
	)
})

# ------------------------------------------------------------------------------
# Test 3: Check that actual_cohort_tes equals the expected cohort effects computed via getActualCohortTes.
# ------------------------------------------------------------------------------
test_that("actual_cohort_tes matches expected output", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res_coefs <- genCoefs(
		G = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	)
	beta <- res_coefs$beta

	res <- getTes(res_coefs)

	# Compute number of treatment effects and determine the indices for treatment dummy coefficients.
	num_treats <- compute_num_treats(T_val, R_val)
	base_cols <- if (d_val > 0) {
		R_val + (T_val - 1) + d_val + d_val * R_val + d_val * (T_val - 1)
	} else {
		R_val + (T_val - 1)
	}
	treat_inds <- seq(from = base_cols + 1, length.out = num_treats)

	# Use getFirstInds() to get the starting indices.
	first_inds <- getFirstInds(G = R_val, T = T_val)

	actual_cohort_tes_expected <- getActualCohortTes(
		G = R_val,
		first_inds = first_inds,
		treat_inds = treat_inds,
		coefs = beta,
		num_treats = num_treats
	)

	expect_equal(res$actual_cohort_tes, actual_cohort_tes_expected)
	expect_equal(
		res$att_true,
		as.numeric(mean(actual_cohort_tes_expected)),
		tolerance = 1e-6
	)
})

# ------------------------------------------------------------------------------
# Test 4: Check that getTes errors when beta has an incorrect length.
# ------------------------------------------------------------------------------
test_that("getTes errors with invalid beta length", {
	R_val <- 5
	T_val <- 30
	d_val <- 12

	res_coefs <- genCoefs(
		G = R_val,
		T = T_val,
		d = d_val,
		density = 0.1,
		eff_size = 2,
		seed = 234
	)

	# Create a faulty coefs object with incorrect beta length.
	faulty_coefs <- res_coefs
	faulty_coefs$beta <- faulty_coefs$beta[-1]

	expect_error(
		getTes(faulty_coefs),
		regexp = "length(beta) == p is not TRUE",
		fixed = TRUE
	)
})

# ------------------------------------------------------------------------------
# Test 5: getTes returns an object of class FETWFE_tes with carried params.
# ------------------------------------------------------------------------------
test_that("getTes returns an object of class FETWFE_tes with carried params", {
	coefs <- genCoefs(
		G = 5,
		T = 30,
		d = 12,
		density = 0.1,
		eff_size = 2,
		seed = 234
	)
	res <- getTes(coefs)
	expect_s3_class(res, "FETWFE_tes")
	expect_equal(res$R, 5)
	expect_equal(res$T, 30)
	expect_equal(res$d, 12)
	expect_equal(res$seed, 234)
})

# ------------------------------------------------------------------------------
# Test 6: print method writes the expected lines.
# NOTE: the regexes below intentionally lock the print-output labels
# ("Cohorts (G)", "Time periods (T)", "Covariates (d)"). Any cosmetic
# relabeling of print.FETWFE_tes requires updating this test.
# ------------------------------------------------------------------------------
test_that("print.FETWFE_tes writes the expected sections", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	)
	out <- capture.output(print(getTes(coefs)))
	joined <- paste(out, collapse = "\n")
	expect_match(joined, "Overall true ATT")
	# Cohorts are labeled by calendar adoption time (g + 1), so for G = 3 the
	# labels are Cohort 2/3/4, not the 1-based index (#261).
	expect_match(joined, "Cohort 2")
	expect_match(joined, "Cohort 3")
	expect_match(joined, "Cohort 4")
	expect_match(joined, "Cohorts \\(G\\)")
	expect_match(joined, "Time periods \\(T\\)")
	expect_match(joined, "Covariates \\(d\\)")
	expect_match(joined, "Seed")
})

# ------------------------------------------------------------------------------
# #261: print / summary / tidy.FETWFE_tes must label cohorts identically, by
# calendar adoption time (cohort g adopts at time g + 1), matching the fitted
# estimators' catt_df$cohort -- never the 1-based loop index.
# ------------------------------------------------------------------------------
test_that("FETWFE_tes cohort labels agree across print/summary/tidy (adoption time, #261)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	)
	tes <- getTes(coefs)
	expect_equal(tes$cohort_times, as.integer(c(2, 3, 4)))

	p <- paste(capture.output(print(tes)), collapse = "\n")
	s <- paste(capture.output(print(summary(tes))), collapse = "\n")
	td <- broom::tidy(tes)

	# All three label cohorts by adoption time (2, 3, 4).
	for (lab in c("Cohort 2:", "Cohort 3:", "Cohort 4:")) {
		expect_match(p, lab, fixed = TRUE)
		expect_match(s, lab, fixed = TRUE)
	}
	expect_setequal(td$term, c("ATT_true", "Cohort 2", "Cohort 3", "Cohort 4"))

	# The old 1-based index labeling is gone (cohort 1 is now "Cohort 2").
	expect_no_match(p, "Cohort 1:", fixed = TRUE)
	expect_no_match(s, "Cohort 1:", fixed = TRUE)
	expect_false("Cohort 1" %in% td$term)
})

# ------------------------------------------------------------------------------
# Test 7: summary method returns the documented fields and dispersion stats.
# ------------------------------------------------------------------------------
test_that("summary.FETWFE_tes returns expected fields and dispersion stats", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	)
	s <- summary(getTes(coefs))
	expect_s3_class(s, "summary.FETWFE_tes")
	expect_named(s$cohort_te_stats, c("min", "max", "median", "sd"))
	expect_equal(
		unname(s$cohort_te_stats["min"]),
		min(s$actual_cohort_tes)
	)
	expect_equal(
		unname(s$cohort_te_stats["max"]),
		max(s$actual_cohort_tes)
	)
	expect_equal(
		unname(s$cohort_te_stats["median"]),
		stats::median(s$actual_cohort_tes)
	)
	expect_equal(
		unname(s$cohort_te_stats["sd"]),
		stats::sd(s$actual_cohort_tes)
	)

	out <- capture.output(print(s))
	joined <- paste(out, collapse = "\n")
	expect_match(joined, "Cohort effect dispersion")
})

# ------------------------------------------------------------------------------
# Test 8 (#189): the "Cohort assignment DGP" section renders only for
# non-marginal objects, and getTes() carries assignment_type /
# assignment_strength so print()/summary() can describe the DGP.
# ------------------------------------------------------------------------------
test_that("print/summary.FETWFE_tes show the assignment-DGP section only when non-marginal", {
	# Marginal: slots present, but the new section is NOT rendered, so the
	# common-case output is unchanged.
	tm <- getTes(genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 1
	))
	expect_identical(tm$assignment_type, "marginal")
	pm <- paste(capture.output(print(tm)), collapse = "\n")
	sm <- paste(capture.output(print(summary(tm))), collapse = "\n")
	expect_no_match(pm, "Cohort assignment DGP", fixed = TRUE)
	expect_no_match(sm, "Cohort assignment DGP", fixed = TRUE)

	# Non-marginal: getTes() carries the slots, and both print() and summary()
	# render the type, strength, and cohort weights.
	cn <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		assignment_type = "multinomial",
		assignment_strength = 1.0,
		seed = 42
	)
	tn <- getTes(cn)
	expect_identical(tn$assignment_type, "multinomial")
	expect_equal(tn$assignment_strength, 1.0)
	expect_length(tn$cohort_weights, 3L)

	pn <- paste(capture.output(print(tn)), collapse = "\n")
	expect_match(pn, "Cohort assignment DGP", fixed = TRUE)
	expect_match(pn, "multinomial", fixed = TRUE)
	expect_match(pn, "Cohort weights", fixed = TRUE)

	sn <- paste(capture.output(print(summary(tn))), collapse = "\n")
	expect_match(sn, "Cohort assignment DGP", fixed = TRUE)
	expect_match(sn, "multinomial", fixed = TRUE)
})
