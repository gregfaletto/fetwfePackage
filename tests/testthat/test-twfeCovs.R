library(testthat)
library(fetwfe)

# Panel-data fixture helpers (generate_panel_data, generate_bad_panel_data,
# generate_minimal_panel_data) are defined in tests/testthat/helper-panel-fixture.R
# and sourced by testthat before this file runs (issue #91).

# ------------------------------------------------------------------------------
# Test 1: Check that valid input produces a list with the expected output elements.
# ------------------------------------------------------------------------------
test_that("twfeCovs returns expected output structure with valid input", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	expect_type(result, "list")

	# Check that the returned list contains the main expected names.
	expected_names <- c(
		"att_hat",
		"att_se",
		"catt_hats",
		"catt_ses",
		"catt_df",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"X_ints",
		"y",
		"X_final",
		"y_final",
		"N",
		"T",
		"R",
		"d",
		"p"
	)
	for (nm in expected_names) {
		expect_true(nm %in% names(result))
	}

	# Also check that catt_df is a data frame with the expected column names.
	expect_s3_class(result$catt_df, "data.frame")
	expect_true(all(
		c("cohort", "estimate", "se", "ci_low", "ci_high") %in%
			colnames(result$catt_df)
	))
})

# ------------------------------------------------------------------------------
# Test 2: Error when pdata is not a data.frame
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when pdata is not a data.frame", {
	expect_error(
		twfeCovs(
			pdata = "not a data.frame",
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"pdata must be a data frame"
	)
})

# ------------------------------------------------------------------------------
# Test 3: Error when time_var column is not integer.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when time_var is not integer", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$time <- as.character(df$time) # convert time to character

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"time_var column 'time' must be integer"
	)
})

# ------------------------------------------------------------------------------
# Test 4: Error when unit_var column is not character.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when unit_var is not character", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$unit <- as.factor(df$unit) # make unit a factor

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"unit_var column 'unit' must be character"
	)
})

# ------------------------------------------------------------------------------
# Test 5: Error when treatment column is not integer.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when treatment is not integer", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$treatment <- as.character(df$treatment) # wrong type

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"treatment column 'treatment' must be integer"
	)
})

# ------------------------------------------------------------------------------
# Test 6: Error when treatment column contains values other than 0 and 1.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when treatment has values other than 0/1", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	# Force an invalid value into treatment (e.g. 2)
	df$treatment[1] <- 2L

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"treatment column 'treatment' must contain only 0/1"
	)
})

# ------------------------------------------------------------------------------
# Test 7: Error when a specified covariate column is missing.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when a covariate column is missing", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	# Remove one of the covariate columns
	df <- df[, !(names(df) %in% "cov2")]

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"covs contains name\\(s\\) not in pdata"
	)
})

# ------------------------------------------------------------------------------
# Test 8: Error when response column is not numeric.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when response is not numeric", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$y <- as.character(df$y)

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"response column 'y' must be numeric or integer"
	)
})

# ------------------------------------------------------------------------------
# Test 9: Warning when alpha > 0.5.
# ------------------------------------------------------------------------------
test_that("twfeCovs warns when alpha > 0.5", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	expect_warning(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			alpha = 0.6,
			verbose = FALSE
		),
		"Provided alpha > 0.5"
	)
})

# ------------------------------------------------------------------------------
# Test 10: Error when all units are treated in the first period.
# Create a panel in which every unit is treated in the first period.
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when all units are treated in the first period", {
	# Create a panel data frame with 5 units and 5 periods,
	# where treatment is 1 for all observations.
	df <- data.frame(
		time = rep(1:5, times = 5),
		unit = rep(sprintf("unit%02d", 1:5), each = 5),
		treatment = rep(1L, 25),
		cov1 = rnorm(25),
		cov2 = runif(25),
		y = rnorm(25)
	)

	expect_error(
		suppressWarnings(twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		)),
		"All units were treated in the first time period or no units remain after filtering; estimating treatment effects is not possible"
	)
})

# ------------------------------------------------------------------------------
# Test 11: Error when there are no never-treated units.
# ------------------------------------------------------------------------------

test_that("twfeCovs errors with a panel having no never-treated units when allow_no_never_treated = FALSE", {
	df_bad <- generate_bad_panel_data(N = 30, T = 10, seed = 123)

	expect_error(
		suppressWarnings(twfeCovs(
			pdata = df_bad,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE,
			allow_no_never_treated = FALSE
		)),
		"allow_no_never_treated"
	)
})

# ------------------------------------------------------------------------------
# Test 12: Error when the data has too few rows (less than 4).
# ------------------------------------------------------------------------------
test_that("twfeCovs errors when data has fewer than 4 rows", {
	df <- data.frame(
		time = as.integer(1:3),
		unit = as.character(c("u1", "u2", "u3")),
		treatment = as.integer(c(0, 0, 0)),
		cov1 = rnorm(3),
		cov2 = runif(3),
		y = rnorm(3)
	)

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y"
		),
		"pdata must have at least 4 rows"
	)
})

# ------------------------------------------------------------------------------
# Test 13: Test that the function returns a standard error
# ------------------------------------------------------------------------------
test_that("twfeCovs returns att_se", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	result1 <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	expect_true(result1$calc_ses)
	expect_false(is.na(result1$att_se))
})

# ------------------------------------------------------------------------------
# Test 14: Test that a warning is thrown when all covariates are constant.
# ------------------------------------------------------------------------------
test_that("twfeCovs warns when all covariates are removed due to constant values", {
	df_const <- generate_panel_data(N = 30, T = 5, R = 2, seed = 101)
	# Set both covariates to a constant value
	df_const$cov1 <- 1
	df_const$cov2 <- 1

	expect_warning(
		twfeCovs(
			pdata = df_const,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		),
		regexp = "All remaining covariates were removed because they were constant across units",
		fixed = TRUE, # literal match
		all = FALSE # accept extra warnings
	)
})

# ------------------------------------------------------------------------------
# Test 15: Test that the overall ATT (att_hat) is numeric and non-missing.
# ------------------------------------------------------------------------------
test_that("twfeCovs returns a valid numeric overall ATT", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 202)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})


# ------------------------------------------------------------------------------
# Test 16: Test that function works on a minimal data set
# ------------------------------------------------------------------------------
test_that("twfeCovs works on a minimal valid dataset", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df_min <- generate_minimal_panel_data(include_cov1 = FALSE)

	result <- twfeCovs(
		pdata = df_min,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		# covs = c("cov1"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1,
	)

	expect_type(result$att_hat, "double")
	expect_false(is.na(result$att_hat))
})

# ------------------------------------------------------------------------------
# Test 17: Test that minimal data set requires provided noise variance
# ------------------------------------------------------------------------------
test_that("minimal data set requires provided noise variance", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df_min <- generate_minimal_panel_data(include_cov1 = FALSE)

	expect_error(
		twfeCovs(
			pdata = df_min,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			# covs = c("cov1"),
			response = "y",
			verbose = FALSE
		),
		"Not enough units available to estimate the noise variance."
	)
})

# ------------------------------------------------------------------------------
# Test 18: Test that function works on only two cohorts
# ------------------------------------------------------------------------------
test_that("twfeCovs works on only two cohorts", {
	# Create a balanced panel with N = 6 units and T = 5 time periods.
	df <- generate_panel_data(N = 30, T = 10, R = 2, seed = 202)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	expect_type(result$att_hat, "double")
	expect_false(is.na(result$att_hat))
})

# ------------------------------------------------------------------------------
# Test 19: Test that at least 2 treated cohorts requires
# ------------------------------------------------------------------------------
test_that("at least two treated cohorts required", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df <- generate_panel_data(N = 30, T = 10, R = 1, seed = 202)

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		),
		"Only one treated cohort detected in data. Currently fetwfe and etwfe only support data sets with at least two treated cohorts."
	)
})


# ------------------------------------------------------------------------------
# Test 20: Overall ATT standard error is computed (non‐zero)
# ------------------------------------------------------------------------------
test_that("Overall ATT standard error is non-negative", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 303)

	# q < 1: expect an overall standard error to be computed
	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1, # provided noise variances
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$att_se))
	expect_true(result$att_se >= 0)
})

# ------------------------------------------------------------------------------
# Test 21: Cohort-specific treatment effect standard errors are computed,
# nonnegative, and match the number of treated cohorts.
# ------------------------------------------------------------------------------
test_that("Cohort-specific standard errors are computed correctly", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 505)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$catt_ses))
	expect_equal(length(result$catt_ses), result$R)
	expect_true(all(result$catt_ses >= 0))
})

# ------------------------------------------------------------------------------
# Test 22: The estimator works when no covariates are provided.
# (processCovs() should issue a warning but continue.)
# ------------------------------------------------------------------------------
test_that("Estimator works with no covariates", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 606)

	# Call twfeCovs with no covariate argument.
	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		verbose = TRUE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
	expect_true(is.numeric(result$att_se))
	expect_true(result$att_se >= 0)
	# Also, the returned 'd' (number of covariates) should be 0.
	expect_equal(result$d, 0)
})

# ------------------------------------------------------------------------------
# Test 23: processCovs() properly handles an empty covariate vector.
# ------------------------------------------------------------------------------
test_that("processCovs handles empty covariate vector correctly", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 707)

	res <- processCovs(
		df = df,
		units = unique(df$unit),
		unit_var = "unit",
		times = sort(unique(df$time)),
		time_var = "time",
		covs = character(0),
		resp_var = "y",
		T = 5,
		verbose = TRUE
	)

	expect_true(is.data.frame(res$df))
	expect_equal(length(res$covs), 0)
})

# ------------------------------------------------------------------------------
# Test 24: Overall ATT and cohort-specific estimates are finite and numeric.
# ------------------------------------------------------------------------------
test_that("Overall and cohort-specific treatment effects are valid", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 808)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
	expect_true(is.numeric(result$catt_hats))
	expect_equal(length(result$catt_hats), result$R)
	expect_true(all(!is.na(result$catt_hats)))
})

# ------------------------------------------------------------------------------
# Test 25: Second test when no covariates are provided.
# (processCovs() should issue a warning but continue.)
# ------------------------------------------------------------------------------
test_that("Estimator works with no covariates again", {
	set.seed(123456L)

	# 5 time periods, 30 individuals, and 4 waves of treatment
	tmax = 5
	imax = 30
	nlvls = 5

	dat =
		expand.grid(time = 1:tmax, id = 1:imax) |>
		within({
			# Initialize columns
			cohort = rep(1:nlvls, each = imax / nlvls)[id]
			effect = NA
			first_treat = NA

			for (lvls in 1:nlvls) {
				effect = ifelse(cohort == lvls, sample(2:10, 1), effect)
				first_treat = ifelse(cohort == lvls, lvls + 1, first_treat)
			}

			first_treat = ifelse(first_treat > tmax, Inf, first_treat)
			treat = time >= first_treat
			rel_time = time - first_treat
			y = id +
				time +
				ifelse(treat, effect * rel_time, 0) +
				rnorm(imax * tmax)

			rm(lvls, cohort, effect)
		})

	# Convert `dat` to `pdata` based on the specified requirements

	# Specify column names for the pdata format
	time_var <- "time" # Column for the time period
	unit_var <- "unit" # Column for the unit identifier
	treatment <- "treated" # Column for the treatment dummy indicator
	covs <- c() # Columns for covariates
	response <- "response" # Column for the response variable

	# Convert the dataset
	pdata <- dat |>
		dplyr::mutate(
			# Rename id to unit and convert to character
			{{ unit_var }} := as.character(id),
			# Ensure treatment dummy is 0/1
			{{ treatment }} := as.integer(treat),
			# Rename y to response
			{{ response }} := y
		) |>
		dplyr::select(
			{{ time_var }},
			{{ unit_var }},
			{{ treatment }},
			{{ response }}
		)

	result <- twfeCovs(
		pdata = pdata, # The panel dataset
		time_var = "time", # The time variable
		unit_var = "unit", # The unit identifier
		treatment = "treated", # The treatment dummy indicator
		response = "response", # The response variable
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
	expect_true(result$att_hat != 0)
	expect_true(is.numeric(result$att_se))
	expect_true(result$att_se > 0)
	# Also, the returned 'd' (number of covariates) should be 0.
	expect_equal(result$d, 0)
})

# ------------------------------------------------------------------------------
# Test 26: Test that adding ridge regularization works
# ------------------------------------------------------------------------------
test_that("adding ridge regularization to twfeCovs works", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 202)

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		add_ridge = TRUE
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		add_ridge = TRUE
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})

test_that("data application works with ridge penalty", {
	set.seed(23451)

	library(bacondecomp)

	data(divorce)

	# sig_eps_sq and sig_eps_c_sq, calculated in a separate run of `etwfe(),
	# are provided to speed up the computation of the example
	res <- suppressWarnings(twfeCovs(
		pdata = divorce[divorce$sex == 2, ],
		time_var = "year",
		unit_var = "st",
		treatment = "changed",
		covs = c("murderrate", "lnpersinc", "afdcrolls"),
		response = "suiciderate_elast_jag",
		sig_eps_sq = 0.0344,
		sig_eps_c_sq = 0.1507,
		verbose = TRUE,
		add_ridge = TRUE
	))

	expect_true(is.numeric(res$att_hat))
	expect_false(is.na(res$att_hat))
})

# ------------------------------------------------------------------------------
# Test 27: Test that processFactors() converts factor covariates to dummies
# ------------------------------------------------------------------------------
test_that("processFactors processes factor covariates correctly", {
	set.seed(987)
	df <- data.frame(
		cov1 = rnorm(10),
		cov2 = factor(sample(c("A", "B", "C"), 10, replace = TRUE))
	)
	res <- processFactors(df, c("cov1", "cov2"))

	# The original factor "cov2" should be removed from pdata.
	expect_true(!("cov2" %in% colnames(res$pdata)))

	# The new covariate names (in res$covs) should not include "cov2" but should include dummy names.
	expect_true(!("cov2" %in% res$covs))
	expect_true(any(grepl("cov2_", res$covs)))

	# Also, check that the number of covariates is: 1 (for cov1) plus (nlevels(cov2)-1).
	expected_dummies <- nlevels(df$cov2) - 1
	expect_equal(length(res$covs), 1 + expected_dummies)
})

# ------------------------------------------------------------------------------
# Test 28: Test that twfeCovs() handles factor covariates appropriately
# ------------------------------------------------------------------------------
test_that("twfeCovs handles factor covariates appropriately", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	# Convert cov2 into a factor with 3 levels.
	set.seed(123)
	df$cov2 <- factor(sample(c("A", "B", "C"), nrow(df), replace = TRUE))

	result <- twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	# Since cov1 is numeric (1 column) and cov2 is a factor with 3 levels,
	# processFactors() should convert cov2 into 2 dummy columns.
	# Thus, the total number of covariates 'd' should be 1 + 2 = 3.
	expect_equal(result$d, 3)

	# (Optional) Run processFactors() separately and verify that the new dummy names are returned.
	pf_res <- processFactors(df, c("cov1", "cov2"))
	expect_true(!("cov2" %in% pf_res$covs))
	expect_true(any(grepl("cov2_", pf_res$covs)))
})

# ------------------------------------------------------------------------------
# Test 29: Check that tibbles work
# ------------------------------------------------------------------------------
test_that("tibbles work as input to twfeCovs", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123) # uses N = 30, T = 10 by default

	result <- twfeCovs(
		pdata = tibble::as_tibble(df),
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)

	expect_type(result, "list")

	# Check that the returned list contains the main expected names.
	expected_names <- c(
		"att_hat",
		"att_se",
		"catt_hats",
		"catt_ses",
		"catt_df",
		"beta_hat",
		"treat_inds",
		"treat_int_inds",
		"sig_eps_sq",
		"sig_eps_c_sq",
		"X_ints",
		"y",
		"X_final",
		"y_final",
		"N",
		"T",
		"R",
		"d",
		"p"
	)
	for (nm in expected_names) {
		expect_true(nm %in% names(result))
	}

	# Also check that catt_df is a data frame with the expected column names.
	expect_s3_class(result$catt_df, "data.frame")
	expect_true(all(
		c("cohort", "estimate", "se", "ci_low", "ci_high") %in%
			colnames(result$catt_df)
	))
})

# ------------------------------------------------------------------------------
# Test 30: Error when a cohort contains fewer than d + 1 units
# ------------------------------------------------------------------------------
test_that("twfeCovs throws error when a cohort contains fewer than d + 1 units", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 123)

	expect_error(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		),
		"At least one cohort contains fewer than d \\+ 1 units\\. The design matrix is rank-deficient\\. Calculating standard errors will not be possible, and estimating treatment effects is only possible using add_ridge = TRUE\\."
	)

	expect_warning(
		twfeCovs(
			pdata = df,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE,
			add_ridge = TRUE
		),
		"At least one cohort contains fewer than d \\+ 1 units\\. The design matrix is rank-deficient\\. Calculating standard errors will not be possible, and estimating treatment effects is only possible using add_ridge = TRUE\\."
	)

	res <- suppressWarnings(twfeCovs(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		add_ridge = TRUE
	))

	expect_true(!is.na(res$att_hat))
	expect_true(res$att_hat != 0)
})

# ------------------------------------------------------------------------------
# Test: twfeCovs surfaces p_value but not selected in catt_df
# ------------------------------------------------------------------------------
test_that("twfeCovs surfaces p_value but not selected in catt_df", {
	set.seed(2026)
	sim <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(
		sim,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- twfeCovsWithSimulatedData(dat, verbose = FALSE)

	expect_true("p_value" %in% colnames(res$catt_df))
	expect_false("selected" %in% colnames(res$catt_df))

	expect_true("att_p_value" %in% names(res))
	expect_false("att_selected" %in% names(res))
	expect_length(res$att_p_value, 1)

	non_na <- !is.na(res$catt_df$p_value)
	expect_true(all(res$catt_df$p_value[non_na] >= 0))
	expect_true(all(res$catt_df$p_value[non_na] <= 1))
})

# ------------------------------------------------------------------------------
# Test: twfeCovs auto-truncates a no-never-treated panel (default)
# ------------------------------------------------------------------------------
test_that("twfeCovs auto-truncates a panel with no never-treated units (default)", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)

	expect_warning(
		res <- twfeCovs(
			pdata = df_bad,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		),
		"auto-truncated"
	)
	expect_true(is.list(res))
	expect_true(is.finite(res$att_hat))
	expect_s3_class(res$catt_df, "data.frame")
})

# ------------------------------------------------------------------------------
# Issue #116 Gap 2: se_type = "cluster" and the auto-truncation of an
# all-treated panel (allow_no_never_treated) are each tested in isolation, but
# never together. This exercises the interaction end-to-end: a cluster SE
# computed on an auto-truncated panel must still be finite and strictly
# positive. The > 0 assertion is load-bearing -- on an all-treated panel the
# bridge estimators select every coefficient out and return att_se = 0
# exactly, which a finiteness-only check would pass vacuously. Parallel to the
# etwfe Gap 2 test in test-etwfe.R.
# ------------------------------------------------------------------------------
test_that("twfeCovs cluster SE is finite and positive on an auto-truncated panel", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)

	expect_warning(
		res <- twfeCovs(
			pdata = df_bad,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			se_type = "cluster",
			allow_no_never_treated = TRUE,
			verbose = FALSE
		),
		"auto-truncated"
	)
	expect_s3_class(res, "twfeCovs")
	expect_true(is.finite(res$att_se))
	expect_false(is.na(res$att_se))
	expect_gt(res$att_se, 0)
})

# ------------------------------------------------------------------------------
# Test: twfeCovs errors cleanly when truncation would yield < 2 retained cohorts
# ------------------------------------------------------------------------------
test_that("twfeCovs errors cleanly when no-never-treated truncation would yield < 2 cohorts", {
	N <- 5
	T_periods <- 6
	r_single <- 4
	df_single <- data.frame(
		time = rep(seq_len(T_periods), times = N),
		unit = rep(paste0("u", seq_len(N)), each = T_periods),
		treatment = rep(
			c(rep(0L, r_single - 1L), rep(1L, T_periods - r_single + 1L)),
			times = N
		),
		cov1 = rep(rnorm(N, mean = 0), each = T_periods),
		cov2 = rep(rnorm(N, mean = 0), each = T_periods),
		y = rnorm(N * T_periods)
	)

	expect_error(
		suppressWarnings(twfeCovs(
			pdata = df_single,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1", "cov2"),
			response = "y",
			verbose = FALSE
		)),
		"treated cohort"
	)
})

# ------------------------------------------------------------------------------
# Helpers for se_type tests: small F1-conforming and AR(1)-corrupted panels.
# simulateData() produces F1-conforming shocks only; for the serially-correlated
# fixture we add an AR(1) shock per unit on top of the simulated response.
# ------------------------------------------------------------------------------
make_se_type_panel <- function(seed = 2026, N = 120) {
	set.seed(seed)
	sim_coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	simulateData(sim_coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

make_ar1_panel <- function(seed = 7, N = 150, rho = 0.85, sd_e = 1) {
	set.seed(seed)
	sim_coefs <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	sim <- simulateData(
		sim_coefs,
		N = N,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5
	)
	pdata <- sim$pdata
	pdata <- pdata[order(pdata$unit, pdata$time), ]
	units <- unique(pdata$unit)
	T_ <- length(unique(pdata$time))
	set.seed(seed + 100)
	for (u in units) {
		idx <- which(pdata$unit == u)
		ar <- numeric(T_)
		ar[1] <- rnorm(1, sd = sd_e / sqrt(1 - rho^2))
		for (t in 2:T_) {
			ar[t] <- rho * ar[t - 1] + rnorm(1, sd = sd_e)
		}
		pdata[idx, sim$response] <- pdata[idx, sim$response] + ar
	}
	sim$pdata <- pdata
	sim
}

# ------------------------------------------------------------------------------
# Test: twfeCovs with se_type omitted matches se_type = "default" identically
# ------------------------------------------------------------------------------
test_that("twfeCovs se_type omitted matches se_type = 'default' identically", {
	sim <- make_se_type_panel()
	res_omit <- twfeCovsWithSimulatedData(sim)
	res_def <- twfeCovsWithSimulatedData(sim, se_type = "default")
	expect_identical(res_omit$att_se, res_def$att_se)
	expect_identical(res_omit$catt_ses, res_def$catt_ses)
	expect_identical(res_omit$att_hat, res_def$att_hat)
	expect_identical(res_omit$catt_hats, res_def$catt_hats)
})

# ------------------------------------------------------------------------------
# Test: twfeCovs with se_type = "cluster" returns finite SEs on a clean panel
# ------------------------------------------------------------------------------
test_that("twfeCovs se_type = 'cluster' returns finite SEs on a clean panel", {
	sim <- make_se_type_panel()
	res <- twfeCovsWithSimulatedData(sim, se_type = "cluster")
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	expect_true(all(is.finite(res$catt_ses)))
	expect_true(all(res$catt_ses > 0))
})

# ------------------------------------------------------------------------------
# Test: under deliberately serially-correlated DGP, cluster SE > default SE
# ------------------------------------------------------------------------------
test_that("twfeCovs cluster-robust SE exceeds default SE under AR(1) shocks", {
	mk <- make_ar1_panel()
	res_def <- twfeCovsWithSimulatedData(mk, se_type = "default")
	res_cls <- twfeCovsWithSimulatedData(mk, se_type = "cluster")
	expect_true(is.finite(res_def$att_se))
	expect_true(is.finite(res_cls$att_se))
	expect_gt(res_def$att_se, 0)
	expect_gt(res_cls$att_se, res_def$att_se)
})

# ------------------------------------------------------------------------------
# Test: $se_type slot on output reflects the argument value
# ------------------------------------------------------------------------------
test_that("twfeCovs $se_type slot reflects the argument value", {
	sim <- make_se_type_panel()
	res_def <- twfeCovsWithSimulatedData(sim)
	res_cls <- twfeCovsWithSimulatedData(sim, se_type = "cluster")
	expect_identical(res_def$se_type, "default")
	expect_identical(res_cls$se_type, "cluster")
})

# ------------------------------------------------------------------------------
# Test (#76 Item 1): twfeCovs() returns a classed object with print + coef methods.
# Prior to this PR the return was an unclassed named list; `print()` fell
# through to `print.default` and dumped the full design matrix.
# ------------------------------------------------------------------------------
test_that("twfeCovs returns a classed object with working print + coef methods (#76 Item 1)", {
	set.seed(42)
	sim <- simulateData(
		genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 42),
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	fit <- twfeCovsWithSimulatedData(sim)
	expect_s3_class(fit, "twfeCovs")
	expect_no_error(capture.output(print(fit)))
	expect_equal(coef(fit), fit$beta_hat)
})
