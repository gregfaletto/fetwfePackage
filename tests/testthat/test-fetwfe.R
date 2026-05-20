library(testthat)
library(fetwfe)

# Panel-data fixture helpers (generate_panel_data, generate_bad_panel_data,
# generate_minimal_panel_data) are defined in tests/testthat/helper-panel-fixture.R
# and sourced by testthat before this file runs (issue #91).

# ------------------------------------------------------------------------------
# Test 1: Check that valid input produces a list with the expected output elements.
# ------------------------------------------------------------------------------
test_that("fetwfe returns expected output structure with valid input", {
	df <- generate_panel_data() # uses N = 30, T = 10 by default

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
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
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"internal",
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
		c("Cohort", "Estimated TE", "SE", "ConfIntLow", "ConfIntHigh") %in%
			colnames(result$catt_df)
	))
})

# ------------------------------------------------------------------------------
# Test 2: Error when pdata is not a data.frame
# ------------------------------------------------------------------------------
test_that("fetwfe errors when pdata is not a data.frame", {
	expect_error(
		fetwfe(
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
test_that("fetwfe errors when time_var is not integer", {
	df <- generate_panel_data()
	df$time <- as.character(df$time) # convert time to character

	expect_error(
		fetwfe(
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
test_that("fetwfe errors when unit_var is not character", {
	df <- generate_panel_data()
	df$unit <- as.factor(df$unit) # make unit a factor

	expect_error(
		fetwfe(
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
test_that("fetwfe errors when treatment is not integer", {
	df <- generate_panel_data()
	df$treatment <- as.character(df$treatment) # wrong type

	expect_error(
		fetwfe(
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
test_that("fetwfe errors when treatment has values other than 0/1", {
	df <- generate_panel_data()
	# Force an invalid value into treatment (e.g. 2)
	df$treatment[1] <- 2L

	expect_error(
		fetwfe(
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
test_that("fetwfe errors when a covariate column is missing", {
	df <- generate_panel_data()
	# Remove one of the covariate columns
	df <- df[, !(names(df) %in% "cov2")]

	expect_error(
		fetwfe(
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
test_that("fetwfe errors when response is not numeric", {
	df <- generate_panel_data()
	df$y <- as.character(df$y)

	expect_error(
		fetwfe(
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
test_that("fetwfe warns when alpha > 0.5", {
	df <- generate_panel_data()

	expect_warning(
		fetwfe(
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
test_that("fetwfe errors when all units are treated in the first period", {
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
		suppressWarnings(fetwfe(
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

test_that("fetwfe errors with a panel having no never-treated units when allow_no_never_treated = FALSE", {
	df_bad <- generate_bad_panel_data(N = 30, T = 10, seed = 123)

	expect_error(
		suppressWarnings(fetwfe(
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
test_that("fetwfe errors when data has fewer than 4 rows", {
	df <- data.frame(
		time = as.integer(1:3),
		unit = as.character(c("u1", "u2", "u3")),
		treatment = as.integer(c(0, 0, 0)),
		cov1 = rnorm(3),
		cov2 = runif(3),
		y = rnorm(3)
	)

	expect_error(
		fetwfe(
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
# Test 13: Test that optional lambda parameters are handled.
# ------------------------------------------------------------------------------
test_that("fetwfe returns expected output when lambda parameters are supplied", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 456)

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		lambda.max = 5,
		lambda.min = 0.5,
		nlambda = 50,
		verbose = FALSE
	)

	expect_true(is.numeric(result$lambda.max))
	expect_true(is.numeric(result$lambda.min))
	# Also check that some model size was recorded.
	expect_true(result$lambda.max_model_size >= 0)
})

# ------------------------------------------------------------------------------
# Test 14: Test that the function returns a standard error when q < 1 and not
# when q >= 1.
# ------------------------------------------------------------------------------
test_that("fetwfe returns att_se for q < 1 and att_se is NA for q >= 1", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 789)

	result1 <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
		verbose = FALSE
	)
	expect_true(result1$internal$calc_ses)
	expect_false(is.na(result1$att_se))

	result2 <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 2, # ridge-like penalty, where SE is not computed
		verbose = FALSE
	)
	expect_false(result2$internal$calc_ses)
	expect_true(is.na(result2$att_se))
	# Issue #84 item 3: q >= 1 ridge contract surface — `att_se` is NA AND
	# every dependent SE / p-value / CI is NA. Tightening the contract here
	# locks the surface against a future regression that surfaces a stale
	# bridge-regime SE through any of these slots.
	expect_true(is.na(result2$att_p_value))
	expect_true(all(is.na(result2$catt_ses)))
	expect_true(all(is.na(result2$catt_df$SE)))
	expect_true(all(is.na(result2$catt_df$ConfIntLow)))
	expect_true(all(is.na(result2$catt_df$ConfIntHigh)))
	expect_true(all(is.na(result2$catt_df$P_value)))
})

# ------------------------------------------------------------------------------
# Test 15: Test that a warning is thrown when all covariates are constant.
# ------------------------------------------------------------------------------
test_that("fetwfe warns when all covariates are removed due to constant values", {
	df_const <- generate_panel_data(N = 30, T = 10, R = 9, seed = 101)
	# Set both covariates to a constant value
	df_const$cov1 <- 1
	df_const$cov2 <- 1

	expect_warning(
		fetwfe(
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
# Test 16: Test that the overall ATT (att_hat) is numeric and non-missing.
# ------------------------------------------------------------------------------
test_that("fetwfe returns a valid numeric overall ATT", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 202)

	result <- fetwfe(
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
# Test 17: Test that function works on a minimal data set
# ------------------------------------------------------------------------------
test_that("fetwfe works on a minimal valid dataset", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df_min <- generate_minimal_panel_data()

	result <- fetwfe(
		pdata = df_min,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1"),
		response = "y",
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1,
	)

	expect_type(result$att_hat, "double")
	expect_false(is.na(result$att_hat))
})

# ------------------------------------------------------------------------------
# Test 18: Test that minimal data set requires provided noise variance
# ------------------------------------------------------------------------------
test_that("minimal data set requires provided noise variance", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df_min <- generate_minimal_panel_data()

	expect_error(
		fetwfe(
			pdata = df_min,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			covs = c("cov1"),
			response = "y",
			verbose = FALSE
		),
		"Not enough units available to estimate the noise variance."
	)
})

# ------------------------------------------------------------------------------
# Test 19: Test that function works on only two cohorts
# ------------------------------------------------------------------------------
test_that("fetwfe works on only two cohorts", {
	# Create a balanced panel with N = 6 units and T = 5 time periods.
	df <- generate_panel_data(N = 30, T = 10, R = 2, seed = 202)

	result <- fetwfe(
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
# Test 20: Test that at least 2 treated cohorts requires
# ------------------------------------------------------------------------------
test_that("at least two treated cohorts required", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df <- generate_panel_data(N = 30, T = 10, R = 1, seed = 202)

	expect_error(
		fetwfe(
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
# Test 21: Overall ATT standard error is computed (non‐zero) when q < 1,
# and is NA when q >= 1.
# ------------------------------------------------------------------------------
test_that("Overall ATT standard error is non-negative for q < 1", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 303)

	# q < 1: expect an overall standard error to be computed
	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1, # provided noise variances
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$att_se))
	expect_true(result$att_se >= 0)
})

test_that("Overall ATT standard error is NA for q >= 1", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 404)

	# q >= 1 (ridge-like): standard error is not computed
	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 2,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	expect_true(is.na(result$att_se))
	# Issue #84 item 3: q >= 1 ridge contract surface — same tightening as
	# the bridge-regime test in test 14 above.
	expect_true(is.na(result$att_p_value))
	expect_true(all(is.na(result$catt_ses)))
	expect_true(all(is.na(result$catt_df$SE)))
	expect_true(all(is.na(result$catt_df$ConfIntLow)))
	expect_true(all(is.na(result$catt_df$ConfIntHigh)))
	expect_true(all(is.na(result$catt_df$P_value)))
})

# ------------------------------------------------------------------------------
# Test 22: Cohort-specific treatment effect standard errors are computed,
# nonnegative, and match the number of treated cohorts.
# ------------------------------------------------------------------------------
test_that("Cohort-specific standard errors are computed correctly", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 505)

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 1
	)

	expect_true(is.numeric(result$catt_ses))
	expect_equal(length(result$catt_ses), result$R)
	expect_true(all(result$catt_ses >= 0))
})

# ------------------------------------------------------------------------------
# Test 23: The estimator works when no covariates are provided.
# (processCovs() should issue a warning but continue.)
# ------------------------------------------------------------------------------
test_that("Estimator works with no covariates", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 606)

	# Call fetwfe with no covariate argument.
	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		q = 0.5,
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
# Test 24: processCovs() properly handles an empty covariate vector.
# ------------------------------------------------------------------------------
test_that("processCovs handles empty covariate vector correctly", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 707)

	res <- processCovs(
		df = df,
		units = unique(df$unit),
		unit_var = "unit",
		times = sort(unique(df$time)),
		time_var = "time",
		covs = character(0),
		resp_var = "y",
		T = 10,
		verbose = TRUE
	)

	expect_true(is.data.frame(res$df))
	expect_equal(length(res$covs), 0)
})

# ------------------------------------------------------------------------------
# Test 25: Overall ATT and cohort-specific estimates are finite and numeric.
# ------------------------------------------------------------------------------
test_that("Overall and cohort-specific treatment effects are valid", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 808)

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
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
# Test 26: Second test when no covariates are provided.
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

	result <- fetwfe(
		pdata = pdata, # The panel dataset
		time_var = "time", # The time variable
		unit_var = "unit", # The unit identifier
		treatment = "treated", # The treatment dummy indicator
		response = "response", # The response variable
		q = 0.5 # The L_q penalty for fusion regularization
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
# Test 27: Test that adding ridge regularization works
# ------------------------------------------------------------------------------
test_that("adding ridge regularization to fetwfe works", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 202)

	result <- fetwfe(
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

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		add_ridge = TRUE,
		q = 0.5
	)

	expect_true(is.numeric(result$att_hat))
	expect_false(is.na(result$att_hat))
})

test_that("data application works with ridge penalty", {
	set.seed(23451)

	library(bacondecomp)

	data(divorce)

	# sig_eps_sq and sig_eps_c_sq, calculated in a separate run of `fetwfe(),
	# are provided to speed up the computation of the example
	res <- suppressWarnings(fetwfe(
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
# Test 28: Test that processFactors() converts factor covariates to dummies
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
# Test 29: Test that fetwfe() handles factor covariates appropriately
# ------------------------------------------------------------------------------
test_that("fetwfe handles factor covariates appropriately", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 123)

	# Convert cov2 into a factor with 3 levels.
	set.seed(123)
	df$cov2 <- factor(sample(c("A", "B", "C"), nrow(df), replace = TRUE))

	result <- fetwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
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
# Test 30: Check that tibbles work
# ------------------------------------------------------------------------------
test_that("tibbles work as input to fewtfe", {
	df <- generate_panel_data() # uses N = 30, T = 10 by default

	result <- fetwfe(
		pdata = tibble::as_tibble(df),
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		q = 0.5,
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
		"lambda.max",
		"lambda.max_model_size",
		"lambda.min",
		"lambda.min_model_size",
		"lambda_star",
		"lambda_star_model_size",
		"internal",
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
		c("Cohort", "Estimated TE", "SE", "ConfIntLow", "ConfIntHigh") %in%
			colnames(result$catt_df)
	))
})

# ------------------------------------------------------------------------------
# Test: fetwfe surfaces P_value and selected in catt_df
# ------------------------------------------------------------------------------
test_that("fetwfe surfaces P_value and selected in catt_df", {
	set.seed(2026)
	sim <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(
		sim,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- fetwfeWithSimulatedData(dat, verbose = FALSE)

	expect_true("P_value" %in% colnames(res$catt_df))
	expect_true("selected" %in% colnames(res$catt_df))
	expect_type(res$catt_df$P_value, "double")
	expect_type(res$catt_df$selected, "logical")

	expect_true("att_p_value" %in% names(res))
	expect_true("att_selected" %in% names(res))
	expect_type(res$att_selected, "logical")
	expect_length(res$att_p_value, 1)
	expect_length(res$att_selected, 1)

	# P_value semantics: in [0, 1] when not NA; NA exactly when the cohort
	# is selected out (Estimated TE == 0).
	non_na <- !is.na(res$catt_df$P_value)
	expect_true(all(res$catt_df$P_value[non_na] >= 0))
	expect_true(all(res$catt_df$P_value[non_na] <= 1))
	expect_identical(
		res$catt_df$selected,
		res$catt_df[["Estimated TE"]] != 0
	)
	expect_true(all(is.na(res$catt_df$P_value[!res$catt_df$selected])))
})

# ------------------------------------------------------------------------------
# Test: fetwfe produces at least one selected-out cohort in a sparse simulation
# (sanity assertion: if the underlying grpreg solver ever stops returning
# exact zeros, this test trips and we need to switch to a tolerance-based
# `selected` definition.)
# ------------------------------------------------------------------------------
test_that("fetwfe produces at least one selected-out cohort in a sparse simulation", {
	set.seed(2026)
	sim <- genCoefs(R = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(
		sim,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- fetwfeWithSimulatedData(dat, verbose = FALSE)

	expect_gt(sum(res$catt_df[["Estimated TE"]] == 0), 0)
})

# ------------------------------------------------------------------------------
# Test: order_by = "pvalue" sorts by ascending P_value with NAs last
# ------------------------------------------------------------------------------
test_that("order_by = 'pvalue' sorts CATT by ascending P_value with NAs last", {
	# Use a simulation that yields a mix of selected (non-NA P_value) and
	# selected-out (NA P_value) cohorts so we can verify both the
	# ascending sort AND the NA-last placement, not just one or the other.
	set.seed(2026)
	sim <- genCoefs(R = 4, T = 6, d = 2, density = 0.3, eff_size = 2)
	dat <- simulateData(
		sim,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- fetwfeWithSimulatedData(dat, verbose = FALSE)

	# Guard: the test is only meaningful when both non-NA and NA P_values
	# are present.
	expect_gte(sum(!is.na(res$catt_df$P_value)), 2)
	expect_gte(sum(is.na(res$catt_df$P_value)), 1)

	output_lines <- capture.output(print(res, order_by = "pvalue"))

	header_idx <- grep(
		"Cohort Average Treatment Effects \\(CATT\\):",
		output_lines
	)
	expect_length(header_idx, 1)

	end_idx <- grep("^$|Model Details:", output_lines)
	end_idx <- end_idx[end_idx > header_idx][1]
	cohort_rows <- output_lines[(header_idx + 2):(end_idx - 1)]

	first_tokens <- vapply(
		strsplit(trimws(cohort_rows), "\\s+"),
		`[`,
		character(1),
		1
	)

	expected_order <- order(
		res$catt_df$P_value,
		na.last = TRUE
	)
	expected_cohorts <- as.character(
		res$catt_df$Cohort[expected_order]
	)
	expect_equal(first_tokens, expected_cohorts)

	# Explicit assertion that the ascending-sort property holds across the
	# printed non-NA rows, not just that the overall order matches one
	# specific permutation. This catches a future regression where the
	# sort logic accidentally became descending or unstable.
	non_na_idx <- which(!is.na(res$catt_df$P_value))
	first_non_na_cohort <- as.character(
		res$catt_df$Cohort[non_na_idx][which.min(
			res$catt_df$P_value[non_na_idx]
		)]
	)
	expect_equal(first_tokens[1], first_non_na_cohort)
	# Printed non-NA P_values appear in non-decreasing order.
	printed_non_na_cohorts <- first_tokens[
		first_tokens %in% as.character(res$catt_df$Cohort[non_na_idx])
	]
	printed_p_values <- res$catt_df$P_value[match(
		printed_non_na_cohorts,
		as.character(res$catt_df$Cohort)
	)]
	expect_equal(printed_p_values, sort(printed_p_values))
})

# ------------------------------------------------------------------------------
# Test: fetwfe auto-truncates a no-never-treated panel (default)
# ------------------------------------------------------------------------------
test_that("fetwfe auto-truncates a panel with no never-treated units (default)", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)

	expect_warning(
		res <- fetwfe(
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
	expect_s3_class(res, "fetwfe")
	expect_true(is.finite(res$att_hat))
	expect_s3_class(res$catt_df, "data.frame")
	expect_true(res$T < 10)
	expect_true(res$R < 9)
})

# ------------------------------------------------------------------------------
# Test: fetwfe errors cleanly when truncation would yield < 2 retained cohorts
# ------------------------------------------------------------------------------
test_that("fetwfe errors cleanly when no-never-treated truncation would yield < 2 cohorts", {
	# All units share a single cohort start time, no never-treated.
	# Five units, six time periods, every unit starts treatment at t = 4.
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
		suppressWarnings(fetwfe(
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
# Test: fetwfe with se_type omitted matches se_type = "default" identically
# ------------------------------------------------------------------------------
test_that("fetwfe se_type omitted matches se_type = 'default' identically", {
	sim <- make_se_type_panel()
	res_omit <- fetwfeWithSimulatedData(sim)
	res_def <- fetwfeWithSimulatedData(sim, se_type = "default")
	expect_identical(res_omit$att_se, res_def$att_se)
	expect_identical(res_omit$catt_ses, res_def$catt_ses)
	expect_identical(res_omit$att_hat, res_def$att_hat)
	expect_identical(res_omit$catt_hats, res_def$catt_hats)
})

# ------------------------------------------------------------------------------
# Test: fetwfe with se_type = "cluster" returns finite SEs on a clean panel
# ------------------------------------------------------------------------------
test_that("fetwfe se_type = 'cluster' returns finite SEs on a clean panel", {
	sim <- make_se_type_panel()
	res <- fetwfeWithSimulatedData(sim, se_type = "cluster")
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	selected_ses <- res$catt_ses[res$catt_ses > 0]
	expect_true(length(selected_ses) >= 1)
	expect_true(all(is.finite(selected_ses)))
})

# ------------------------------------------------------------------------------
# Test: under deliberately serially-correlated DGP, cluster SE > default SE
# ------------------------------------------------------------------------------
test_that("fetwfe cluster-robust SE exceeds default SE under AR(1) shocks", {
	mk <- make_ar1_panel()
	res_def <- fetwfeWithSimulatedData(mk, se_type = "default")
	res_cls <- fetwfeWithSimulatedData(mk, se_type = "cluster")
	expect_true(is.finite(res_def$att_se))
	expect_true(is.finite(res_cls$att_se))
	expect_gt(res_def$att_se, 0)
	expect_gt(res_cls$att_se, res_def$att_se)
})

# ------------------------------------------------------------------------------
# Item 7 (issue #84): end-to-end REML × cluster-SE interaction. Most existing
# tests use `*WithSimulatedData()` which passes the simulator's true sigmas
# through, bypassing REML in `estOmegaSqrtInv()`. This test crosses the two
# orthogonal estimation paths — REML variance estimation (triggered by
# `sig_eps_sq = NA, sig_eps_c_sq = NA`) AND cluster-robust SE — to lock in
# the lack of conflict between them. If a future regression in either route
# changes the surface, this assertion fires.
# ------------------------------------------------------------------------------
test_that("fetwfe runs end-to-end with REML estimation + cluster SE", {
	testthat::skip_if_not_installed("lme4")

	set.seed(20260519)
	sim_coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 20260519
	)
	sim <- simulateData(
		sim_coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	# Crucially: leave sig_eps_sq / sig_eps_c_sq at NA so REML runs.
	res <- fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		covs = sim$covs,
		response = sim$response,
		sig_eps_sq = NA,
		sig_eps_c_sq = NA,
		se_type = "cluster",
		verbose = FALSE
	)

	expect_s3_class(res, "fetwfe")
	expect_true(is.finite(res$att_hat))
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	# REML estimated something non-trivial for both variance components.
	expect_gt(res$sig_eps_sq, 0)
	expect_gt(res$sig_eps_c_sq, 0)
	# The cluster route was actually selected (the surface label persists).
	expect_identical(res$se_type, "cluster")
})

# ------------------------------------------------------------------------------
# Test: q >= 1 returns NA SE under both se_type values (oracle property required)
# ------------------------------------------------------------------------------
test_that("fetwfe q >= 1 returns NA SE under both se_type values", {
	sim <- make_se_type_panel()
	res_def <- fetwfeWithSimulatedData(sim, q = 1.5, se_type = "default")
	res_cls <- fetwfeWithSimulatedData(sim, q = 1.5, se_type = "cluster")
	expect_true(is.na(res_def$att_se))
	expect_true(is.na(res_cls$att_se))
})

# ------------------------------------------------------------------------------
# Test: $se_type slot on output reflects the argument value
# ------------------------------------------------------------------------------
test_that("fetwfe $se_type slot reflects the argument value", {
	sim <- make_se_type_panel()
	res_def <- fetwfeWithSimulatedData(sim)
	res_cls <- fetwfeWithSimulatedData(sim, se_type = "cluster")
	expect_identical(res_def$se_type, "default")
	expect_identical(res_cls$se_type, "cluster")
})

# ------------------------------------------------------------------------------
# Test: fetwfe runs end-to-end with zero covariates (covs = c() path).
# The package's d = 0 path is exercised by test-genCoefs.R and
# test-simulateData.R at the simulator level; this test surfaces the same
# capability via the public fetwfeWithSimulatedData wrapper, matching the
# referee request (Report.pdf points 12-13) for a no-covariate
# demonstration alongside the with-covariates tests in this file.
# ------------------------------------------------------------------------------
test_that("fetwfe runs and gives sensible output with zero covariates", {
	set.seed(2026)
	sim_coefs <- genCoefs(
		R = 3,
		T = 6,
		d = 0,
		density = 0.5,
		eff_size = 2,
		seed = 2026
	)
	sim_data <- simulateData(
		sim_coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	res <- fetwfeWithSimulatedData(sim_data)

	expect_s3_class(res, "fetwfe")
	expect_true(is.finite(res$att_hat))
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	expect_length(res$catt_hats, 3)
	expect_length(res$catt_ses, 3)
	expect_s3_class(res$catt_df, "data.frame")
	expect_true(all(
		c(
			"Cohort",
			"Estimated TE",
			"SE",
			"ConfIntLow",
			"ConfIntHigh",
			"P_value",
			"selected"
		) %in%
			colnames(res$catt_df)
	))
	# The covariate-related slots should reflect d = 0:
	expect_identical(res$d, 0L)
})
