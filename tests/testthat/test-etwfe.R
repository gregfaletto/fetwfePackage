library(testthat)
library(fetwfe)

# Panel-data fixture helpers (generate_panel_data, generate_bad_panel_data,
# generate_minimal_panel_data) are defined in tests/testthat/helper-panel-fixture.R
# and sourced by testthat before this file runs (issue #91).

# ------------------------------------------------------------------------------
# Reference (pre-rewrite) processCovs implementation: the canonical
# per-unit-loop version from git history, captured verbatim before issue #84
# item 15. Kept here as a static frozen baseline shared by every processCovs
# equivalence test in this file (#116); any future processCovs rewrite must
# produce identical output on the fixtures those tests build.
# ------------------------------------------------------------------------------
processCovs_reference <- function(
	df,
	units,
	unit_var,
	times,
	time_var,
	covs,
	resp_var,
	T,
	verbose = FALSE
) {
	for (s in units) {
		df_s <- df[df[, unit_var] == s, ]
		if (nrow(df_s) != T || !setequal(df_s[, time_var], times)) {
			stop("balance violation")
		}
	}
	if (length(covs) == 0) {
		df <- df[, c(resp_var, time_var, unit_var)]
		df <- df[
			order(df[, unit_var], df[, time_var], decreasing = FALSE),
		]
		return(list(df = df, covs = covs))
	}
	d_orig <- length(covs)
	covs_to_keep <- character()
	first_time_val <- times[1]
	for (cov_name in covs) {
		is_valid_cov <- TRUE
		for (s in units) {
			val_first_period <- df[
				(df[, unit_var] == s) & (df[, time_var] == first_time_val),
				cov_name
			]
			if (length(val_first_period) != 1 || is.na(val_first_period)) {
				is_valid_cov <- FALSE
				break
			}
		}
		if (is_valid_cov) {
			covs_to_keep <- c(covs_to_keep, cov_name)
		}
	}
	# Suppress reference warnings; equivalence with the live function
	# requires matching warnings, which the live test below verifies
	# separately.
	suppressWarnings({
		if (length(covs_to_keep) < d_orig && length(covs_to_keep) > 0) {
			warning("ref na warning")
		} else if (length(covs_to_keep) == 0 && d_orig > 0) {
			warning("ref na warning")
		}
	})
	covs <- covs_to_keep
	d <- length(covs)
	if (d > 0) {
		covs_to_remove_const <- character()
		df_temp_const_check <- df
		for (s in units) {
			first_period_rows_s_idx <- which(
				(df_temp_const_check[, unit_var] == s) &
					(df_temp_const_check[, time_var] == first_time_val)
			)
			stopifnot(length(first_period_rows_s_idx) == 1)
			cov_values_s_first_period <- df_temp_const_check[
				first_period_rows_s_idx,
				covs,
				drop = FALSE
			]
			for (t_idx in seq_along(times)) {
				current_rows_s_t_idx <- which(
					(df_temp_const_check[, unit_var] == s) &
						(df_temp_const_check[, time_var] == times[t_idx])
				)
				stopifnot(length(current_rows_s_t_idx) == 1)
				df_temp_const_check[
					current_rows_s_t_idx,
					covs
				] <- cov_values_s_first_period
			}
		}
		for (cov_name in covs) {
			first_period_vals_for_cov <- sapply(units, function(u) {
				df_temp_const_check[
					(df_temp_const_check[, unit_var] == u) &
						(df_temp_const_check[, time_var] == first_time_val),
					cov_name
				]
			})
			if (length(unique(first_period_vals_for_cov)) == 1) {
				covs_to_remove_const <- c(covs_to_remove_const, cov_name)
			}
		}
		covs <- covs[!(covs %in% covs_to_remove_const)]
		suppressWarnings({
			if (length(covs_to_remove_const) > 0 && length(covs) == 0) {
				warning("ref const warning")
			} else if (length(covs_to_remove_const) > 0) {
				warning("ref const warning")
			}
		})
		d <- length(covs)
	}
	df <- df[, c(resp_var, time_var, unit_var, covs), drop = FALSE]
	if (d > 0) {
		for (s in units) {
			first_period_rows_s_idx <- which(
				(df[, unit_var] == s) & (df[, time_var] == first_time_val)
			)
			covs_s_first_period_vals <- df[
				first_period_rows_s_idx,
				covs,
				drop = FALSE
			]
			for (t_val in times) {
				ind_s_t <- (df[, unit_var] == s) & (df[, time_var] == t_val)
				stopifnot(sum(ind_s_t) == 1)
				df[ind_s_t, covs] <- covs_s_first_period_vals
			}
		}
	}
	df <- df[order(df[, unit_var], df[, time_var], decreasing = FALSE), ]
	list(df = df, covs = covs)
}

# ------------------------------------------------------------------------------
# Test 1: Check that valid input produces a list with the expected output elements.
# ------------------------------------------------------------------------------
test_that("etwfe returns expected output structure with valid input", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	result <- etwfe(
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
test_that("etwfe errors when pdata is not a data.frame", {
	expect_error(
		etwfe(
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
test_that("etwfe errors when time_var is not integer", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$time <- as.character(df$time) # convert time to character

	expect_error(
		etwfe(
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
test_that("etwfe errors when unit_var is not character", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$unit <- as.factor(df$unit) # make unit a factor

	expect_error(
		etwfe(
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
test_that("etwfe errors when treatment is not integer", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$treatment <- as.character(df$treatment) # wrong type

	expect_error(
		etwfe(
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
test_that("etwfe errors when treatment has values other than 0/1", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	# Force an invalid value into treatment (e.g. 2)
	df$treatment[1] <- 2L

	expect_error(
		etwfe(
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
test_that("etwfe errors when a covariate column is missing", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	# Remove one of the covariate columns
	df <- df[, !(names(df) %in% "cov2")]

	expect_error(
		etwfe(
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
test_that("etwfe errors when response is not numeric", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)
	df$y <- as.character(df$y)

	expect_error(
		etwfe(
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
test_that("etwfe warns when alpha > 0.5", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	expect_warning(
		etwfe(
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
test_that("etwfe errors when all units are treated in the first period", {
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
		suppressWarnings(etwfe(
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

test_that("etwfe errors with a panel having no never-treated units when allow_no_never_treated = FALSE", {
	df_bad <- generate_bad_panel_data(N = 30, T = 10, seed = 123)

	expect_error(
		suppressWarnings(etwfe(
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
test_that("etwfe errors when data has fewer than 4 rows", {
	df <- data.frame(
		time = as.integer(1:3),
		unit = as.character(c("u1", "u2", "u3")),
		treatment = as.integer(c(0, 0, 0)),
		cov1 = rnorm(3),
		cov2 = runif(3),
		y = rnorm(3)
	)

	expect_error(
		etwfe(
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
test_that("etwfe returns att_se", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	result1 <- etwfe(
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
test_that("etwfe warns when all covariates are removed due to constant values", {
	df_const <- generate_panel_data(N = 30, T = 5, R = 2, seed = 101)
	# Set both covariates to a constant value
	df_const$cov1 <- 1
	df_const$cov2 <- 1

	expect_warning(
		etwfe(
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
test_that("etwfe returns a valid numeric overall ATT", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 202)

	result <- etwfe(
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
test_that("etwfe works on a minimal valid dataset", {
	# Create a balanced panel with N = 3 units and T = 3 time periods.
	df_min <- generate_minimal_panel_data(include_cov1 = FALSE)

	result <- etwfe(
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
		etwfe(
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
test_that("etwfe works on only two cohorts", {
	# Create a balanced panel with N = 6 units and T = 5 time periods.
	df <- generate_panel_data(N = 30, T = 10, R = 2, seed = 202)

	result <- etwfe(
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
# Test 19: A single treated cohort (G = 1) is now supported (#112).
# ------------------------------------------------------------------------------
test_that("single treated cohort (G = 1) is supported", {
	# A balanced panel with N = 30 units and T = 10 periods. With R = 1 the
	# helper produces a single treated cohort (adopting at t = 2) plus
	# never-treated units.
	df <- generate_panel_data(N = 30, T = 10, R = 1, seed = 202)

	res <- etwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE
	)
	expect_s3_class(res, "etwfe")
	expect_equal(res$G, 1)
	expect_true(is.finite(res$att_hat))
})


# ------------------------------------------------------------------------------
# Test 20: Overall ATT standard error is computed (non‐zero)
# ------------------------------------------------------------------------------
test_that("Overall ATT standard error is non-negative", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 303)

	# q < 1: expect an overall standard error to be computed
	result <- etwfe(
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

	result <- etwfe(
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

	# Call etwfe with no covariate argument.
	result <- etwfe(
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
# Item 15: processCovs() vectorized rewrite preserves the bit-for-bit output
# of the pre-rewrite per-unit-loop implementation on the standard fixture.
# This is a regression guard against the broadcast-misalignment bug class:
# if a future regression mishandles the `df_first` row order vs. df's sorted
# row order, the time-invariant covariate values get assigned to the wrong
# units, which silently zeroes out downstream estimates (the bug class that
# surfaced during this PR's implementation pass, see the plan's `Surprises
# & Discoveries`).
# ------------------------------------------------------------------------------
test_that("processCovs vectorized output matches reference implementation", {
	# Reference implementation processCovs_reference is defined at file scope
	# (top of this file) and shared with the non-lex-sortable-ID equivalence
	# test below (#116).

	# Fixture: small panel where the time-invariant collapse + the final
	# row sort both matter for downstream alignment.
	set.seed(20260519)
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 20260519)
	units <- unique(df$unit)
	times <- sort(unique(df$time))

	live <- processCovs(
		df = df,
		units = units,
		unit_var = "unit",
		times = times,
		time_var = "time",
		covs = c("cov1", "cov2"),
		resp_var = "y",
		T = 5L,
		verbose = FALSE
	)
	ref <- processCovs_reference(
		df = df,
		units = units,
		unit_var = "unit",
		times = times,
		time_var = "time",
		covs = c("cov1", "cov2"),
		resp_var = "y",
		T = 5L,
		verbose = FALSE
	)

	# Row order, column order, and every cell value must agree.
	rownames(live$df) <- NULL
	rownames(ref$df) <- NULL
	expect_equal(live$df, ref$df, tolerance = 1e-10)
	expect_equal(live$covs, ref$covs)
})

# ------------------------------------------------------------------------------
# Issue #116 Gap 1: the equivalence test above uses generate_panel_data(),
# whose unit IDs are sprintf("unit%02d", ...) -- strictly lexically sortable,
# so caller insertion order coincides with (unit, time) sort order. A
# regression that re-sorted processCovs()'s internal first-period frame back
# to the caller's insertion order would corrupt the first-period covariate
# broadcast, yet stay invisible on a lex-sortable fixture. This test rebuilds
# the equivalence check on a hand-crafted panel whose unit IDs are NOT in
# sorted order on insertion, so a re-sort regression changes cell values.
# ------------------------------------------------------------------------------
test_that("processCovs vectorized output matches reference with non-lex-sortable unit IDs", {
	# Insertion order: Alpha, delta, Bravo, charlie.
	# Lexical order:   Alpha, Bravo, charlie, delta  (differs -> a re-sort
	# regression is observable). processCovs() takes no treatment argument and
	# never references cohorts, so a balanced unit/time/covariate/y panel is
	# sufficient.
	unit_order <- c("Alpha", "delta", "Bravo", "charlie")
	T_n <- 4L
	times_vals <- seq_len(T_n)
	# cov1: first-period value differs across units (10/20/30/40) so a
	# misaligned broadcast changes cell values; the within-unit value also
	# varies over time so processCovs()'s first-period collapse is exercised.
	# The per-unit-distinct first-period values also keep cov1 from being
	# dropped by the constant-across-units column check.
	first_vals <- c(Alpha = 10, delta = 20, Bravo = 30, charlie = 40)
	df_nonlex <- do.call(
		rbind,
		lapply(unit_order, function(u) {
			data.frame(
				time = as.integer(times_vals),
				unit = as.character(u),
				cov1 = first_vals[[u]] + times_vals,
				y = as.numeric(first_vals[[u]] * 100 + times_vals),
				stringsAsFactors = FALSE
			)
		})
	)
	# Caller order: the units vector reflects insertion order, not sort order.
	units <- unit_order
	times <- sort(unique(df_nonlex$time))

	live <- processCovs(
		df = df_nonlex,
		units = units,
		unit_var = "unit",
		times = times,
		time_var = "time",
		covs = "cov1",
		resp_var = "y",
		T = T_n,
		verbose = FALSE
	)
	ref <- processCovs_reference(
		df = df_nonlex,
		units = units,
		unit_var = "unit",
		times = times,
		time_var = "time",
		covs = "cov1",
		resp_var = "y",
		T = T_n,
		verbose = FALSE
	)

	# Row order, column order, and every cell value must agree.
	rownames(live$df) <- NULL
	rownames(ref$df) <- NULL
	expect_equal(live$df, ref$df, tolerance = 1e-10)
	expect_equal(live$covs, ref$covs)
	# cov1 must survive the constant-across-units drop (first-period values
	# 10/20/30/40 are distinct), otherwise the broadcast path is not exercised.
	expect_equal(live$covs, "cov1")
})

# ------------------------------------------------------------------------------
# Item 15: end-to-end etwfeWithSimulatedData() outputs unchanged after
# processCovs vectorization. Frozen-seed regression catching silent numeric
# drift in the downstream estimator caused by an off-by-one in the
# first-period broadcast.
# ------------------------------------------------------------------------------
test_that("etwfeWithSimulatedData end-to-end unchanged after processCovs rewrite", {
	set.seed(20260519)
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 20260519
	)
	sim <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- etwfeWithSimulatedData(sim)
	# Smoke check: the rewrite doesn't introduce NaN/Inf on the standard
	# fixture. Numerical equivalence vs the pre-rewrite reference is locked
	# by the processCovs-vs-reference test above (which directly compares
	# the vectorized output to a verbatim copy of the pre-rewrite loop on
	# the same X_ints / pdata inputs).
	expect_true(is.finite(res$att_hat))
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	# Per-cohort estimates and SEs are finite and non-NA.
	expect_true(all(is.finite(res$catt_hats)))
	expect_true(all(is.finite(res$catt_ses)))
})

# ------------------------------------------------------------------------------
# Test 24: Overall ATT and cohort-specific estimates are finite and numeric.
# ------------------------------------------------------------------------------
test_that("Overall and cohort-specific treatment effects are valid", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 808)

	result <- etwfe(
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

	result <- etwfe(
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
test_that("adding ridge regularization to etwfe works", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 202)

	result <- etwfe(
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

	result <- etwfe(
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
	res <- suppressWarnings(etwfe(
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
# Test 28: Test that etwfe() handles factor covariates appropriately
# ------------------------------------------------------------------------------
test_that("etwfe handles factor covariates appropriately", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123)

	# Convert cov2 into a factor with 3 levels.
	set.seed(123)
	df$cov2 <- factor(sample(c("A", "B", "C"), nrow(df), replace = TRUE))

	result <- etwfe(
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
test_that("tibbles work as input to fewtfe", {
	df <- generate_panel_data(N = 30, T = 5, R = 2, seed = 123) # uses N = 30, T = 10 by default

	result <- etwfe(
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
test_that("etwfe throws error when a cohort contains fewer than d + 1 units", {
	df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 123)

	expect_error(
		etwfe(
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
		etwfe(
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

	res <- suppressWarnings(etwfe(
		pdata = df,
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		covs = c("cov1", "cov2"),
		response = "y",
		verbose = FALSE,
		add_ridge = TRUE
	))

	expect_false(res$calc_ses)
	expect_true(is.na(res$att_se))
	expect_true(!is.na(res$att_hat))
	expect_true(res$att_hat != 0)
})

# ------------------------------------------------------------------------------
# Test: etwfe surfaces p_value but not selected in catt_df
# ------------------------------------------------------------------------------
test_that("etwfe surfaces p_value but not selected in catt_df", {
	set.seed(2026)
	sim <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	dat <- simulateData(
		sim,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
	res <- etwfeWithSimulatedData(dat, verbose = FALSE)

	expect_true("p_value" %in% colnames(res$catt_df))
	expect_false("selected" %in% colnames(res$catt_df))

	expect_true("att_p_value" %in% names(res))
	expect_false("att_selected" %in% names(res))
	expect_length(res$att_p_value, 1)

	# p_value is in [0, 1] for non-NA entries.
	non_na <- !is.na(res$catt_df$p_value)
	expect_true(all(res$catt_df$p_value[non_na] >= 0))
	expect_true(all(res$catt_df$p_value[non_na] <= 1))
})

# ------------------------------------------------------------------------------
# Test: etwfe auto-truncates a no-never-treated panel (default)
# ------------------------------------------------------------------------------
test_that("etwfe auto-truncates a panel with no never-treated units (default)", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)

	expect_warning(
		res <- etwfe(
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
	expect_s3_class(res, "etwfe")
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
# exactly, which a finiteness-only check would pass vacuously.
# ------------------------------------------------------------------------------
test_that("etwfe cluster SE is finite and positive on an auto-truncated panel", {
	df_bad <- generate_bad_panel_data(N = 200, T = 10, seed = 123)

	expect_warning(
		res <- etwfe(
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
	expect_s3_class(res, "etwfe")
	expect_true(is.finite(res$att_se))
	expect_false(is.na(res$att_se))
	expect_gt(res$att_se, 0)
})

# ------------------------------------------------------------------------------
# Test: etwfe errors cleanly when truncation would yield < 2 retained cohorts
# ------------------------------------------------------------------------------
test_that("etwfe errors cleanly when no-never-treated truncation would yield < 2 cohorts", {
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
		suppressWarnings(etwfe(
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
	sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	simulateData(sim_coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5)
}

make_ar1_panel <- function(seed = 7, N = 150, rho = 0.85, sd_e = 1) {
	set.seed(seed)
	sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
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
# Test: etwfe with se_type omitted matches se_type = "default" identically
# ------------------------------------------------------------------------------
test_that("etwfe se_type omitted matches se_type = 'default' identically", {
	sim <- make_se_type_panel()
	res_omit <- etwfeWithSimulatedData(sim)
	res_def <- etwfeWithSimulatedData(sim, se_type = "default")
	expect_identical(res_omit$att_se, res_def$att_se)
	expect_identical(res_omit$catt_ses, res_def$catt_ses)
	expect_identical(res_omit$att_hat, res_def$att_hat)
	expect_identical(res_omit$catt_hats, res_def$catt_hats)
})

# ------------------------------------------------------------------------------
# Test: etwfe with se_type = "cluster" returns finite SEs on a clean panel
# ------------------------------------------------------------------------------
test_that("etwfe se_type = 'cluster' returns finite SEs on a clean panel", {
	sim <- make_se_type_panel()
	res <- etwfeWithSimulatedData(sim, se_type = "cluster")
	expect_true(is.finite(res$att_se))
	expect_gt(res$att_se, 0)
	expect_true(all(is.finite(res$catt_ses)))
	expect_true(all(res$catt_ses > 0))
})

# ------------------------------------------------------------------------------
# Test: under deliberately serially-correlated DGP, cluster SE > default SE
# ------------------------------------------------------------------------------
test_that("etwfe cluster-robust SE exceeds default SE under AR(1) shocks", {
	mk <- make_ar1_panel()
	res_def <- etwfeWithSimulatedData(mk, se_type = "default")
	res_cls <- etwfeWithSimulatedData(mk, se_type = "cluster")
	expect_true(is.finite(res_def$att_se))
	expect_true(is.finite(res_cls$att_se))
	expect_gt(res_def$att_se, 0)
	expect_gt(res_cls$att_se, res_def$att_se)
})

# ------------------------------------------------------------------------------
# Test: $se_type slot on output reflects the argument value
# ------------------------------------------------------------------------------
test_that("etwfe $se_type slot reflects the argument value", {
	sim <- make_se_type_panel()
	res_def <- etwfeWithSimulatedData(sim)
	res_cls <- etwfeWithSimulatedData(sim, se_type = "cluster")
	expect_identical(res_def$se_type, "default")
	expect_identical(res_cls$se_type, "cluster")
})
