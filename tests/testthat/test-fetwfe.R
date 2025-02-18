library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Helper function: generate a synthetic balanced panel data set.
# We now use N=30 units and T=10 periods so that the degrees‐of‐freedom condition
# (N * (T - 1) - p > 0) holds.
# For each unit, we randomly decide whether the unit is never treated (first_treat = Inf)
# or, if treated, assign a first treatment time (ensuring that no unit is treated in time 1).
# ------------------------------------------------------------------------------
generate_panel_data <- function(N = 30, T = 10, R = 9, seed = 123) {
  set.seed(seed)
  stopifnot(R <= T - 1)
  # Create a vector of unit IDs (as characters)
  unit_ids <- sprintf("unit%02d", 1:N)
  time_vals <- 1:T
  
  # For each unit, decide the first treatment time:
  # With probability 0.5 the unit is never treated (first_treat = Inf);
  # Otherwise, sample a first treatment time uniformly between 2 and T.
  first_treat <- sapply(unit_ids, function(u) {
    if (runif(1) < 0.5) {
      Inf
    } else {
      if(R > 1){
        sample(2:(R + 1), 1)
      } else{
        2
      }
      
    }
  })
  
  # Build panel data (each unit appears for every time period)
  df <- do.call(rbind, lapply(seq_along(unit_ids), function(i) {
    unit <- unit_ids[i]
    ft <- first_treat[i]
    data.frame(
      time      = as.integer(time_vals),            # must be integer
      unit      = as.character(unit),               # must be character
      treatment = as.integer(ifelse(time_vals >= ft, 1, 0)),  # must be integer 0/1
      cov1      = rnorm(T),
      cov2      = runif(T),
      y         = rnorm(T)                          # outcome (numeric)
    )
  }))
  
  # The idCohorts function automatically removes any unit that was treated in period 1.
  # (This is simulated by removing any row with time == 1 and treatment == 1.)
  df <- df[!(df$time == 1 & df$treatment == 1), ]
  
  # Order rows by unit then time
  df <- df[order(df$unit, df$time), ]
  rownames(df) <- NULL
  return(df)
}

generate_bad_panel_data <- function(N = 30, T = 10, seed = 123) {
  set.seed(seed)
  # Create a vector of unit IDs (as characters)
  unit_ids <- sprintf("unit%02d", 1:N)
  time_vals <- 1:T
  
  # For each unit, decide the first treatment time:
  # Always sample a first treatment time uniformly between 2 and T
  # No never-treated units.
  first_treat <- sapply(unit_ids, function(u) {
      sample(2:T, 1)
  })
  
  # Build panel data (each unit appears for every time period)
  df <- do.call(rbind, lapply(seq_along(unit_ids), function(i) {
    unit <- unit_ids[i]
    ft <- first_treat[i]
    data.frame(
      time      = as.integer(time_vals),            # must be integer
      unit      = as.character(unit),               # must be character
      treatment = as.integer(ifelse(time_vals >= ft, 1, 0)),  # must be integer 0/1
      cov1      = rnorm(T),
      cov2      = runif(T),
      y         = rnorm(T)                          # outcome (numeric)
    )
  }))
  
  # Order rows by unit then time
  df <- df[order(df$unit, df$time), ]
  rownames(df) <- NULL
  return(df)
}

generate_minimal_panel_data <- function(seed = 123) {
  set.seed(seed)

  N <- 3
  T <- 3
  # Create a vector of unit IDs (as characters)
  unit_ids <- sprintf("unit%02d", 1:N)
  time_vals <- 1:T
  
  # For each unit, decide the first treatment time:
  # one treated, one untreated unit.
  first_treat <- sapply(unit_ids, function(u) {
    if (u  == "unit01") {
      Inf
    } else if (u  == "unit02") {
      2
    } else{
      3
    }
  })
  
  # Build panel data (each unit appears for every time period)
  df <- do.call(rbind, lapply(seq_along(unit_ids), function(i) {
    unit <- unit_ids[i]
    ft <- first_treat[i]
    data.frame(
      time      = as.integer(time_vals),            # must be integer
      unit      = as.character(unit),               # must be character
      treatment = as.integer(ifelse(time_vals >= ft, 1, 0)),  # must be integer 0/1
      cov1      = rnorm(T),
      y         = rnorm(T)                          # outcome (numeric)
    )
  }))
  
  # The idCohorts function automatically removes any unit that was treated in period 1.
  # (This is simulated by removing any row with time == 1 and treatment == 1.)
  df <- df[!(df$time == 1 & df$treatment == 1), ]
  
  # Order rows by unit then time
  df <- df[order(df$unit, df$time), ]
  rownames(df) <- NULL
  return(df)
}

# ------------------------------------------------------------------------------
# Test 1: Check that valid input produces a list with the expected output elements.
# ------------------------------------------------------------------------------
test_that("fetwfe returns expected output structure with valid input", {
  df <- generate_panel_data()  # uses N = 30, T = 10 by default
  
  result <- fetwfe(
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    q         = 0.5,
    verbose   = FALSE
  )
  
  expect_type(result, "list")
  
  # Check that the returned list contains the main expected names.
  expected_names <- c("att_hat", "att_se", "catt_hats", "catt_ses", "catt_df",
                      "beta_hat", "treat_inds", "treat_int_inds",
                      "sig_eps_sq", "sig_eps_c_sq",
                      "lambda.max", "lambda.max_model_size",
                      "lambda.min", "lambda.min_model_size",
                      "lambda_star", "lambda_star_model_size",
                      "X_ints", "y", "X_final", "y_final",
                      "N", "T", "R", "d", "p")
  for (nm in expected_names) {
    expect_true(nm %in% names(result))
  }
  
  # Also check that catt_df is a data frame with the expected column names.
  expect_s3_class(result$catt_df, "data.frame")
  expect_true(all(c("Cohort", "Estimated TE", "SE", "ConfIntLow", "ConfIntHigh") %in% colnames(result$catt_df)))
})

# ------------------------------------------------------------------------------
# Test 2: Error when pdata is not a data.frame
# ------------------------------------------------------------------------------
test_that("fetwfe errors when pdata is not a data.frame", {
  expect_error(
    fetwfe(
      pdata     = "not a data.frame",
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "is.data.frame"
  )
})

# ------------------------------------------------------------------------------
# Test 3: Error when time_var column is not integer.
# ------------------------------------------------------------------------------
test_that("fetwfe errors when time_var is not integer", {
  df <- generate_panel_data()
  df$time <- as.character(df$time)  # convert time to character
  
  expect_error(
    fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "is.integer"
  )
})

# ------------------------------------------------------------------------------
# Test 4: Error when unit_var column is not character.
# ------------------------------------------------------------------------------
test_that("fetwfe errors when unit_var is not character", {
  df <- generate_panel_data()
  df$unit <- as.factor(df$unit)  # make unit a factor
  
  expect_error(
    fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "is.character"
  )
})

# ------------------------------------------------------------------------------
# Test 5: Error when treatment column is not integer.
# ------------------------------------------------------------------------------
test_that("fetwfe errors when treatment is not integer", {
  df <- generate_panel_data()
  df$treatment <- as.character(df$treatment)  # wrong type
  
  expect_error(
    fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "is.integer"
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
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "all\\(pdata\\[, treatment\\] %in% c\\(0, 1\\)\\)"
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
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "all\\(covs %in% colnames\\(pdata\\)\\)"
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
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "is.numeric"
  )
})

# ------------------------------------------------------------------------------
# Test 9: Warning when alpha > 0.5.
# ------------------------------------------------------------------------------
test_that("fetwfe warns when alpha > 0.5", {
  df <- generate_panel_data()
  
  expect_warning(
    fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y",
      alpha     = 0.6,
      verbose   = FALSE
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
    time      = rep(1:5, times = 5),
    unit      = rep(sprintf("unit%02d", 1:5), each = 5),
    treatment = rep(1L, 25),
    cov1      = rnorm(25),
    cov2      = runif(25),
    y         = rnorm(25)
  )
  
  expect_error(
    suppressWarnings(fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y",
      verbose   = FALSE
    )),
    "All units were treated in the first time period; estimating treatment effects is not possible"
  )
})

# ------------------------------------------------------------------------------
# Test 11: Error when there are no never-treated units.
# ------------------------------------------------------------------------------

test_that("fetwfe errors with a panel having no never-treated units", {
  df_bad <- generate_bad_panel_data(N = 30, T = 10, seed = 123)
  
  expect_error(
    suppressWarnings(fetwfe(
      pdata     = df_bad,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y",
      verbose   = FALSE
    )),
    "No never-treated units detected in data to fit model; estimating treatment effects is not possible"
  )
})

# ------------------------------------------------------------------------------
# Test 12: Error when the data has too few rows (less than 4).
# ------------------------------------------------------------------------------
test_that("fetwfe errors when data has fewer than 4 rows", {
  df <- data.frame(
    time      = as.integer(1:3),
    unit      = as.character(c("u1", "u2", "u3")),
    treatment = as.integer(c(0, 0, 0)),
    cov1      = rnorm(3),
    cov2      = runif(3),
    y         = rnorm(3)
  )
  
  expect_error(
    fetwfe(
      pdata     = df,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y"
    ),
    "nrow\\(pdata\\) >= 4"
  )
})

# ------------------------------------------------------------------------------
# Test 13: Test that optional lambda parameters are handled.
# ------------------------------------------------------------------------------
test_that("fetwfe returns expected output when lambda parameters are supplied", {
  df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 456)
  
  result <- fetwfe(
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    lambda.max = 5,
    lambda.min = 0.5,
    nlambda    = 50,
    verbose   = FALSE
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
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    q         = 0.5,
    verbose   = FALSE
  )
  expect_false(is.na(result1$att_se))
  
  result2 <- fetwfe(
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    q         = 2,  # ridge-like penalty, where SE is not computed
    verbose   = FALSE
  )
  expect_true(is.na(result2$att_se))
})

# ------------------------------------------------------------------------------
# Test 15: Test that an error is thrown when all covariates are constant.
# In this case, processCovs should remove all covariates and then stop.
# ------------------------------------------------------------------------------
test_that("fetwfe errors when all covariates are removed due to constant values", {
  df_const <- generate_panel_data(N = 30, T = 10, R = 9, seed = 101)
  # Set both covariates to a constant value
  df_const$cov1 <- 1
  df_const$cov2 <- 1
  
  expect_error(
    fetwfe(
      pdata     = df_const,
      time_var  = "time",
      unit_var  = "unit",
      treatment = "treatment",
      covs      = c("cov1", "cov2"),
      response  = "y",
      verbose   = FALSE
    ),
    "All covariates were removed"  # expecting an error message including this substring
  )
})

# ------------------------------------------------------------------------------
# Test 16: Test that the overall ATT (att_hat) is numeric and non-missing.
# ------------------------------------------------------------------------------
test_that("fetwfe returns a valid numeric overall ATT", {
  df <- generate_panel_data(N = 30, T = 10, R = 9, seed = 202)
  
  result <- fetwfe(
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    verbose   = FALSE
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
    pdata     = df_min,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1"),
    response  = "y",
    verbose   = FALSE,
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
    pdata     = df_min,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1"),
    response  = "y",
    verbose   = FALSE
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
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    verbose   = FALSE
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
    pdata     = df,
    time_var  = "time",
    unit_var  = "unit",
    treatment = "treatment",
    covs      = c("cov1", "cov2"),
    response  = "y",
    verbose   = FALSE
  ),
  "Only one treated cohort detected in data. Currently fetwfe only supports data sets with at least two treated cohorts."
  )
})
