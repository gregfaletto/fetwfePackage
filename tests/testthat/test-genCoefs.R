library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# Helper functions: expected dimensions for genCoefs() output.
# ------------------------------------------------------------------------------
compute_num_treats <- function(T, R) {
  T * R - (R * (R + 1)) / 2
}

compute_p <- function(T, R, d) {
  num_treats <- compute_num_treats(T, R)
  R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
}

# ------------------------------------------------------------------------------
# Test 1: Check that genCoefs returns a list with components beta and theta,
#         and that both are numeric vectors.
# ------------------------------------------------------------------------------
test_that("genCoefs returns expected output structure", {
  res <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2)
  
  expect_type(res, "list")
  expect_true("beta" %in% names(res))
  expect_true("theta" %in% names(res))
  
  expect_true(is.numeric(res$beta))
  expect_true(is.numeric(res$theta))
})

# ------------------------------------------------------------------------------
# Test 2: Check that the length of beta and theta equals the expected number.
# ------------------------------------------------------------------------------
test_that("genCoefs returns beta and theta of correct length", {
  R <- 5; T <- 30; d <- 12
  expected_length <- compute_p(T, R, d)
  
  res <- genCoefs(R = R, T = T, d = d, density = 0.1, eff_size = 2)
  expect_equal(length(res$beta), expected_length)
  expect_equal(length(res$theta), expected_length)
})

# ------------------------------------------------------------------------------
# Test 3: Reproducibility test: With the same seed, the output should be identical.
# ------------------------------------------------------------------------------
test_that("genCoefs is reproducible with the same seed", {
  set.seed(123)
  res1 <- genCoefs(R = 5, T = 30, density = 0.1, eff_size = 2, d = 12)
  
  set.seed(123)
  res2 <- genCoefs(R = 5, T = 30, density = 0.1, eff_size = 2, d = 12)
  
  expect_equal(res1$beta, res2$beta)
  expect_equal(res1$theta, res2$theta)
})

# ------------------------------------------------------------------------------
# Test 4: Verify that the beta vector returned by genCoefs() is a valid input
#         for genRandomData().
# ------------------------------------------------------------------------------
test_that("beta from genCoefs is a valid input for genRandomData", {
  R <- 5; T <- 30; d <- 12
  res_coefs <- genCoefs(R = R, T = T, d = d, density = 0.1, eff_size = 2)
  
  N <- 120
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  
  # Use res_coefs$beta as the beta argument for genRandomData().
  sim_data <- genRandomData(N = N, T = T, R = R, d = d,
                            sig_eps_sq = sig_eps_sq, sig_eps_c_sq = sig_eps_c_sq,
                            beta = res_coefs$beta, seed = 123, gen_ints = TRUE)
  
  expect_type(sim_data, "list")
  # Check that the returned list contains key elements (e.g., pdata, X, y)
  expected_names <- c("pdata", "X", "y", "coefs", "N", "T", "R", "d", "p")
  for (nm in expected_names) {
    expect_true(nm %in% names(sim_data))
  }
  
  # Check that the reported p equals the expected dimension.
  expect_equal(sim_data$p, compute_p(T, R, d))
})

# ------------------------------------------------------------------------------
# Test 5: Check that the returned theta vector has approximate sparsity equal to density.
# ------------------------------------------------------------------------------
test_that("genCoefs produces theta with approximately correct sparsity", {
  R <- 5; T <- 30; d <- 12; density <- 0.1
  res <- genCoefs(R = R, T = T, d = d, density = density, eff_size = 2)
  theta <- res$theta
  total <- length(theta)
  nonzero <- sum(theta != 0)
  prop_nonzero <- nonzero / total
  
  # Allow some tolerance given randomness; here we allow a deviation of 0.03.
  expect_true(abs(prop_nonzero - density) < 0.03)
})


# ------------------------------------------------------------------------------
# Test 6: d = 0 case.
# ------------------------------------------------------------------------------
test_that("genRandomData and fetwfe work when d = 0", {
  # Set parameters for simulation
  N <- 120
  T_val <- 30
  R_val <- 5
  d_val <- 0
  sig_eps_sq <- 5
  sig_eps_c_sq <- 5
  
  # For d = 0, the expected number of columns in the design matrix is:
  # p = R + (T - 1) + num_treats, where num_treats = T * R - (R*(R+1))/2.
  num_treats <- T_val * R_val - (R_val * (R_val + 1)) / 2
  expected_p <- R_val + (T_val - 1) + num_treats
  
  # Generate coefficients with no covariates using genCoefs.
  coefs <- genCoefs(R = R_val, T = T_val, density = 0.1, eff_size = 2, d = d_val)
  
  # Use genRandomData() with gen_ints = FALSE.
  sim_data <- genRandomData(
    N = N, T = T_val, R = R_val, d = d_val, 
    sig_eps_sq = sig_eps_sq, sig_eps_c_sq = sig_eps_c_sq,
    beta = coefs$beta, seed = 123, gen_ints = FALSE
  )
  
  # Check that cov_names is empty since d = 0.
  expect_equal(length(sim_data$cov_names), 0)
  
  # Check that the design matrix X has the expected dimensions.
  expect_equal(dim(sim_data$X), c(N * T_val, expected_p))
  
  # Now run fetwfe() on the simulated panel data.
  result <- fetwfe(
    pdata     = sim_data$pdata,
    time_var  = sim_data$time_var,
    unit_var  = sim_data$unit_var,
    treatment = sim_data$treatment,
    covs      = sim_data$cov_names,  # should be empty
    response  = sim_data$response,
    sig_eps_sq = sig_eps_sq,
    sig_eps_c_sq = sig_eps_c_sq,
    verbose   = FALSE
  )
  
  # Verify that the returned d (number of covariates) is 0.
  expect_equal(result$d, 0)
  
  # Also, overall ATT (att_hat) should be numeric and non-missing.
  expect_true(is.numeric(result$att_hat))
  expect_false(is.na(result$att_hat))
})