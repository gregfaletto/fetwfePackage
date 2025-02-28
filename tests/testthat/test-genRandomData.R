library(testthat)
library(fetwfe)

# Helper functions to compute expected dimensions:
compute_num_treats <- function(T, R) {
  T * R - (R * (R + 1)) / 2
}

# For full interactions:
compute_p_int <- function(T, R, d) {
  num_treats <- compute_num_treats(T, R)
  R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
}

# For no interactions:
compute_p_no_int <- function(T, R, d) {
  num_treats <- compute_num_treats(T, R)
  R + (T - 1) + d + num_treats
}

# ------------------------------------------------------------------------------
# Test 1: Output structure and dimensions when gen_ints = TRUE
# ------------------------------------------------------------------------------
test_that("genRandomData (with interactions) returns expected output", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_int <- compute_p_int(T_val, R_val, d_val)
  beta_int <- rnorm(p_int)
  
  res <- genRandomData(N=N,
    T=T_val,
    R=R_val,
    d=d_val, 
    sig_eps_sq=sig_eps_sq,
    sig_eps_c_sq=sig_eps_c_sq,
    beta=beta_int,
    seed = 123,
    gen_ints = TRUE)
  
  # Check list names (X now instead of X_int)
  expected_names <- c("X", "y", "coefs", "first_inds", "N_UNTREATED",
                      "assignments", "indep_assignments", "actual_cohort_tes",
                      "att_true", "p", "N", "T", "R", "d", "sig_eps_sq", "sig_eps_c_sq")
  for (nm in expected_names) {
    expect_true(nm %in% names(res))
  }
  
  # Check dimensions of X and y
  expect_equal(dim(res$X), c(N * T_val, p_int))
  expect_equal(length(res$y), N * T_val)
  
  # Check that att_true equals the mean of actual_cohort_tes
  expect_equal(res$att_true, as.numeric(mean(res$actual_cohort_tes)), tolerance = 1e-6)
})

# ------------------------------------------------------------------------------
# Test 2: Output structure and dimensions when gen_ints = FALSE
# ------------------------------------------------------------------------------
test_that("genRandomData (without interactions) returns expected output", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_no_int <- compute_p_no_int(T_val, R_val, d_val)

  p_int <- compute_p_int(T_val, R_val, d_val)
  beta_int <- rnorm(p_int)
  
  res <- genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_int, seed = 123, gen_ints = FALSE)
  
  expect_equal(dim(res$X), c(N * T_val, p_no_int))
  expect_equal(length(res$y), N * T_val)
  expect_equal(res$p, p_no_int)
  
  # Check that the treatment dummy columns are at the end of the design matrix.
  num_treats <- compute_num_treats(T_val, R_val)
  base_cols <- if (d_val > 0) {
    R_val + (T_val - 1) + d_val + d_val * R_val + d_val * (T_val - 1)
  } else {
    R_val + (T_val - 1)
  }
  treat_inds <- seq(from = base_cols + 1, length.out = num_treats)
  # Use getActualCohortTes() to compute actual cohort effects.
  actual_cohort_tes <- getActualCohortTes(R_val, res$first_inds, treat_inds, beta_int, num_treats)
  expect_equal(res$actual_cohort_tes, actual_cohort_tes)
  
  # Check that overall ATT equals mean of actual cohort effects.
  expect_equal(res$att_true, as.numeric(mean(actual_cohort_tes)), tolerance = 1e-6)
})

# ------------------------------------------------------------------------------
# Test 3: Reproducibility test (seed)
# ------------------------------------------------------------------------------
test_that("genRandomData is reproducible with same seed", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_no_int <- compute_p_no_int(T_val, R_val, d_val)

  p_int <- compute_p_int(T_val, R_val, d_val)
  beta_int <- rnorm(p_int)
  
  res1 <- genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_int, seed = 456, gen_ints = FALSE)
  res2 <- genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_int, seed = 456, gen_ints = FALSE)
  
  expect_equal(res1$X, res2$X)
  expect_equal(res1$y, res2$y)
  expect_equal(res1$first_inds, res2$first_inds)
  expect_equal(res1$assignments, res2$assignments)
  expect_equal(res1$indep_assignments, res2$indep_assignments)
  expect_equal(res1$actual_cohort_tes, res2$actual_cohort_tes)
  expect_equal(res1$att_true, res2$att_true)
})

# ------------------------------------------------------------------------------
# Test 4: Error when beta has incorrect length for gen_ints = TRUE
# ------------------------------------------------------------------------------
test_that("genRandomData errors when beta has wrong length (with interactions)", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_int <- compute_p_int(T_val, R_val, d_val)
  beta_wrong <- rnorm(p_int - 1)
  
  expect_error(
    genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_wrong, seed = 123, gen_ints = TRUE),
    "length\\(beta\\) must be"
  )
})

# ------------------------------------------------------------------------------
# Test 5: Error when beta has incorrect length for gen_ints = FALSE
# ------------------------------------------------------------------------------
test_that("genRandomData errors when beta has wrong length (without interactions)", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_no_int <- compute_p_no_int(T_val, R_val, d_val)
  beta_wrong <- rnorm(p_no_int + 2)
  
  expect_error(
    genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_wrong, seed = 123, gen_ints = FALSE),
    "length\\(beta\\) must be"
  )
})

# ------------------------------------------------------------------------------
# Test 6: Check that outputs (X and y) are non-missing.
# ------------------------------------------------------------------------------
test_that("genRandomData returns non-missing X and y", {
  N <- 120; T_val <- 30; R_val <- 5; d_val <- 12
  sig_eps_sq <- 5; sig_eps_c_sq <- 5
  p_no_int <- compute_p_no_int(T_val, R_val, d_val)

  p_int <- compute_p_int(T_val, R_val, d_val)
  beta_int <- rnorm(p_int)
  
  res <- genRandomData(N, T_val, R_val, d_val, sig_eps_sq, sig_eps_c_sq, beta_int, seed = 789, gen_ints = FALSE)
  expect_true(all(!is.na(res$X)))
  expect_true(all(!is.na(res$y)))
})
