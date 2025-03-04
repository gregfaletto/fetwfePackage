generateBaseEffects <- function(N, d, T, R, distribution = "gaussian"){
  ret <- genCohortTimeFE(N, T, R, d)
  
  if (distribution == "gaussian") {
    X <- matrix(rnorm(N * d), nrow = N, ncol = d)
  } else if (distribution == "uniform") {
    # Generate U(-sqrt(3), sqrt(3)) so that variance = 1 (same as standard normal)
    a <- sqrt(3)
    X <- matrix(runif(N * d, min = -a, max = a), nrow = N, ncol = d)
  } else {
    stop("Unsupported distribution. Please choose 'gaussian' or 'uniform'.")
  }

  stopifnot(ncol(ret$cohort_fe) == R)
  
  X_long <- X[rep(1:N, each = T), ]
  return(list(cohort_fe = ret$cohort_fe, 
              time_fe = ret$time_fe, 
              X_long = X_long,
              assignments = ret$assignments,
              cohort_inds = ret$inds))
}


genCohortTimeFE <- function(N, T, R, d){
    # The observations will be arranged row-wise as blocks of T, one unit at
    # a time. So the first T rows correspond to all observations from the first
    # unit, and so on. Therefore the first N_UNTREATED*T rows contain all of
    # the observations from untreated units, then the next N_PER_COHORT*T
    # rows contain the observations from the first cohort, and so on.

    # Each cohort, as well as the untreated group, will have at least d + 1
    # units. The remaining units will be allocated to each of these R + 1
    # groups uniformly at random.
    stopifnot(N >= (R + 1)*(d + 1))

    # Generate cohort assignments
    assignments <- genAssignments(N, R)
    
    # Cohort fixed effects
    cohort_fe <- matrix(0, N*T, R)
    first_ind_r <- assignments[1]*T + 1

    inds <- list()

    for(r in 1:R){
        stopifnot(all(cohort_fe[, r] == 0))

        last_ind_r <- first_ind_r + assignments[r + 1] * T - 1

        stopifnot(last_ind_r <= N * T)
        stopifnot(length(first_ind_r:last_ind_r) == assignments[r + 1]*T)
        stopifnot(all(cohort_fe[first_ind_r:last_ind_r, ] == 0))

        # Now add cohort fixed effects in the rth column for these rows
        cohort_fe[first_ind_r:last_ind_r, r] <- rep(1, assignments[r + 1]*T)
        inds[[r]] <- first_ind_r:last_ind_r
        first_ind_r <- last_ind_r + 1

        stopifnot(length(inds[[r]]) == assignments[r + 1]*T)
    }

    stopifnot(last_ind_r == N*T)

    # Time fixed effects: only do 2 through T.
    time_fe <- matrix(0, N*T, T - 1)
    # No fixed effect for first time
    for(t in 1:(T - 1)){
        stopifnot(all(time_fe[, t] == 0))
        # For each time (except the first), we have N rows of observations which
        # need to be assigned a time fixed effect. Identify those rows:
        # The (t + 1)st row of every block should have a 1 in the (t + 1)st
        # column
        rows_t <- (0:(N - 1))*T + t + 1
        stopifnot(length(rows_t) == N)
        stopifnot(max(rows_t) <= N * T)
        time_fe[rows_t, t] <- rep(1, N)
    }

    return(list(cohort_fe=cohort_fe, time_fe=time_fe, assignments=assignments,
        inds=inds))
}

genAssignments <- function(N, R){
    # Make sure at least one observation in each cohort
    stopifnot(N >= R + 1)
    pass_condition <- FALSE
    while(!pass_condition){
        assignments <- rmultinom(n=1, size=N, prob=rep(1/(R + 1), R + 1))[, 1]
        pass_condition <- all(assignments >= 1)
    }

    stopifnot(sum(assignments) == N)
    stopifnot(all(assignments >= 1))
    stopifnot(length(assignments) == R + 1)
    return(assignments)
}

genTreatVarsSim <- function(n_treats, N, T, R, assignments, cohort_inds,
    N_UNTREATED, first_inds_test, d){
    # Treatment indicators
    treat_mat_long <- matrix(0, N*T, n_treats)
    treat_ind <- 0
    total_feats_added <- 0
    first_inds <- rep(as.integer(NA), R)

    stopifnot(length(cohort_inds) == R)
    stopifnot(length(assignments) == R + 1)
    stopifnot(min(cohort_inds[[1]]) == N_UNTREATED * T + 1)

    all_inds_so_far <- integer()

    for(r in 1:R){
        n_treats_r <- T - r

        stopifnot(length(cohort_inds[[r]]) == assignments[r + 1]*T)

        counter <- 0L

        for(t in (r + 1):T){
   
            treat_ind <- treat_ind + 1

            stopifnot(treat_ind <= n_treats)

            if(t == r + 1){
                first_inds[r] <- treat_ind

                stopifnot(treat_ind == first_inds_test[r])
            }

            first_ind_r <-  min(cohort_inds[[r]]) + t - (r + 1) + r

            stopifnot(first_ind_r >= min(cohort_inds[[r]]))

            last_ind_r <- first_ind_r + (assignments[r + 1] - 1)*T

            stopifnot(last_ind_r <= max(cohort_inds[[r]]))
            stopifnot((last_ind_r - first_ind_r)/T ==
                round((last_ind_r - first_ind_r)/T))
            stopifnot((last_ind_r - first_ind_r)/T == assignments[r + 1] - 1)

            r_inds <- seq(first_ind_r, last_ind_r, by=T)

            all_inds_so_far <- c(all_inds_so_far, r_inds)

            stopifnot(length(r_inds) == assignments[r + 1])
            stopifnot(all(treat_mat_long[, treat_ind] == 0))

            treat_mat_long[r_inds, treat_ind] <- 1
            total_feats_added <- total_feats_added + 1
            counter <- counter + 1
        }

        stopifnot(counter == n_treats_r)

        treat_inds_r <- first_inds[r]:(first_inds[r] + length((r + 1):T) - 1)

        stopifnot(all(treat_inds_r <= n_treats))
        stopifnot(length(treat_inds_r) == n_treats_r)
        stopifnot(all(colSums(treat_mat_long[, treat_inds_r, drop=FALSE]) ==
            assignments[r + 1]))

        if(r < R){
            stopifnot(last_ind_r == min(cohort_inds[[r + 1]]) - 1)
            stopifnot(last_ind_r == max(cohort_inds[[r]]))
        }
    }

    stopifnot(max(r_inds) == N*T)

    stopifnot(all(colSums(treat_mat_long) >= 1))
    stopifnot(total_feats_added == n_treats)

    stopifnot(all(!is.na(first_inds)))
    stopifnot(length(first_inds) == length(unique(first_inds)))
    stopifnot(first_inds[1] == 1)

    stopifnot(identical(first_inds, first_inds_test))

    return(list(treat_mat_long=treat_mat_long, first_inds=first_inds))
}

getActualCohortTes <- function(R, first_inds, treat_inds, coefs, num_treats){
    actual_cohort_tes <- rep(as.numeric(NA), R)

    for(r in 1:R){
        first_ind_r <- first_inds[r]
        if(r < R){
            last_ind_r <- first_inds[r + 1] - 1
        } else{
            last_ind_r <- num_treats
        }
        
        actual_cohort_tes[r] <- mean(coefs[treat_inds][first_ind_r:last_ind_r])
    }

    return(actual_cohort_tes)
}

# Function for concisely picking a sign randomly
rfunc <- function(n, prob){
    # sample(c(-1, 1), size=n, replace=TRUE)
    vec <- rbinom(n=n, size=1, prob=prob)
    vec[vec == 0] <- -1
    stopifnot(all(vec %in% c(-1, 1)))
    return(vec)
}

testGenRandomDataInputs <- function(beta, R, T, d, N, sig_eps_sq, sig_eps_c_sq){

    stopifnot(R <= T - 1)
    stopifnot(T >= 3)
    stopifnot(R >= 2)
    stopifnot(N >= R + 1)
    stopifnot(sig_eps_sq > 0)
    stopifnot(sig_eps_c_sq > 0)

    # Compute the number of treatment effects (common to both cases)
    num_treats <- getNumTreats(R=R, T=T)
  
    # --- Full design matrix with interactions ---
    # Expected number of columns:
    # If d > 0: p = R + (T - 1) + d + d*R + d*(T - 1) + num_treats + num_treats*d
    # If d == 0: p = R + (T - 1) + num_treats
    p_expected <- getP(R=R, T=T, d=d, num_treats=num_treats)

    if (length(beta) != p_expected) {
      stop(sprintf("For gen_ints = TRUE, length(beta) must be %d", p_expected))
    }

    return(list(num_treats=num_treats, p_expected=p_expected))
}

#' Generate Coefficient Vector for Data Generation
#'
#' This function generates a coefficient vector \code{beta} along with a sparse auxiliary vector
#' \code{theta} for simulation studies of the fused extended two-way fixed effects estimator. The
#' returned \code{beta} is formatted to align with the design matrix created by \code{genRandomData()},
#' and is a valid input for the \code{beta} argument of that function. The vector \code{theta} is sparse,
#' with nonzero entries occurring with probability \code{density} and scaled by \code{eff_size}.
#'
#' @param R Integer. The number of treated cohorts (treatment is assumed to start in periods 2 to \code{R + 1}).
#' @param T Integer. The total number of time periods.
#' @param d Integer. The number of time-invariant covariates. If \code{d > 0}, additional terms corresponding
#'   to covariate main effects and interactions are included in \code{beta}.
#' @param density Numeric in (0,1). The probability that any given entry in the initial sparse coefficient
#'   vector \code{theta} is nonzero.
#' @param eff_size Numeric. The magnitude used to scale nonzero entries in \code{theta}. Each nonzero entry is
#'   set to \code{eff_size} or \code{-eff_size} (with a 60 percent chance for a positive value).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{beta}}{A numeric vector representing the full coefficient vector after the inverse fusion transform.}
#'   \item{\code{theta}}{A numeric vector that is sparse, from which \code{beta} is derived.}
#' }
#'
#' @details
#' The length of \code{beta} is given by
#' \deqn{p = R + (T - 1) + d + dR + d(T - 1) + \mathit{num\_treats} + (\mathit{num\_treats} \times d)}{p = R + (T - 1) + d + dR + d(T - 1) + num_treats + (num_treats * d)},
#' where the number of treatment parameters is defined as
#' \deqn{\mathit{num\_treats} = T \times R - \frac{R(R+1)}{2}}{num_treats = T * R - R(R+1)/2}.
#'
#' The function operates in two steps:
#' \enumerate{
#'   \item It first creates a sparse vector \code{theta} of length \eqn{p}, with nonzero entries occurring
#'   with probability \code{density}. Nonzero entries are set to \code{eff_size} or \code{-eff_size} (with a 60\%
#'   chance for a positive value).
#'   \item The full coefficient vector \code{beta} is then computed by applying an inverse fusion transform to \code{theta}
#'   using internal routines (e.g., \code{genBackwardsInvFusionTransformMat()} and \code{genInvTwoWayFusionTransformMat()}).
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate coefficients for a panel with 5 treated cohorts, 30 time periods,
#'   # a density of 0.1 for nonzero coefficients, effect size 2, and 12 covariates.
#'   coefs <- genCoefs(R = 5, T = 30, density = 0.1, eff_size = 2, d = 12)
#'
#'   # Use the generated beta vector in a simulation study:
#'   simData <- genRandomData(
#'     N = 120, T = 30, R = 5, d = 12,
#'     sig_eps_sq = 5, sig_eps_c_sq = 5,
#'     beta = coefs$beta,
#'     seed = 123
#'   )
#' }
#'
genCoefsCore <- function(R, T, d, density, eff_size){

    # Check that T is a numeric scalar and at least 3.
  if (!is.numeric(T) || length(T) != 1 || T < 3) {
    stop("T must be a numeric value greater than or equal to 3")
  }
  
  # Check that R is a numeric scalar and at least 2.
  if (!is.numeric(R) || length(R) != 1 || R < 2) {
    stop("R must be a numeric value greater than or equal to 2 (currently there is only support for data sets with staggered adoptions, so at least two treated cohorts)")
  }
  
  # Check that R does not exceed T - 1.
  if (R > T - 1) {
    stop("R must be less than or equal to T - 1")
  }
  
  # Check that d is a numeric scalar and is non-negative.
  if (!is.numeric(d) || length(d) != 1 || d < 0) {
    stop("d must be a non-negative numeric value")
  }
  
  # Check that density is a numeric scalar strictly between 0 and 1.
  if (!is.numeric(density) || length(density) != 1 || density <= 0 || density >= 1) {
    stop("density must be numeric and strictly between 0 and 1")
  }
  
  # Check that eff_size is numeric.
  if (!is.numeric(eff_size) || length(eff_size) != 1) {
    stop("eff_size must be a numeric value")
  }

    stopifnot(R >= 2)
    stopifnot(T >= 3)
    stopifnot(R <= T - 1)

    num_treats <- getNumTreats(R=R, T=T)

    p <- getP(R=R, T=T, d=d, num_treats=num_treats)

    theta <- rep(0, p)

    # Make sure at least one feature is selected
    pass_condition <- FALSE
    while(!pass_condition){
        theta_inds <- which(as.logical(rbinom(n=p, size=1, prob=density)))
        pass_condition <- length(theta_inds > 0)
    }
    
    num_coefs <- length(theta_inds)
    # Generate signs of coefficients in transformed space, and bias away from
    # 0.5 (as described in paper)
    signs <- rfunc(num_coefs, prob=0.6)

    theta[theta_inds] <- eff_size*signs

    # Now we have coefficients that are sparse in the appropriate feature space.
    # The last step is to transform them to the original feature space. Since
    # theta = D %*% beta, beta = solve(D) %*% theta.
    beta <- rep(as.numeric(NA), p)

    beta[1:R] <- genBackwardsInvFusionTransformMat(R) %*% theta[1:R]

    stopifnot(all(is.na(beta[(R + 1):(R + T - 1)])))
    beta[(R + 1):(R + T - 1)] <- genBackwardsInvFusionTransformMat(T - 1) %*%
        theta[(R + 1):(R + T - 1)]

    if(d > 0){
        # Coefficients corresponding to X don't need to be transformed
        stopifnot(all(is.na(beta[(R + T - 1 + 1):(R + T - 1 + d)])))
        beta[(R + T - 1 + 1):(R + T - 1 + d)] <- theta[(R + T - 1 + 1):
            (R + T - 1 + d)]

        # Cohort-X interactions (one cohort at a time, with all interactions for
        # X. So R blocks of size d.)

        for(j in 1:d){
            first_ind_j <- R + T - 1 + d + j
            last_ind_j <- R + T - 1 + d + (R - 1)*d + j

            inds_j <- seq(first_ind_j, last_ind_j, by=d)

            stopifnot(length(inds_j) == R)
            stopifnot(all(is.na(beta[inds_j])))

            beta[inds_j] <- genBackwardsInvFusionTransformMat(R) %*% theta[inds_j]
        }

        stopifnot(all(!is.na(beta[1:(R + T - 1 + d + R*d)])))
        stopifnot(all(is.na(beta[(R + T - 1 + d + R*d + 1):p])))

        # Time-X interactions
        for(j in 1:d){
            first_ind_j <- R + T - 1 + d + R*d + j
            last_ind_j <- R + T - 1 + d + R*d + (T - 2)*d + j

            inds_j <- seq(first_ind_j, last_ind_j, by=d)
            stopifnot(length(inds_j) == T - 1)
            stopifnot(all(is.na(beta[inds_j])))

            beta[inds_j] <- genBackwardsInvFusionTransformMat(T - 1) %*% theta[inds_j]
        }

        stopifnot(all(!is.na(beta[1:(R + T - 1 + d + R*d + (T - 1)*d)])))
        stopifnot(all(is.na(beta[(R + T - 1 + d + R*d + (T - 1)*d + 1):p])))
    }

    # Base treatment effects: need to identify indices of first treatment
    # effect for each cohort
    first_inds <- getFirstInds(R=R, T=T)

    treat_inds <- getTreatInds(R=R, T=T, d=d, num_treats=num_treats)

    stopifnot(all(is.na(beta[treat_inds])))

    beta[treat_inds] <- genInvTwoWayFusionTransformMat(num_treats,
        first_inds, R) %*% theta[treat_inds]

    stopifnot(all(!is.na(beta[1:(R + T - 1 + d + R*d + (T - 1)*d +
        num_treats)])))

    if(d > 0){

        stopifnot(all(is.na(beta[(R + T - 1 + d + R*d + (T - 1)*d + num_treats +
            1):p])))

        # Treatment effect-X interactions
        for(j in 1:d){
            first_ind_j <- R + T - 1 + d + R*d + (T - 1)*d + num_treats + j
            last_ind_j <- R + T - 1 + d + R*d + (T - 1)*d + num_treats +
                (num_treats - 1)*d + j

            inds_j <- seq(first_ind_j, last_ind_j, by=d)

            stopifnot(length(inds_j) == num_treats)
            stopifnot(all(is.na(beta[inds_j])))

            beta[inds_j] <- genInvTwoWayFusionTransformMat(num_treats,
                first_inds, R) %*% theta[inds_j]
        }
    }

    stopifnot(all(!is.na(beta)))

    # Confirm beta satisfies input requirements of genRandomData() (make
    # up values for N, sig_eps_sq, and sig_eps_c_sq that meet requirements)

    testGenRandomDataInputs(
        beta=beta,
        R=R,
        T=T,
        d=d,
        N=R+1,
        sig_eps_sq=1,
        sig_eps_c_sq=1)

    return(list(beta=beta, theta=theta))
}


#' Generate Random Panel Data for FETWFE Simulations
#'
#' @description
#' Generates a random panel dataset for simulation studies of the fused extended two-way fixed
#' effects (FETWFE) estimator. The function creates a balanced panel with \eqn{N} units over \eqn{T}
#' time periods, assigns treatment status across \eqn{R} treated cohorts (with equal marginal
#' probabilities for treatment and non-treatment), and constructs a design matrix along with the
#' corresponding outcome. When \code{gen_ints = TRUE} the full design matrix is returned (including
#' interactions between covariates and fixed effects and treatment indicators). When 
#' \code{gen_ints = FALSE} the design matrix is generated in a simpler format (with no interactions)
#' as expected by \code{fetwfe()}. Moreover, the covariates are generated according to the
#' specified \code{distribution}: by default, covariates are drawn from a normal distribution;
#' if \code{distribution = "uniform"}, they are drawn uniformly from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' When \eqn{d = 0} (i.e. no covariates), no covariate-related columns or interactions are generated.
#'
#' @param N Integer. Number of units in the panel.
#' @param T Integer. Number of time periods.
#' @param R Integer. Number of treated cohorts (with treatment starting in periods 2 to T).
#' @param d Integer. Number of time-invariant covariates.
#' @param sig_eps_sq Numeric. Variance of the idiosyncratic (observation-level) noise.
#' @param sig_eps_c_sq Numeric. Variance of the unit-level random effects.
#' @param beta Numeric vector. Coefficient vector for data generation. Its required length depends on
#'   the value of \code{gen_ints}:
#'   \itemize{
#'     \item If \code{gen_ints = TRUE} and \code{d > 0}, the expected length is 
#'       \eqn{p = R + (T-1) + d + dR + d(T-1) + num\_treats + num\_treats \times d}, where 
#'       \eqn{num\_treats = T \times R - \frac{R(R+1)}{2}}.
#'     \item If \code{gen_ints = TRUE} and \code{d = 0}, the expected length is 
#'       \eqn{p = R + (T-1) + num\_treats}.
#'     \item If \code{gen_ints = FALSE}, the expected length is 
#'       \eqn{p = R + (T-1) + d + num\_treats}.
#'   }
#' @param seed (Optional) Integer. Seed for reproducibility.
#' @param gen_ints Logical. If \code{TRUE}, generate the full design matrix with interactions;
#'   if \code{FALSE} (the default), generate a design matrix without any interaction terms.
#' @param distribution Character. Distribution to generate covariates.
#'   Defaults to \code{"gaussian"}. If set to \code{"uniform"}, covariates are drawn uniformly
#'   from \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{pdf}{A dataframe containing generated data that can be passed to \code{fetwfe()}.}
#'   \item{X}{The design matrix. When \code{gen_ints = TRUE}, \eqn{X} has \eqn{p} columns with
#'     interactions; when \code{gen_ints = FALSE}, \eqn{X} has no interactions.}
#'   \item{y}{A numeric vector of length \eqn{N \times T} containing the generated responses.}
#'   \item{covs}{A character vector containing the names of the generated features (if \eqn{d > 0}),
#'          or simply an empty vector (if \eqn{d = 0})}
#'   \item{coefs}{The coefficient vector \eqn{\beta} used for data generation.}
#'   \item{first_inds}{A vector of indices indicating the first treatment effect for each treated cohort.}
#'   \item{N_UNTREATED}{The number of never-treated units.}
#'   \item{assignments}{A vector of counts (of length \eqn{R+1}) indicating how many units fall into
#'         the never-treated group and each of the \eqn{R} treated cohorts.}
#'   \item{indep_counts}{Independent cohort assignments (for auxiliary purposes).}
#'   \item{p}{The number of columns in the design matrix \eqn{X}.}
#'   \item{N}{Number of units.}
#'   \item{T}{Number of time periods.}
#'   \item{R}{Number of treated cohorts.}
#'   \item{d}{Number of covariates.}
#'   \item{sig_eps_sq}{The idiosyncratic noise variance.}
#'   \item{sig_eps_c_sq}{The unit-level noise variance.}
#' }
#'
#' @details
#' When \code{gen_ints = TRUE}, the function constructs the design matrix by first generating
#' base fixed effects and a long-format covariate matrix (via \code{generateBaseEffects()}), then
#' appending interactions between the covariates and cohort/time fixed effects (via \code{generateFEInts()})
#' and finally treatment indicator columns and treatment-covariate interactions
#' (via \code{genTreatVarsSim()} and \code{genTreatInts()}). When \code{gen_ints = FALSE},
#' the design matrix consists only of the base fixed effects, covariates, and treatment indicators.
#'
#' The argument \code{distribution} controls the generation of covariates. For
#' \code{"gaussian"}, covariates are drawn from \code{rnorm}; for \code{"uniform"},
#' they are drawn from \code{runif} on the interval \eqn{[-\sqrt{3}, \sqrt{3}]}.
#'
#' When \eqn{d = 0} (i.e. no covariates), the function omits any covariate-related columns
#' and their interactions.
#'
#' @examples
#' \dontrun{
#' # Full design with interactions (default behavior, with gaussian covariates):
#' N <- 120; T <- 30; R <- 5; d <- 12; sig_eps_sq <- 5; sig_eps_c_sq <- 5
#' num_treats <- getNumTreats(R=R, T=T)
#' p_int <- R + (T - 1) + d + d * R + d * (T - 1) + num_treats + num_treats * d
#' beta_int <- rnorm(p_int)
#' sim_int <- genRandomData(N, T, R, d, sig_eps_sq, sig_eps_c_sq, beta_int,
#'                           seed = 123, gen_ints = TRUE, distribution = "gaussian")
#'
#' # Simple design without interactions using uniform covariates returned:
#' sim_no_int <- genRandomData(N, T, R, d, sig_eps_sq, sig_eps_c_sq, beta_int,
#'                             seed = 123, gen_ints = FALSE, distribution = "uniform")
#'
#' # When d = 0, no covariate or interaction terms are generated.
#' p_no_cov <- R + (T - 1) + num_treats
#' beta_no_cov <- rnorm(p_no_cov)
#' sim_no_cov <- genRandomData(N, T, R, 0, sig_eps_sq, sig_eps_c_sq, beta_no_cov,
#'                             seed = 123, gen_ints = TRUE, distribution = "gaussian")
#' }
#'
simulateDataCore <- function(N, T, R, d, sig_eps_sq, sig_eps_c_sq, beta, seed = NULL, 
                          gen_ints = FALSE, distribution = "gaussian") {
    if (!is.null(seed)) set.seed(seed)

    res <- testGenRandomDataInputs(
        beta=beta,
        R=R,
        T=T,
        d=d,
        N=N,
        sig_eps_sq=sig_eps_sq,
        sig_eps_c_sq=sig_eps_c_sq)

    num_treats <- res$num_treats
    p_expected <- res$p_expected

    rm(res)

    # Generate base effects and covariates (using specified distribution)
    res_base <- generateBaseEffects(
        N=N,
        d=d,
        T=T,
        R=R,
        distribution = distribution
        )
    cohort_fe <- res_base$cohort_fe
    time_fe <- res_base$time_fe
    X_long <- res_base$X_long
    assignments <- res_base$assignments
    cohort_inds <- res_base$cohort_inds

    stopifnot(ncol(cohort_fe) == R)

    # Base matrix: cohort FE, time FE, and covariates (if any)
    X_base <- if (d > 0) {
      cbind(cohort_fe, time_fe, X_long)
    } else {
      cbind(cohort_fe, time_fe)
    }

    indep_assignments <- genAssignments(N, R)

    if (d > 0) {
      res_ints <- generateFEInts(X_long, cohort_fe, time_fe, N, T, R, d)
      X_ints1 <- cbind(X_base, res_ints$X_long_cohort, res_ints$X_long_time)
    } else {
      X_ints1 <- X_base
    }

    first_inds_test <- getFirstInds(R=R, T=T)
    res_treat <- genTreatVarsSim(num_treats, N, T, R, assignments, cohort_inds,
                                 N_UNTREATED = assignments[1],
                                 first_inds_test = first_inds_test, d = d)
    treat_mat_long <- res_treat$treat_mat_long
    first_inds <- res_treat$first_inds

    X_ints2 <- cbind(X_ints1, treat_mat_long)

    if (d > 0) {
        stopifnot(ncol(cohort_fe) == R)
        X_long_treat <- genTreatInts(
            treat_mat_long=treat_mat_long,
            X_long=X_long,
            n_treats=num_treats,
            cohort_fe=cohort_fe,
            N=N,
            T=T,
            R=R,
            d=d,
            N_UNTREATED = assignments[1]
            )

        X_final <- cbind(X_ints2, X_long_treat)
    } else {
      X_final <- X_ints2
    }

    if (ncol(X_final) != p_expected) {
      stop("Constructed design matrix with interactions has incorrect number of columns.")
    }

    unit_res <- rnorm(N, mean = 0, sd = sqrt(sig_eps_c_sq))
    y <- X_final %*% beta + rep(unit_res, each = T) + 
      rnorm(N * T, mean = 0, sd = sqrt(sig_eps_sq))
    y <- y - mean(y)

    if (gen_ints){
        X_ret <- X_final
    } else{
        # Return X with no interactions
        # Expected number of columns: p = R + (T - 1) + d + num_treats
        p_expected <- R + (T - 1) + d + num_treats

        X_ret <- cbind(X_base, treat_mat_long)
        if (ncol(X_ret) != p_expected) {
          stop("Constructed design matrix without interactions has incorrect number of columns.")
        }
    }

    # Prepare dataframe for `fetwfe()`

    # We know that when gen_ints = FALSE, the design matrix X is:
    # X = [cohort_fe, time_fe, X_long, treat_mat_long]
    # The base part (cohort_fe, time_fe, X_long) has (R + (T-1) + d) columns.
    base_cols <- R + (T - 1) + d

    # The treatment dummy block is in columns (base_cols + 1) : (base_cols + num_treats)
    treat_dummy <- cbind(X_base, treat_mat_long)[, (base_cols + 1):(base_cols + num_treats), drop = FALSE]
    # For each row, the observed treatment indicator is 1 if any entry in treat_dummy is 1.
    treatment_vec <- as.integer(apply(treat_dummy, 1, function(x) { any(x == 1) }))

    # Extract the covariate columns from the base part: they are the last d columns of the base part.
    if(d > 0){
        cov_cols <- seq(from = (R + (T - 1) + 1), to = (R + (T - 1) + d))
        covariates <- cbind(X_base, treat_mat_long)[, cov_cols, drop = FALSE]
    }
    

    # Construct a data frame with the panel structure.
    df_panel <- data.frame(
    time = rep(1:T, times = N),
    unit = rep(sprintf("unit%02d", 1:N), each = T),
    treatment = treatment_vec,
    y = as.numeric(y)
    )

    cov_names <- c()

    if(d > 0){
        # Add covariate columns with names "cov1", "cov2", ...
        for (j in seq_len(d)) {
            cov_name_j <- paste0("cov", j)
            df_panel[[cov_name_j]] <- covariates[, j]
            cov_names <- c(cov_names, cov_name_j)
        }
    }

    stopifnot(length(cov_names) == d)

    # Ensure that time is integer and unit is character.
    df_panel$time <- as.integer(df_panel$time)
    df_panel$unit <- as.character(df_panel$unit)
    df_panel$treatment <- as.integer(df_panel$treatment)

    # confirm that outputs satisfy input requirements of fetwfe()
    # (Plug in default values for arguments not generated here)
    checkFetwfeInputs(
        pdata=df_panel,
        time_var="time",
        unit_var="unit",
        treatment="treatment",
        response="y",
        covs=cov_names,
        indep_counts=indep_assignments,
        sig_eps_sq=sig_eps_sq,
        sig_eps_c_sq=sig_eps_c_sq,
        lambda.max=NA,
        lambda.min=NA,
        nlambda=100,
        q=0.5,
        verbose=FALSE,
        alpha=0.05,
        add_ridge=FALSE
    )
    
    return(list(
      pdata = df_panel,
      X = X_ret,
      y = y,
      covs = cov_names,
      time_var = "time",
      unit_var = "unit",
      treatment = "treatment",
      response = "y",
      coefs = beta,
      first_inds = first_inds,
      N_UNTREATED = assignments[1],
      assignments = assignments,
      indep_counts = indep_assignments,
      p = p_expected,
      N = N,
      T = T,
      R = R,
      d = d,
      sig_eps_sq = sig_eps_sq,
      sig_eps_c_sq = sig_eps_c_sq
    ))
}