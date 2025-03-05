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




