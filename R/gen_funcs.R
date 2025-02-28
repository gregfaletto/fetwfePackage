generateBaseEffects <- function(N, d, T, R){
    ret <- genCohortTimeFE(N, T, R, d)

    # Generate time-invariant covariates
    X <- matrix(rnorm(N*d), N, d)

    # Put into matrix form (will have T observations from first unit, then T
    # from next, and so on)
    X_long <- X[rep(1:N, each=T), ]

    return(list(cohort_fe=ret$cohort_fe, time_fe=ret$time_fe, X_long=X_long,
        assignments=ret$assignments, cohort_inds=ret$inds))
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

getFirstInds <- function(n_treats, R, T){
    # Let's identify the indices of the first treatment effects for each cohort.
    # The first one is index 1, then the second one is index (T - 1) + 1 = T,
    # then the third one is (T - 1) + (T - 2) + 1 = 2*T - 2. In general, for
    # r > 1 the rth one will occur at index
    #
    # (T - 1) + (T - 2) + ... + (T - (r - 1)) + 1
    # = 1 + (r - 1)*(T - 1 + T - r + 1)/2
    # = 1 + (r - 1)*(2*T - r)/2.
    #
    # (Looks like the formula works for r = 1 too.)
    f_inds <- integer(R)

    for(r in 1:R){
        f_inds[r] <- 1 + (r - 1)*(2*T - r)/2
    }
    stopifnot(all(f_inds <= n_treats))
    stopifnot(f_inds[1] == 1)
    # Last cohort has T - R treatment effects to estimate. So last first_ind
    # should be at position num_treats - (T - R) + 1 = num_treats - T + R + 1.
    stopifnot(f_inds[R] == n_treats - T + R + 1)

    return(f_inds)
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
        stopifnot(all(colSums(treat_mat_long[, treat_inds_r]) ==
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