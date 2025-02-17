#' @import glmnet

prepXints <- function(
    data,
    time_var,
    unit_var,
    treatment,
    covs,
    response,
    verbose=FALSE
    ){

    # Check inputs
    stopifnot(is.data.frame(data))
    stopifnot(nrow(data) >= 4) # bare minimum, 2 units at 2 times

    stopifnot(is.character(time_var))
    stopifnot(length(time_var) == 1)
    stopifnot(time_var %in% colnames(data))
    stopifnot(is.integer(data[[time_var]]))

    stopifnot(is.character(unit_var))
    stopifnot(length(unit_var) == 1)
    stopifnot(unit_var %in% colnames(data))
    stopifnot(is.character(data[[unit_var]]))

    stopifnot(is.character(treatment))
    stopifnot(length(treatment) == 1)
    stopifnot(treatment %in% colnames(data))
    stopifnot(is.integer(data[[treatment]]))
    stopifnot(all(data[, treatment] %in% c(0, 1)))

    stopifnot(is.character(covs))
    if(length(covs) < 1){
        stop("Must include at least one covariate.")
    }
    stopifnot(all(covs %in% colnames(data)))
    for(cov in covs){
        stopifnot(is.numeric(data[[cov]]) | is.integer(data[[cov]]))
    }

    stopifnot(is.character(response))
    stopifnot(length(response) == 1)
    stopifnot(response %in% colnames(data))
    stopifnot(is.numeric(data[[response]]) | is.integer(data[[response]]))

    stopifnot(is.logical(verbose))
    stopifnot(length(verbose) == 1)

    # Identify cohorts, units, and times; remove units that were treated in
    # first period; remove treatment variable after identifying cohorts
    ret <- idCohorts(
        df=data,
        time_var=time_var,
        unit_var=unit_var,
        treat_var=treatment,
        covs=covs
        )

    data <- ret$df
    cohorts <- ret$cohorts
    units <- ret$units
    times <- ret$times

    rm(ret)

    N <- length(units)
    T <- length(times)
    R <- length(cohorts)

    stopifnot(R <= T - 1)

    N_treated <- length(unlist(cohorts))

    stopifnot(N_treated <= N)

    N_UNTREATED <- N - N_treated

    # Get in-sample counts of units in each cohort
    in_sample_counts <- rep(as.numeric(NA), R + 1)

    for(r in 1:R){
        in_sample_counts[r + 1] <- length(cohorts[[r]])
    }

    in_sample_counts[1] <- N - N_treated
    
    stopifnot(all(in_sample_counts == round(in_sample_counts)))
    in_sample_counts <- as.integer(in_sample_counts)
    stopifnot(all(!is.na(in_sample_counts)))
    stopifnot(sum(in_sample_counts) == N)



    if(in_sample_counts[1] == 0){
        warning("No untreated units found. It will not be possible to estimate treatment effects on this data set.")
    }

    names(in_sample_counts) <- c("Never_treated", names(cohorts))

    # Replace time-varying covariates with value in first (pre-treated) period;
    # remove covariates with any missing values in this period
    ret <- processCovs(
        df=data,
        units=units,
        unit_var=unit_var,
        times=times,
        time_var=time_var,
        covs=covs,
        resp_var=response,
        T=T,
        verbose=verbose
        )

    data <- ret$df
    covs <- ret$covs

    rm(ret)

    d <- length(covs)

    stopifnot(all(!is.na(data)))
    stopifnot(nrow(data) == N*T)

    # Add cohort dummies, time dummies, and treatment variables. (Names of
    # cohorts must be same as times of treatment.) Also, center response.
    ret <- addDummies(
        df=data,
        cohorts=cohorts,
        times=times,
        N=N,
        T=T,
        unit_var=unit_var,
        time_var=time_var,
        resp_var=response,
        n_cohorts=R
        )

    y <- ret$y # response
    cohort_treat_names <- ret$cohort_treat_names # List of names of treatment
    # dummies for each cohort
    time_var_names <- ret$time_var_names # Names of time dummies
    cohort_vars <- ret$cohort_vars # Names of cohort dummies
    first_inds <- ret$first_inds # Among blocks of treatment variables, indices
    # of first treatment effects corresponding to each cohort
    cohort_var_mat <- ret$cohort_var_mat # Cohort dummies; R columns total
    time_var_mat <- ret$time_var_mat # Time dummies; T - 1 columns total
    treat_var_mat <- ret$treat_var_mat # Treatment dummies

    num_treats <- ncol(treat_var_mat)

    stopifnot(num_treats == length(unlist(cohort_treat_names)))
    stopifnot(all(covs %in% colnames(data)))

    stopifnot(is.numeric(ncol(treat_var_mat)) | is.integer(ncol(treat_var_mat)))
    stopifnot(ncol(treat_var_mat) >= 1)

    covariate_mat <- as.matrix(data[, covs]) # Covariates

    rm(ret)
    rm(data)

    #### Create matrix with all interactions
    p <- R + T - 1 + d * (1 + R + T - 1) + (d + 1) * num_treats

    X_ints <- genXintsData(
        cohort_fe=cohort_var_mat,
        time_fe=time_var_mat,
        X_long=covariate_mat,
        treat_mat_long=treat_var_mat,
        N=N,
        R=R,
        T=T,
        d=d,
        N_UNTREATED=N_UNTREATED,
        p=p
        )

    return(list(
        X_ints=X_ints,
        y=y,
        N=N,
        T=T,
        d=d,
        p=p,
        in_sample_counts=in_sample_counts,
        num_treats=num_treats,
        first_inds=first_inds)
    )
}


fetwfe_core <- function(
    X_ints,
    y,
    in_sample_counts,
    N,
    T,
    d,
    p,
    num_treats,
    first_inds,
    indep_counts=NA,
    sig_eps_sq=NA,
    sig_eps_c_sq=NA,
    lambda.max=NA,
    lambda.min=NA,
    nlambda=100,
    q=0.5,
    verbose=FALSE,
    alpha=0.05
    ){

    #
    #
    # Check inputs
    #
    #

    stopifnot(N >= 2) # bare minimum, 2 units at 2 times

    stopifnot(T >= 2) # bare minimum, 2 units at 2 times

    if(any(!is.na(sig_eps_sq))){
        stopifnot(is.numeric(sig_eps_sq) | is.integer(sig_eps_sq))
        stopifnot(length(sig_eps_sq) == 1)
        stopifnot(sig_eps_sq >= 0)
    }

    if(any(!is.na(sig_eps_c_sq))){
        stopifnot(is.numeric(sig_eps_c_sq) | is.integer(sig_eps_c_sq))
        stopifnot(length(sig_eps_c_sq) == 1)
        stopifnot(sig_eps_c_sq >= 0)
    }

    stopifnot(sum(in_sample_counts) == N)
    stopifnot(all(in_sample_counts >= 0))
    if(in_sample_counts[1] == 0){
        stop("No never-treated units detected in data to fit model; estimating treatment effects is not possible")
    }
    if(length(names(in_sample_counts)) != length(in_sample_counts)){
        stop("in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)")
    }

    if(length(names(in_sample_counts)) != length(unique(names(in_sample_counts)))){
        stop("in_sample_counts must have all unique named entries (with names corresponding to the names of each cohort)")
    }

    R <- length(in_sample_counts) - 1
    stopifnot(R >= 1)
    stopifnot(R <= T - 1)

    indep_count_data_available <- FALSE
    if(any(!is.na(indep_counts))){
        if(sum(indep_counts) != N){
            stop("Number of units in independent cohort count data does not equal number of units in data to be used to fit model.")
        }
        if(length(indep_counts) != length(in_sample_counts)){
            stop("Number of counts in independent counts does not match number of cohorts in data to be used to fit model.")
        }
        if(any(indep_counts <= 0)){
            stop("At least one cohort in the independent count data has 0 members")
        }
        indep_count_data_available <- TRUE
    }

    if(any(!is.na(lambda.max))){
        stopifnot(is.numeric(lambda.max) | is.integer(lambda.max))
        stopifnot(length(lambda.max) == 1)
        stopifnot(lambda.max >= 0)
    }

    if(any(!is.na(lambda.min))){
        stopifnot(is.numeric(lambda.min) | is.integer(lambda.min))
        stopifnot(length(lambda.min) == 1)
        stopifnot(lambda.min >= 0)
        if(any(!is.na(lambda.max))){
            stopifnot(lambda.max >= lambda.min)
        }
    }

    stopifnot(is.numeric(q) | is.integer(q))
    stopifnot(length(q) == 1)
    stopifnot(q > 0)
    stopifnot(q <= 2)

    stopifnot(is.logical(verbose))
    stopifnot(length(verbose) == 1)

    stopifnot(is.numeric(alpha))
    stopifnot(length(alpha) == 1)
    stopifnot(alpha > 0)
    stopifnot(alpha < 1)

    #
    #
    # Step 1: change coordinates of data so that regular bridge regression
    # penalty applied to transformed dataframe results in FETWFE penalty applied
    # to original data
    #
    #

    if(verbose){
        message("Transforming matrix...")
    }
    
    # Transform matrix (change of coordinates so that fitting regular bridge
    # regression results in FETWFE fusion penalties)
    X_mod <- transformXintDataApp(
        X_ints,
        N=N,
        T=T,
        R=R,
        d=d,
        num_treats=num_treats,
        first_inds=first_inds
        )

    #
    #
    # Step 2: get (known or estimated) covariance matrix within observations
    # (due to umit-level random effects) and pre-multiply X and y by
    # inverse square root matrix
    #
    #

    if(verbose){
        message("Getting omega sqrt inverse estimate...")
        t0 <- Sys.time()
    }

    if(is.na(sig_eps_sq) | is.na(sig_eps_c_sq)){
        # Get omega_sqrt_inv matrix to multiply y and X_mod by on the left
        omega_res <- estOmegaSqrtInv(
            y,
            X_ints,
            N=N,
            T=T,
            p=p
            )

        sig_eps_sq <- omega_res$sig_eps_sq
        sig_eps_c_sq <- omega_res$sig_eps_c_sq

        rm(omega_res)

        if(verbose){
            message("Done! Time to estimate noise variances:")
            message(Sys.time() - t0)
            t0 <- Sys.time()
        }
    }

    stopifnot(!is.na(sig_eps_sq) & !is.na(sig_eps_c_sq))

    Omega <- diag(rep(sig_eps_sq, T)) + matrix(sig_eps_c_sq, T, T)

    Omega_sqrt_inv <- expm::sqrtm(solve(Omega))

    if(verbose){
        message("Time to get sqrt inverse matrix:")
        message(Sys.time() - t0)
    }

    y_final <- kronecker(diag(N), sqrt(sig_eps_sq)*Omega_sqrt_inv) %*% y
    X_final <- kronecker(diag(N), sqrt(sig_eps_sq)*Omega_sqrt_inv) %*% X_mod

    #
    #
    # Step 3: get cohort-specific sample proportions (estimated treatment
    # probabilities)
    #
    #

    cohort_probs <- in_sample_counts[2:(R + 1)]/
        sum(in_sample_counts[2:(R + 1)])

    stopifnot(all(!is.na(cohort_probs)))
    stopifnot(all(cohort_probs >= 0))
    stopifnot(all(cohort_probs <= 1))
    stopifnot(length(cohort_probs) == R)
    stopifnot(abs(sum(cohort_probs) - 1) < 10^(-6))
   
    cohort_probs_overall <- in_sample_counts[2:(R + 1)]/N

    stopifnot(abs(1 - sum(cohort_probs_overall) - in_sample_counts[1]/N) <
        10^(-6))

    if(indep_count_data_available){
        indep_cohort_probs <- indep_counts[2:(R + 1)]/
            sum(indep_counts[2:(R + 1)])

        stopifnot(all(!is.na(indep_cohort_probs)))
        stopifnot(all(indep_cohort_probs >= 0))
        stopifnot(all(indep_cohort_probs <= 1))
        stopifnot(length(indep_cohort_probs) == R)
        stopifnot(abs(sum(indep_cohort_probs) - 1) < 10^(-6))

        indep_cohort_probs_overall <- indep_counts[2:(R + 1)]/N

        stopifnot(abs(1 - sum(
            indep_cohort_probs_overall) - indep_counts[1] / N) < 10^(-6))
    } else{
        indep_cohort_probs <- NA
        indep_cohort_probs_overall <- NA
    }

    #
    #
    # Step 4: estimate bridge regression and extract fitted coefficients
    #
    #

    # Estimate bridge regression
    if(verbose){
        message("Estimating bridge regression...")
        t0 <- Sys.time()
    }

    if(!is.na(lambda.max) & !is.na(lambda.min)){
        fit <- grpreg::gBridge(
            X=X_final,
            y=y_final,
            gamma=q,
            lambda.max=lambda.max,
            lambda.min=lambda.min,
            nlambda=nlambda
            )
    } else if(!is.na(lambda.max)){
        fit <- grpreg::gBridge(
            X=X_final,
            y=y_final,
            gamma=q,
            lambda.max=lambda.max,
            nlambda=nlambda
            )
    } else if(!is.na(lambda.min)){
        fit <- grpreg::gBridge(
            X=X_final,
            y=y_final,
            gamma=q,
            lambda.min=lambda.min,
            nlambda=nlambda
            )
    } else{
        fit <- grpreg::gBridge(
            X=X_final,
            y=y_final,
            gamma=q,
            nlambda=nlambda
            )
    }
    
    if(verbose){
        message("Done! Time for estimation:")
        message(Sys.time() - t0)
    }

    # For diagnostics later, store largest and smallest lambda, as well as 
    # corresponding smallest and largest model sizes, to return.
    lambda.max <- max(fit$lambda)
    lambda.max_model_size <- sum(fit$beta[, ncol(beta)] != 0)

    lambda.min <- min(fit$lambda)
    lambda.min_model_size <- sum(fit$beta[, 1] != 0)

    # Select a single set of fitted coefficients by using BIC to choose among
    # the penalties that were fitted
    res <- getBetaBIC(
        fit,
        N=N,
        T=T,
        p=p,
        X_mod=X_mod,
        y=y
        )

    theta_hat <- res$theta_hat
    lambda_star_ind <- res$lambda_star_ind
    lambda_star_model_size <- res$lambda_star_model_size

    lambda_star <- fit$lambda[lambda_star_ind]

    c_names <- names(in_sample_counts)[2:(R + 1)]
    stopifnot(length(c_names) == R)

    # Handle edge case where no features are selected
    if(lambda_star_model_size == 0){
        if(verbose){
            message("No features selected; all treatment effects estimated to be 0.")
        }
        
        if(q < 1){
            ret_se <- 0
        } else{
            ret_se <- NA
        }

        catt_df_to_ret <- data.frame(c_names, rep(0, R), rep(ret_se, R),
                rep(ret_se, R), rep(ret_se, R))

        colnames(cat_df_to_ret) <- c("Cohort", "Estimated TE", "SE",
            "ConfIntLow", "ConfIntHigh")
        return(list(
            in_sample_att_hat=0,
            in_sample_att_se=ret_se,
            in_sample_att_se_no_prob=ret_se,
            indep_att_hat=0,
            indep_att_se=ret_se,
            # indep_att_se_no_prob=indep_att_se_no_prob,
            catt_hats=rep(0, R),
            catt_ses=rep(ret_se, R),
            catt_df=catt_df_to_ret,
            theta_hat=theta_hat,
            beta_hat=theta_hat,
            treat_inds=treat_inds,
            treat_int_inds=treat_int_inds,
            cohort_probs=cohort_probs,
            indep_cohort_probs=indep_cohort_probs,
            sig_eps_sq=sig_eps_sq,
            sig_eps_c_sq=sig_eps_c_sq,
            lambda.max=lambda.max,
            lambda.max_model_size=lambda.max_model_size,
            lambda.min=lambda.min,
            lambda.min_model_size=lambda.min_model_size,
            lambda_star=lambda_star,
            lambda_star_model_size=lambda_star_model_size,
            X_ints=X_ints,
            y=y,
            X_final=X_final,
            y_final=y_final,
            N=N,
            T=T,
            R=R,
            d=d,
            p=p
            )
        )
    }

    # intercept
    eta_hat <- theta_hat[1]

    # estimated coefficients
    theta_hat <- theta_hat[2:(p + 1)]

    # Indices of selected features in transformed feature space
    sel_feat_inds <- which(theta_hat != 0)

    # Indices corresponding to base treatment effects
    treat_inds <- (R + T - 1 + d + R*d + (T - 1)*d + 1):
        (R + T - 1 + d + R*d + (T - 1)*d + num_treats)

    stopifnot(length(treat_inds) == num_treats)
    stopifnot(max(treat_inds) + 1 <= p)
    stopifnot(max(treat_inds) == R + T - 1 + d + R*d + (T - 1)*d + num_treats)

    treat_int_inds <- (max(treat_inds) + 1):p

    stopifnot(length(treat_int_inds) == num_treats * d)

    # treat_int_ind_first <- R + T - 1 + d + R*d + (T - 1)*d + num_treats + 1
    # stopifnot(treat_int_ind_first == max(treat_inds) + 1)

    sel_treat_inds <- sel_feat_inds[sel_feat_inds %in% treat_inds]

    stopifnot(length(sel_treat_inds) == length(unique(sel_treat_inds)))
    stopifnot(length(sel_treat_inds) <= length(sel_feat_inds))
    stopifnot(length(sel_treat_inds) <= length(treat_inds))
    stopifnot(is.integer(sel_treat_inds) | is.numeric(sel_treat_inds))

    sel_treat_inds_shifted <- sel_treat_inds - (R + T - 1 + d + R*d + (T - 1)*d)

    stopifnot(all(sel_treat_inds_shifted %in% 1:num_treats))
    stopifnot(all(sel_treat_inds_shifted >= 1))
    stopifnot(all(sel_treat_inds_shifted <= num_treats))

    # Handle edgge case where no treatment features selected
    if(length(sel_treat_inds_shifted) == 0){
        if(verbose){
            message("No features selected; all treatment effects estimated to be 0.")
        }
        
        if(q < 1){
            ret_se <- 0
        } else{
            ret_se <- NA
        }

        catt_df_to_ret <- data.frame(c_names, rep(0, R), rep(ret_se, R),
                rep(ret_se, R), rep(ret_se, R))

        colnames(catt_df_to_ret) <- c("Cohort", "Estimated TE", "SE",
            "ConfIntLow", "ConfIntHigh")
        return(list(
            in_sample_att_hat=0,
            in_sample_att_se=ret_se,
            in_sample_att_se_no_prob=ret_se,
            indep_att_hat=0,
            indep_att_se=ret_se,
            # indep_att_se_no_prob=indep_att_se_no_prob,
            catt_hats=rep(0, R),
            catt_ses=rep(ret_se, R),
            catt_df=catt_df_to_ret,
            theta_hat=theta_hat,
            beta_hat=theta_hat,
            treat_inds=treat_inds,
            treat_int_inds=treat_int_inds,
            cohort_probs=cohort_probs,
            indep_cohort_probs=indep_cohort_probs,
            sig_eps_sq=sig_eps_sq,
            sig_eps_c_sq=sig_eps_c_sq,
            lambda.max=lambda.max,
            lambda.max_model_size=lambda.max_model_size,
            lambda.min=lambda.min,
            lambda.min_model_size=lambda.min_model_size,
            lambda_star=lambda_star,
            lambda_star_model_size=lambda_star_model_size,
            X_ints=X_ints,
            y=y,
            X_final=X_final,
            y_final=y_final,
            N=N,
            T=T,
            R=R,
            d=d,
            p=p
            )
        )
    }

    #
    #
    # Step 5: transform estimated coefficients back to original feature
    # space
    #
    #

    beta_hat <- untransformCoefsDataApp(
        theta_hat,
        first_inds,
        T=T,
        R=R,
        p=p,
        d=d,
        num_treats=num_treats
        )

    # TODO(gregfaletto): add labels to these estimates to make them easier to 
    # interpret

    # Get actual estimated treatment effects and standard errors
    tes <- beta_hat[treat_inds]

    stopifnot(length(tes) == num_treats)
    stopifnot(all(sel_treat_inds_shifted %in% 1:num_treats))

    stopifnot(length(tes) == num_treats)
    stopifnot(all(theta_hat[treat_inds][sel_treat_inds_shifted] != 0))
    stopifnot(all(theta_hat[treat_inds][setdiff(1:num_treats,
        sel_treat_inds_shifted)] == 0))

    stopifnot(length(first_inds) == R)
    stopifnot(max(first_inds) <= num_treats)

    stopifnot(length(sel_feat_inds) > 0)
    stopifnot(length(sel_treat_inds_shifted) > 0)

    #
    #
    # Step 6: calculate cohort-specific treatment effects and standard
    # errors
    #
    #

    res <- getCohortATTsFinal(
        X_final,
        sel_feat_inds,
        treat_inds,
        num_treats,
        first_inds,
        sel_treat_inds_shifted,
        c_names=c_names,
        tes=tes,
        sig_eps_sq=sig_eps_sq,
        R=R,
        N=N,
        T=T,
        fused=TRUE,
        calc_ses = q < 1,
        alpha=alpha
        )

    cohort_te_df <- res$cohort_te_df
    cohort_tes <- res$cohort_tes
    cohort_te_ses <- res$cohort_te_ses
    psi_mat <- res$psi_mat
    gram_inv <- res$gram_inv
    d_inv_treat_sel <- res$d_inv_treat_sel
    calc_ses <- res$calc_ses

    rm(res)

    if(calc_ses){
        stopifnot(nrow(d_inv_treat_sel) == num_treats)
        stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
    }

    

    #
    #
    # Step 7: calculate overall average treatment effect on treated units
    #
    #

    # Get overal estimated ATT!

    in_sample_te_results <- getTeResults2(
        sig_eps_sq=sig_eps_sq,
        N=N,
        T=T,
        R=R,
        num_treats=num_treats,
        cohort_tes=cohort_tes,
        cohort_probs=cohort_probs,
        psi_mat=psi_mat,
        gram_inv=gram_inv,
        sel_treat_inds_shifted=sel_treat_inds_shifted,
        tes=tes,
        d_inv_treat_sel=d_inv_treat_sel,
        cohort_probs_overall=cohort_probs_overall,
        first_inds=first_inds,
        theta_hat_treat_sel=theta_hat[sel_treat_inds],
        calc_ses=calc_ses,
        indep_probs=FALSE
        )

    in_sample_att_hat <- in_sample_te_results$att_hat
    in_sample_att_se <- in_sample_te_results$att_te_se
    in_sample_att_se_no_prob <- in_sample_te_results$att_te_se_no_prob

    if(indep_count_data_available){
        indep_te_results <- getTeResults2(
            sig_eps_sq=sig_eps_sq,
            N=N,
            T=T,
            R=R,
            num_treats=num_treats,
            cohort_tes=cohort_tes,
            cohort_probs=indep_cohort_probs,
            psi_mat=psi_mat,
            gram_inv=gram_inv,
            sel_treat_inds_shifted=sel_treat_inds_shifted,
            tes=tes,
            d_inv_treat_sel=d_inv_treat_sel,
            cohort_probs_overall=indep_cohort_probs_overall,
            first_inds=first_inds,
            theta_hat_treat_sel=theta_hat[sel_treat_inds],
            calc_ses=calc_ses,
            indep_probs=TRUE
            )
        indep_att_hat <- indep_te_results$att_hat
        indep_att_se <- indep_te_results$att_te_se
        indep_att_se_no_prob <- indep_te_results$att_te_se_no_prob
    } else{
        indep_att_hat <- NA
        indep_att_se <- NA
        indep_att_se_no_prob <- NA
    }



    return(list(
        in_sample_att_hat=in_sample_att_hat,
        in_sample_att_se=in_sample_att_se,
        in_sample_att_se_no_prob=in_sample_att_se_no_prob,
        indep_att_hat=indep_att_hat,
        indep_att_se=indep_att_se,
        # indep_att_se_no_prob=indep_att_se_no_prob,
        catt_hats=cohort_tes,
        catt_ses=cohort_te_ses,
        catt_df=cohort_te_df,
        theta_hat=theta_hat,
        beta_hat=beta_hat,
        treat_inds=treat_inds,
        treat_int_inds=treat_int_inds,
        cohort_probs=cohort_probs,
        indep_cohort_probs=indep_cohort_probs,
        sig_eps_sq=sig_eps_sq,
        sig_eps_c_sq=sig_eps_c_sq,
        lambda.max=lambda.max,
        lambda.max_model_size=lambda.max_model_size,
        lambda.min=lambda.min,
        lambda.min_model_size=lambda.min_model_size,
        lambda_star=lambda_star,
        lambda_star_model_size=lambda_star_model_size,
        X_ints=X_ints,
        y=y,
        X_final=X_final,
        y_final=y_final,
        N=N,
        T=T,
        R=R,
        d=d,
        p=p
        )
    )
}




###### Input: a design matrix with columns for response, time period variable
# (a single categorical or numeric/integer variable indicating the time period
# for each observation), unit variable (a single categorical variable 
# indicating which unit each observation is), treatment (a binary variable;
# 1 if unit is treated at that time and 0 if not; treatment must be an
# absorbing state in this model); covariates (which are fixed over time)

## Output: a named/identified list of cohorts (units that were treated at the
# same time); a vector of all identified units; a vector of all identified time
# periods; also the design dataframe will be returned with any units that were
# treated in the first time period removed (must only have units treated after
# first time period) and also treatment variable removed

idCohorts <- function(df, time_var, unit_var, treat_var, covs){
    stopifnot(time_var %in% colnames(df))
    stopifnot(unit_var %in% colnames(df))
    stopifnot(treat_var %in% colnames(df))
    stopifnot(all(covs %in% colnames(df)))

    # Form design matrix
    units <- unique(df[, unit_var])
    N <- length(units)
    times <- sort(unique(df[, time_var]))
    T <- length(times)

    # Variable to identify cohorts
    cohorts <- list()
    for(t in 1:T){
        cohorts[[t]] <- character()
    }
    names(cohorts) <- times

    for(s in units){
        df_s <- df[df[, unit_var] == s, ]
        # Assume this is a balanced panel
        if(nrow(df_s) != T){
            stop(paste("Panel does not appear to be balanced (unit", s,
                "does not have exactly T observations for T =", T))
        }

        if(any(df_s[, treat_var] == 1)){
            # Identify first year of treatment (and cohort)
            treat_year_s_ind <- min(which(df_s[, treat_var] == 1))
            # Make sure treatment is absorbing
            if(any(df_s[treat_year_s_ind:T, treat_var] != 1)){
                stop(paste(
                    "Treatment does not appear to be an absorbing state for unit",
                    s))
            }

            cohorts[[treat_year_s_ind]] <- c(cohorts[[treat_year_s_ind]], s)
        }
    }

    stopifnot(length(unlist(cohorts)) <= N)
    cohorts <- cohorts[lengths(cohorts) > 0]
    
    # Need at least one untreated period, so have to omit units that were
    # treated in the very first time period
    first_year_cohort <- cohorts[[as.character(times[1])]]
    df <- df[!(df[, unit_var] %in% first_year_cohort), ]

    units <- unique(df[, unit_var])

    if(length(units) == 0) {
        stop("All units were treated in the first time period; estimating treatment effects is not possible")
    }

    if(length(units) < N){
        warning(paste(N - length(units),
            "units were removed because they were treated in the first time period"))
    }
    N <- length(units)


    # Treatment no longer needed
    df <- df[, colnames(df) != treat_var]

    # Make sure there is an empty cohort for the first time
    cohorts[[as.character(times[1])]] <- character()

    # Order cohorts in order of times
    cohorts <- cohorts[order(as.numeric(names(cohorts)))]

    # This should have been the first cohort
    stopifnot(length(cohorts[[1]]) == 0)

    cohorts <- cohorts[-1]

    stopifnot(all(lengths(cohorts) >= 1))
    stopifnot(length(cohorts) <= T)
    stopifnot(length(unlist(cohorts)) <= N)

    return(list(df=df, cohorts=cohorts, units=units, times=times))
}


###### Replace time-varying variables with pre-treatment values for all times;


processCovs <- function(df, units, unit_var, times, time_var, covs, resp_var,
    T, verbose=FALSE){

    d <- length(covs)

    for(s in units){
      
        df_s <- df[df[, unit_var] == s, ]
        # Assume a balanced panel
        stopifnot(nrow(df_s) == T)

        # Identify first time covariates
        covs_s <- df_s[df_s[, time_var] == times[1], covs]
        stopifnot(length(covs_s) == length(covs))
        # Remove covariates with mising values
        covs <- covs[!is.na(covs_s)]
        if(length(covs) == 0){
            stop("All covariates were removed because for all of them, at least one unit had a missing value during the first time period")
        }
    }

    if(length(covs) < d){
        warning(paste(d - length(covs),
            "covariates were removed because they contained missing values."))
    }

    d <- length(covs)

    ###### Check if any covariates only have one unique value, and if so drop
    covs_to_remove <- character()
    for(cov in covs){
        if(length(unique(df[, cov])) == 1){
            if(verbose){
                message("Removing covariate because all units have the same value:")
                message(cov)
            }
            covs_to_remove <- c(covs_to_remove, cov)
        }
    }

    covs <- covs[!(covs %in% covs_to_remove)]
    if(length(covs) == 0){
        stop("All covariates were removed after screening out missing values in first time period and checking whether all units has the same value of a covariate")
    }
    if(length(covs) < d){
        warning(paste(d - length(covs),
            "covariates were removed because all units had the same value."))
    }

    d <- length(covs)

    # No longer need discarded covariates
    df <- df[, c(resp_var, time_var, unit_var, covs)]

    # Finally, replace all remaining covariates for each unit with value from
    # first time period
    if(verbose){
        message("For any time-varying covariates, replacing all values with value in first (pre-treatment) period")
    }
    for(s in units){
        df_s <- df[df[, unit_var] == s, ]
        stopifnot(nrow(df_s) == T)

        # Identify first year covariates
        covs_s <- df_s[df_s[, time_var] == times[1], covs]
        # Replace all covariates with first year covariates
        for(t in 1:T){
            ind_s_t <- (df[, unit_var] == s) & (df[, time_var] == times[t])
            stopifnot(sum(ind_s_t) == 1)
            df[ind_s_t, covs] <- covs_s
        }
    }

    # Sort rows: first T rows should be the observations for the first unit, and
    # so on
    df <- df[order(df[, unit_var], df[, time_var], decreasing=FALSE), ]

    return(list(df=df, covs=covs))
}

######################### Add cohort dummies and treatment variables

genTreatVarsRealData <- function(cohort_name, c_t_names, N, T, n_treated_times,
    unit_vars, time_vars, cohort, treated_times){

    # Create matrix of variables for cohort indicator and cohort/treatment
    # variables to append to df
    treat_vars <- matrix(rep(0L, N*T*(1 + n_treated_times)), N*T,
        1 + n_treated_times)

    colnames(treat_vars) <- c(cohort_name, c_t_names)

    c_i_inds <- unit_vars %in% cohort
    stopifnot(length(c_i_inds) == N*T)

    # Add cohort dummy
    treat_vars[c_i_inds, 1] <- 1

    for(t in 1:n_treated_times){
        c_i_t_inds <- c_i_inds & (time_vars == treated_times[t])

        stopifnot(length(c_i_t_inds) == N*T)
        stopifnot(sum(c_i_t_inds) == length(cohort))

        treat_vars[c_i_t_inds, t + 1] <- 1L
    }

    return(treat_vars)
}

# Names of cohorts must be the same as time of treatment

addDummies <- function(df, cohorts, times, N, T, unit_var, time_var,
    resp_var, n_cohorts){
    # Add cohort dummies, treatment dummies (for each cohort and each time after
    # treatment), and time dummies to df

    # Total number of treated times for all cohorts (later, this will be the
    # total number of treatment effects to estimate)
    num_treats <- 0

    # A list of names of indicator variables for cohorts
    cohort_vars <- character()

    # Matrix of cohort variable indicators
    cohort_var_mat <- matrix(as.integer(NA), N*T, n_cohorts)

    # Matrix of treatment indicators (just creating one column now to 
    # initalize; will delete first column later)
    treat_var_mat <- matrix(as.integer(NA), N*T, 1)

    # A list of names of variables for cohort/time treatments
    cohort_treat_names <- list()

    # Indices of first treatment time for each cohort
    first_inds <- 1

    for(i in 1:n_cohorts){

        # Time of first treatment for this cohort
        y1_treat_i <- as.integer(names(cohorts)[i])
        stopifnot(y1_treat_i <= max(times))

        treated_times_i <- y1_treat_i:max(times)

        # How many treated times are there?
        n_treated_times <- length(treated_times_i)
        stopifnot(n_treated_times <= T)
        stopifnot(n_treated_times == max(times) - y1_treat_i + 1)

        num_treats <- num_treats + n_treated_times
        first_inds <- c(first_inds, num_treats + 1)

        # Cohort/treatment time variable names
        c_i_names <- paste("c", i, treated_times_i, sep="_")
        stopifnot(length(c_i_names) == n_treated_times)

        cohort_treat_names[[i]] <- c_i_names

        cohort_vars <- c(cohort_vars, paste("c_", i, sep=""))
        stopifnot(length(cohort_vars) == i)

        treat_vars_i <- genTreatVarsRealData(cohort_name=cohort_vars[i],
            c_t_names=c_i_names, N=N, T=T, n_treated_times=n_treated_times,
            unit_vars=df[, unit_var], time_vars=df[, time_var],
            cohort=cohorts[[i]], treated_times=treated_times_i)

        # First column is cohort dummy and its name is cohort_vars[i].
        # Remaining columns are treatment dummies with names c_i_names.

        stopifnot(all(is.na(cohort_var_mat[, i])))
        stopifnot(ncol(treat_vars_i) == n_treated_times + 1)

        cohort_var_mat[, i] <- treat_vars_i[, 1]
        treat_var_mat <- cbind(treat_var_mat, treat_vars_i[,
            2:(n_treated_times + 1)])

        df <- data.frame(df, treat_vars_i)
    }

    stopifnot(length(first_inds) == n_cohorts + 1)
    first_inds <- first_inds[1:n_cohorts]
    stopifnot(all(first_inds %in% 1:num_treats))

    stopifnot(length(cohort_treat_names) == n_cohorts)
    names(cohort_treat_names) <- names(cohorts)

    # treat_var_mat should have one extra column at the beginning that we used
    # to initalize it
    stopifnot(is.numeric(ncol(treat_var_mat)) | is.integer(ncol(treat_var_mat)))
    stopifnot(ncol(treat_var_mat) == num_treats + 1)
    treat_var_mat <- treat_var_mat[, 2:(num_treats + 1), drop=FALSE]

    # Add time dummies for all but first time
    time_var_mat <- matrix(0L, N*T, T - 1)
    time_var_names <- paste("t", times[2:T], sep="_")
    colnames(time_var_mat) <- time_var_names
    for(t in 2:T){
        c_t_inds <- df[, time_var] %in% times[t]
        stopifnot(length(c_t_inds) == N*T)
        stopifnot(sum(c_t_inds) == N)

        # Add cohort dummy
        time_var_mat[c_t_inds, t - 1] <- 1
    }

    # Assume balanced panel
    stopifnot(all(colSums(time_var_mat) == N))

    stopifnot(length(unlist(cohort_treat_names)) == num_treats)
    # stopifnot(all(unlist(cohort_treat_names) %in% colnames(X_df)))
    stopifnot(ncol(cohort_var_mat) == n_cohorts)

    stopifnot(length(time_var_names) == T - 1)
    stopifnot(ncol(time_var_mat) == T - 1)
    # stopifnot(all(time_var_names %in% colnames(X_df)))

    stopifnot(ncol(treat_var_mat) == num_treats)
    stopifnot(ncol(treat_var_mat) >= 1)

    stopifnot(length(cohort_vars) == n_cohorts)
    # stopifnot(all(cohort_vars %in% colnames(X_df)))

    # Center response
    y <- df[, resp_var] - mean(df[, resp_var])

    return(list(time_var_mat=time_var_mat, cohort_var_mat=cohort_var_mat,
        treat_var_mat=treat_var_mat, y=y, cohort_treat_names=cohort_treat_names,
        time_var_names=time_var_names, cohort_vars=cohort_vars,
        first_inds=first_inds))
}

generateFEInts <- function(X_long, cohort_fe, time_fe, N, T, R, d){
    # Interact with cohort effects
    X_long_cohort <- matrix(as.numeric(NA), nrow=N*T, ncol=R*d)

    stopifnot(ncol(cohort_fe) == R)
    stopifnot(nrow(cohort_fe) == N*T)
    stopifnot(ncol(time_fe) == T - 1)
    stopifnot(ncol(X_long) == d)
    stopifnot(is.matrix(X_long))
    stopifnot(!is.data.frame(X_long))

    for(r in 1:R){
        # Notice that these are arranged one cohort at a time, interacted with
        # all covariates
        first_col_r <- (r - 1)*d + 1
        last_col_r <- r*d

        stopifnot(last_col_r - first_col_r + 1 == d)

        stopifnot(all(is.na(X_long_cohort[, first_col_r:last_col_r])))

        stopifnot(ncol(cohort_fe[, r]*X_long) == length(first_col_r:last_col_r))

        X_long_cohort[, first_col_r:last_col_r] <- cohort_fe[, r]*X_long

    }

    stopifnot(all(!is.na(X_long_cohort)))
    stopifnot(nrow(X_long_cohort) == N*T)
    stopifnot(ncol(X_long_cohort) == R*d)

    # Interact with treatment effects
    X_long_time <- matrix(as.numeric(NA), nrow=N*T, ncol=(T - 1)*d)

    for(t in 1:(T - 1)){
        # Notice that these are arranged one time at a time, interacted with all
        # covariates
        first_col_t <- (t - 1)*d + 1
        last_col_t <- t*d

        stopifnot(last_col_t - first_col_t + 1 == d)
        stopifnot(all(is.na(X_long_time[, first_col_t:last_col_t])))

        X_long_time[, first_col_t:last_col_t] <- time_fe[, t]*X_long
    }

    stopifnot(all(!is.na(X_long_time)))
    stopifnot(nrow(X_long_time) == N*T)
    stopifnot(ncol(X_long_time) == (T - 1)*d)

    return(list(X_long_cohort=X_long_cohort, X_long_time=X_long_time))
}


genTreatInts <- function(treat_mat_long, X_long, n_treats, cohort_fe, N, T, R,
    d, N_UNTREATED){

    stopifnot(ncol(cohort_fe) == R)
    stopifnot(ncol(treat_mat_long) == n_treats)

    stopifnot(is.numeric(d) | is.integer(d))
    stopifnot(d >= 1)

    stopifnot(is.numeric(n_treats) | is.integer(n_treats))
    stopifnot(n_treats >= 1)

    stopifnot(is.numeric(N) | is.integer(N))
    stopifnot(N >= 2)

    stopifnot(is.numeric(T) | is.integer(T))
    stopifnot(T >= 2)
    # Interact with covariates
    X_long_treat <- matrix(as.numeric(NA), nrow=N*T, ncol=d*n_treats)

    # First need to center covariates with respect to cohort means
    X_long_centered <- matrix(as.numeric(NA), nrow=N*T, ncol=d)

    # Never treated group
    unt_row_inds <- which(rowSums(cohort_fe) == 0)

    # Assume balanced panel
    stopifnot(length(unt_row_inds) == N_UNTREATED*T)

    stopifnot(all(is.na(X_long_centered[unt_row_inds, ])))

    X_long_centered[unt_row_inds, ] <- scale(X_long[unt_row_inds, ],
        center=TRUE, scale=FALSE)

    for(r in 1:R){
        # Notice that these are arranged column-wise one cohort at a time,
        # interacted with all treatment effects

        cohort_r_inds <- which(cohort_fe[, r] == 1)

        new_mat <- scale(X_long[cohort_r_inds, ], center=TRUE, scale=FALSE)

        stopifnot(all(!is.na(new_mat)))
        stopifnot(all(is.na(X_long_centered[cohort_r_inds, ])))

        X_long_centered[cohort_r_inds, ] <- new_mat
            
    }
    # Final index should correspond to last row
    stopifnot(all(!is.na(X_long_centered)))

    for(i in 1:n_treats){
        first_col_i <- (i - 1)*d + 1
        last_col_i <- i*d

        stopifnot(all(is.na(X_long_treat[, first_col_i:last_col_i])))

        X_long_treat[, first_col_i:last_col_i] <- treat_mat_long[, i]*
            X_long_centered
    }

    stopifnot(all(!is.na(X_long_treat)))
    stopifnot(nrow(X_long_treat) == N*T)
    stopifnot(ncol(X_long_treat) == n_treats*d)

    return(X_long_treat)
}

genXintsData <- function(cohort_fe, time_fe, X_long, treat_mat_long, N, R, T,
    d, N_UNTREATED, p){

    X_int <- cbind(cohort_fe, time_fe, X_long)
    stopifnot(ncol(X_int) == R + T - 1 + d)

    # Generate interactions of X with time and cohort fixed effects
    # Note: columns of X_long_cohort are arranged in blocks of size d for one
    # cohort at a time (that is, the first d columns are the interactions of
    # all d features with the indicator variables for the first block, and
    # so on). Similarly, the columns of X_long_time are arranged in T - 1
    # blocks of size d.
    res <- generateFEInts(X_long, cohort_fe, time_fe, N, T, R, d)

    X_long_cohort <- res$X_long_cohort
    X_long_time <- res$X_long_time

    X_int <- cbind(X_int, X_long_cohort, X_long_time, treat_mat_long)

    rm(res)

    stopifnot(is.integer(ncol(treat_mat_long)) | is.numeric(ncol(treat_mat_long)))
    stopifnot(ncol(treat_mat_long) >= 1)

    # Generate interactions of treatment effects with X
    X_long_treat <- genTreatInts(
        treat_mat_long=treat_mat_long, 
        X_long=X_long, 
        n_treats=ncol(treat_mat_long),
        cohort_fe,
        N=N,
        T=T,
        R=R,
        d=d,
        N_UNTREATED=N_UNTREATED
        )

    # The first R columns are the cohort fixed effects in order
    # Next T are time fixed effects in order
    # Next d are X
    # Next d*R are cohort effects interacted with X. (The first d of these
    # are first cohort effects interacted with X, and so on until the Rth
    # cohort.)
    # Next d*T are time effects interacted with X (similarly to the above,
    # the first d are first time effects interacted with X, and so on until
    # the Tth time).
    # Next num_treats = R*T - R*(R + 1)/2 columns are base treatment effects
    # (for each cohort and time)
    # Finally, the next num_treats*d are interactions of all of these
    # treatment effects over time--first d are the first column of
    # treat_mat_long interacted with all of the columns of X, and so on.
    X_int <- cbind(X_int, X_long_treat)

    stopifnot(ncol(X_int) == p)
    stopifnot(nrow(X_int) == N * T)


    return(X_int)
}


getFirstIndsDataApp <- function(n_treats, R){
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

transformXintDataApp <- function(X_int, N, T, R, d, num_treats, first_inds=NA){
    

    p = ncol(X_int)

    X_mod <- matrix(as.numeric(NA), nrow=N*T, ncol=p)

    stopifnot(nrow(X_int) == N*T)

    # Handle cohort fixed effects:
    X_mod[, 1:R] <- X_int[, 1:R] %*% genBackwardsInvFusionTransformMatDataApp(R)

    # Time fixed effects
    stopifnot(all(is.na(X_mod[, (R + 1):(R + T - 1)])))
    X_mod[, (R + 1):(R + T - 1)] <- X_int[, (R + 1):(R + T - 1)] %*%
        genBackwardsInvFusionTransformMatDataApp(T - 1)

    # X goes in unpenalized
    stopifnot(all(is.na(X_mod[, (R + T - 1 + 1):(R + T - 1 + d)])))
    X_mod[, (R + T - 1 + 1):(R + T - 1 + d)] <- X_int[, (R + T - 1 + 1):
        (R + T - 1 + d)]

    stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d)])))
    stopifnot(all(is.na(X_mod[, (R + T - 1 + d + 1):p])))

    # For cohort effects interacted with X: we have d*R columns to deal with.
    # For each individual feature, this will be handled using
    # genTransformedMatFusion.
    if(any(is.na(first_inds))){
        first_inds <- getFirstIndsDataApp(num_treats, R)
    }

    for(j in 1:d){
        # Get indices corresponding to interactions between feature j and cohort
        # fixed effects--these are the first feature, the (1 + d)th feature, and
        # so on R times
        feat_1 <- R + T - 1 + d + j
        feat_R <- R + T - 1 + d + (R - 1)*d + j
        feat_inds_j <- seq(feat_1, feat_R, by=d)
        stopifnot(length(feat_inds_j) == R)

        stopifnot(all(is.na(X_mod[, feat_inds_j])))

        X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
            genBackwardsInvFusionTransformMatDataApp(R)

    }
    stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R*d)])))
    stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R*d + 1):p])))

    # Similar for time effects interacted with X

    for(j in 1:d){
        # Get indices corresponding to interactions between feature j and time
        # fixed effects--these are the first feature, the (1 + d)th feature, and
        # so on T - 1 times
        feat_1 <- R + T - 1 + d + R*d + j
        feat_T_minus_1 <- R + T - 1 + d + R*d + (T - 2)*d + j
        feat_inds_j <- seq(feat_1, feat_T_minus_1, by=d)
        stopifnot(length(feat_inds_j) == T - 1)

        stopifnot(all(is.na(X_mod[, feat_inds_j])))
        X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
            genBackwardsInvFusionTransformMatDataApp(T - 1)

    }
    stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R*d + (T - 1)*d)])))
    stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R*d + (T - 1)*d + 1):p])))

    # Now base treatment effects. For each cohort, will penalize base term, then
    # fuse remaining terms toward it. Also, for each cohort, will penalize base
    # treatment effect of this cohort to base of previous cohort. New function
    # genTransformedMatTwoWayFusion does this.

    feat_inds <- (R + T - 1 + d + R*d + (T - 1)*d + 1):
        (R + T - 1 + d + R*d + (T - 1)*d + num_treats)

    # Now ready to generate the appropriate transformed matrix
    stopifnot(all(is.na(X_mod[, feat_inds])))

    X_mod[, feat_inds] <- X_int[, feat_inds] %*%
        genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

    stopifnot(all(!is.na(X_mod[, 1:(R + T - 1 + d + R*d + (T - 1)*d + num_treats)])))
    stopifnot(all(is.na(X_mod[, (R + T - 1 + d + R*d + (T - 1)*d + num_treats + 1):p])))

    # Lastly, penalize interactions between each treatment effect and each feature.
    # Feature-wise, we can do this with genTransformedMatTwoWayFusion, in the same
    # way that we did for previous interactions with X.
    for(j in 1:d){
        # Recall that we have arranged the last d*num_Feats features in X_int
        # as follows: the first d are the first column of treat_mat_long interacted
        # with all of the columns of X, and so on. So, the columns that interact
        # the jth feature with all of the treatment effects are columns j, j + 1*d,
        # j + 2*d, ..., j + (num_treats - 1)*d.
        inds_j <- seq(j, j + (num_treats - 1)*d, by=d)
        stopifnot(length(inds_j) == num_treats)
        inds_j <- inds_j + R + T - 1 + d + R*d + (T - 1)*d + num_treats

        # Now ready to generate the appropriate transformed matrix
        stopifnot(all(is.na(X_mod[, inds_j])))

        X_mod[, inds_j] <- X_int[, inds_j] %*%
            genInvTwoWayFusionTransformMat(num_treats, first_inds, R)

        stopifnot(all(!is.na(X_mod[, inds_j])))
    }

    stopifnot(all(!is.na(X_mod)))
    stopifnot(ncol(X_mod) == p)
    stopifnot(nrow(X_mod) == N * T)

    return(X_mod)
}


untransformCoefsDataApp <- function(beta_hat_mod, T, R, p, d, num_treats,
    first_inds=NA){

    stopifnot(length(beta_hat_mod) == p)
    beta_hat <- rep(as.numeric(NA), p)

    if(any(is.na(first_inds))){
        first_inds <- getFirstIndsDataApp(num_treats, R)
    }
    
    # First handle R cohort fixed effects effects
    beta_hat[1:R] <- genBackwardsInvFusionTransformMatDataApp(R) %*% beta_hat_mod[1:R] 

    stopifnot(all(!is.na(beta_hat[1:R])))
    stopifnot(all(is.na(beta_hat[(R + 1):p])))

    # Next, T - 1 time fixed effects
    beta_hat[(R + 1):(R + T - 1)] <- genBackwardsInvFusionTransformMatDataApp(T - 1) %*%
        beta_hat_mod[(R + 1):(R + T - 1)]

    stopifnot(all(!is.na(beta_hat[1:(R + T - 1)])))
    stopifnot(all(is.na(beta_hat[(R + T):p])))

    # Next, coefficients for X stay the same

    beta_hat[(R + T):(R + T - 1 + d)] <- beta_hat_mod[(R + T):(R + T - 1 + d)]

    stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d)])))
    stopifnot(all(is.na(beta_hat[(R + T + d):p])))

    # Next, coefficients for cohort effects interacted with X. For each individual
    # feature, this will be handled using untransformVecFusion.
    for(j in 1:d){
        # Get indices corresponding to interactions between feature j and cohort
        # fixed effects--these are the first feature, the (1 + d)th feature, and
        # so on R times
        feat_1 <- R + T - 1 + d + j
        feat_R <- R + T - 1 + d + (R - 1)*d + j
        feat_inds_j <- seq(feat_1, feat_R, by=d)
        stopifnot(length(feat_inds_j) == R)

        stopifnot(all(is.na(beta_hat[feat_inds_j])))

        beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMatDataApp(R) %*% 
            beta_hat_mod[feat_inds_j]
        stopifnot(all(!is.na(beta_hat[feat_inds_j])))
    }

    stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d + R*d)])))
    stopifnot(all(is.na(beta_hat[(R + T - 1 + d + R*d + 1):p])))

    # Similar for time effects interacted with X

    for(j in 1:d){
        # Get indices corresponding to interactions between feature j and time
        # fixed effects--these are the first feature, the (1 + d)th feature, and
        # so on T - 1 times
        feat_1 <- R + T - 1 + d + R*d + j
        feat_T_minus_1 <- R + T - 1 + d + R*d + (T - 2)*d + j
        feat_inds_j <- seq(feat_1, feat_T_minus_1, by=d)
        stopifnot(length(feat_inds_j) == T - 1)

        stopifnot(all(is.na(beta_hat[feat_inds_j])))

        beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMatDataApp(T - 1) %*%
            beta_hat_mod[feat_inds_j]
        stopifnot(all(!is.na(beta_hat[feat_inds_j])))
    }
    stopifnot(all(!is.na(beta_hat[1:(R + T - 1 + d + R*d + (T - 1)*d)])))
    stopifnot(all(is.na(beta_hat[(R + T - 1 + d + R*d + (T - 1)*d + 1):p])))

    # Now base treatment effects.

    feat_inds <- (R + T - 1 + d + R*d + (T - 1)*d + 1):
        (R + T - 1 + d + R*d + (T - 1)*d + num_treats)

    stopifnot(all(is.na(beta_hat[feat_inds])))

    beta_hat[feat_inds] <- genInvTwoWayFusionTransformMat(num_treats,
        first_inds, R) %*% beta_hat_mod[feat_inds]

    stopifnot(all(!is.na(beta_hat[1:
        (R + T - 1 + d + R*d + (T - 1)*d + num_treats)])))
    stopifnot(all(is.na(beta_hat[(R + T - 1 + d + R*d + (T - 1)*d + num_treats +
        1):p])))

    # Lastly, interactions between each treatment effect and each feature.
    # Feature-wise, we can do this with untransformTwoWayFusionCoefs, in the same
    # way that we did for previous interactions with X.
    for(j in 1:d){
        # Recall that we have arranged the last d*num_Feats features in X_int
        # as follows: the first d are the first column of treat_mat_long interacted
        # with all of the columns of X, and so on. So, the columns that interact
        # the jth feature with all of the treatment effects are columns j, j + 1*d,
        # j + 2*d, ..., j + (num_treats - 1)*d.
        inds_j <- seq(j, j + (num_treats - 1)*d, by=d)
        stopifnot(length(inds_j) == num_treats)
        inds_j <- inds_j + R + T - 1 + d + R*d + (T - 1)*d + num_treats

        # Now ready to untransform the estimated coefficients
        stopifnot(all(is.na(beta_hat[inds_j])))

        beta_hat[inds_j] <- genInvTwoWayFusionTransformMat(num_treats,
            first_inds, R) %*% beta_hat_mod[inds_j]

        stopifnot(all(!is.na(beta_hat[inds_j])))
    }

    stopifnot(all(!is.na(beta_hat)))

    return(beta_hat)
}

sse_bridge <- function(eta_hat, beta_hat, y, X_mod, N, T){
    stopifnot(length(eta_hat) == 1)
    stopifnot(length(beta_hat) == ncol(X_mod))

    y_hat <- X_mod %*% beta_hat + eta_hat
    stopifnot(length(y_hat) == N*T)
    stopifnot(all(!is.na(y_hat)))

    ret <- sum((y - y_hat)^2)/(N*T)

    stopifnot(!is.na(ret))
    stopifnot(ret >= 0)

    return(ret)
}



genBackwardsFusionTransformMat <- function(n_vars){
    # Generates D matrix in relation theta = D beta, where D beta is what
    # we want to penalize (for a single set of coefficients where we want to
    # penalize last coefficient directly and penalize remaining coefficints
    # towards the next coefficient)
    D <- matrix(0, n_vars, n_vars)


    for(i in 1:n_vars){
        for(j in 1:n_vars){
            if(i == j){
                D[i, j] <- 1
            }
            if(j == i + 1){
                D[i, j] <- -1
            }
        }
    }

    return(D)
}

genBackwardsInvFusionTransformMatDataApp <- function(n_vars){
    # Generates inverse of D matrix in relation theta = D beta, where D beta is
    # what we want to penalize (for a single set of coefficients where we want
    # to penalize last coefficient directly and penalize remaining coefficints
    # towards the next coefficient)
    D_inv <- matrix(0, n_vars, n_vars)

    diag(D_inv) <- 1

    D_inv[upper.tri(D_inv)] <- 1

    return(D_inv)
}


genInvFusionTransformMat <- function(n_vars){
    # Generates inverse of D matrix in relation theta = D beta, where D beta is
    # what we want to penalize (for a single set of coefficients where we want
    # to penalize first coefficient directly and penalize remaining coefficints
    # towards the previous coefficient)
    D_inv <- matrix(0, n_vars, n_vars)

    diag(D_inv) <- 1

    D_inv[lower.tri(D_inv)] <- 1

    return(D_inv)
}



estOmegaSqrtInv <- function(y, X_ints, N, T, p){

    if(N*(T - 1) - p <= 0){
        stop("Not enough units available to estimate the noise variance.")
    }
    stopifnot(N > 1)

    # Estimate standard deviations
    y_fe <- y
    X_ints_fe <- X_ints
    for(i in 1:N){
        inds_i <- ((i - 1)*T + 1):(i*T)
        stopifnot(length(inds_i) == T)
        y_fe[inds_i] <- y[inds_i] - mean(y[inds_i])
        X_ints[inds_i, ] <- X_ints[inds_i, ] - colMeans(X_ints[inds_i, ])
    }

    lin_mod_fe <- glmnet::cv.glmnet(x=X_ints_fe, y=y_fe, alpha=0)

    # Get residuals
    y_hat_fe <- stats::predict(lin_mod_fe, s="lambda.min", newx=X_ints_fe)

    # Get coefficients
    beta_hat_fe <- stats::coef(lin_mod_fe, s="lambda.min")[2:(p + 1)]

    resids <- y_fe - y_hat_fe

    tau <- rep(1, T)
    M <- diag(rep(1, T)) - outer(tau, tau)/T

    sigma_hat_sq <- 0
    alpha_hat <- rep(as.numeric(NA), N)

    for(i in 1:N){
        inds_i <- ((i - 1)*T + 1):(i*T)
        sigma_hat_sq <- sigma_hat_sq + as.numeric(resids[inds_i] %*% M %*%
            resids[inds_i])
        alpha_hat[i] <- mean(y_hat_fe[inds_i]) -
            colMeans(X_ints_fe[inds_i, ] %*% beta_hat_fe)
    }

    stopifnot(all(!is.na(alpha_hat)))

    sigma_hat_sq <- sigma_hat_sq/(N*(T - 1) - p)

    sigma_c_sq_hat <- sum((alpha_hat - mean(alpha_hat))^2)/(N - 1)

    return(list(sig_eps_sq=sigma_hat_sq, sig_eps_c_sq=sigma_c_sq_hat))

}


getSecondVarTermDataApp <- function(cohort_probs, psi_mat,
    sel_treat_inds_shifted, tes, d_inv_treat_sel, cohort_probs_overall,
    first_inds, theta_hat_treat_sel, num_treats, N, T, R){

    stopifnot(ncol(d_inv_treat_sel) == length(sel_treat_inds_shifted))
    stopifnot(length(theta_hat_treat_sel) == length(sel_treat_inds_shifted))
    Sigma_pi_hat <- -outer(cohort_probs_overall[1:(R)],
        cohort_probs_overall[1:(R)])
    diag(Sigma_pi_hat) <- cohort_probs_overall[1:(R)]*
        (1 - cohort_probs_overall[1:(R)])

    stopifnot(nrow(Sigma_pi_hat) == R)
    stopifnot(ncol(Sigma_pi_hat) == R)

    # Jacobian
    jacobian_mat <- matrix(as.numeric(NA), nrow=R,
        ncol=length(sel_treat_inds_shifted))

    sel_inds <- list()

    for(r in 1:R){
        first_ind_r <- first_inds[r]

        if(r < R){
            last_ind_r <- first_inds[r + 1] - 1
        }
        else{
            last_ind_r <- num_treats
        }
        stopifnot(last_ind_r >= first_ind_r)
        sel_inds[[r]] <- first_ind_r:last_ind_r
        if(r > 1){
            stopifnot(min(sel_inds[[r]]) > max(sel_inds[[r - 1]]))
            stopifnot(length(sel_inds[[r]]) < length(sel_inds[[r - 1]]))
        }
    }
    stopifnot(all.equal(unlist(sel_inds), 1:num_treats))

    for(r in 1:R){
        cons_r <- (sum(cohort_probs_overall) -
            cohort_probs_overall[r]) / sum(cohort_probs_overall)^2 

        if(length(sel_treat_inds_shifted) > 1){
            jacobian_mat[r, ] <- cons_r *
                colMeans(d_inv_treat_sel[sel_inds[[r]], ])
        } else{
            jacobian_mat[r, ] <- cons_r *
                mean(d_inv_treat_sel[sel_inds[[r]], ])
        }
        for(r_double_prime in setdiff(1:R, r)){
            cons_r_double_prime <- (sum(cohort_probs_overall) -
                cohort_probs_overall[r_double_prime]) / sum(cohort_probs_overall)^2 

            if(length(sel_treat_inds_shifted) > 1){
                jacobian_mat[r, ] <- jacobian_mat[r, ] - cons_r_double_prime *
                    colMeans(d_inv_treat_sel[sel_inds[[r_double_prime]], ])
            } else{
                jacobian_mat[r, ] <- jacobian_mat[r, ] - cons_r_double_prime *
                    mean(d_inv_treat_sel[sel_inds[[r_double_prime]], ])
            }
        }
    }

    stopifnot(all(!is.na(jacobian_mat)))

    att_var_2 <- T * as.numeric(t(theta_hat_treat_sel) %*%
        t(jacobian_mat) %*% Sigma_pi_hat %*% jacobian_mat %*%
        theta_hat_treat_sel) / (N * T)

    return(att_var_2)
}

getCohortATTsFinal <- function(
    X_final,
    sel_feat_inds,
    treat_inds,
    num_treats,
    first_inds,
    sel_treat_inds_shifted,
    c_names,
    tes,
    sig_eps_sq,
    R,
    N,
    T,
    fused,
    calc_ses,
    alpha=0.05
    ){

    stopifnot(max(sel_treat_inds_shifted) <= num_treats)
    stopifnot(min(sel_treat_inds_shifted) >= 1)
    stopifnot(length(tes) == num_treats)
    stopifnot(all(!is.na(tes)))

    # Start by getting Gram matrix needed for standard errors
    if(calc_ses){
        res <- getGramInv(N, T, X_final, sel_feat_inds, treat_inds, num_treats,
            sel_treat_inds_shifted, calc_ses)

        gram_inv <- res$gram_inv
        calc_ses <- res$calc_ses
    } else{
        gram_inv <- NA
    }
    

    if(fused){
        # Get the parts of D_inv that have to do with treatment effects
        d_inv_treat <- genInvTwoWayFusionTransformMat(num_treats, first_inds, R)
    }
    
    # First, each cohort
    cohort_tes <- rep(as.numeric(NA), R)
    cohort_te_ses <- rep(as.numeric(NA), R)

    psi_mat <- matrix(0, length(sel_treat_inds_shifted), R)

    d_inv_treat_sel <- matrix(0, nrow=0, ncol=length(sel_treat_inds_shifted))

    for(r in 1:R){
        # Get indices corresponding to rth treatment
        first_ind_r <- first_inds[r]
        if(r < R){
            last_ind_r <- first_inds[r + 1] - 1
        } else{
            last_ind_r <- num_treats
        }

        stopifnot(last_ind_r >= first_ind_r)
        stopifnot(all(first_ind_r:last_ind_r %in% 1:num_treats))
        
        cohort_tes[r] <- mean(tes[first_ind_r:last_ind_r])

        if(calc_ses){
            # Calculate standard errors

            if(fused){
                res_r <- getPsiRFused(first_ind_r, last_ind_r,
                    sel_treat_inds_shifted, d_inv_treat)

                psi_r <- res_r$psi_r 

                stopifnot(nrow(res_r$d_inv_treat_sel) == last_ind_r -
                    first_ind_r + 1)
                stopifnot(ncol(res_r$d_inv_treat_sel) ==
                    length(sel_treat_inds_shifted))

                stopifnot(is.matrix(res_r$d_inv_treat_sel))

                d_inv_treat_sel <- rbind(d_inv_treat_sel, res_r$d_inv_treat_sel)

                if(nrow(d_inv_treat_sel) != last_ind_r){
                    err_mes <- paste("nrow(d_inv_treat_sel) == last_ind_r is not TRUE. ",
                        "nrow(d_inv_treat_sel): ", nrow(d_inv_treat_sel), ". num_treats: ",
                        num_treats, ". R: ", R, ". first_inds: ",
                        paste(first_inds, collapse=", "), ". r: ", r,
                        ". first_ind_r: ", first_ind_r, ". last_ind_r: ",
                        last_ind_r, ". nrow(res_r$d_inv_treat_sel):",
                        nrow(res_r$d_inv_treat_sel))
                    stop(err_mes)
                }
                
                rm(res_r)
            } else{
                psi_r <- getPsiRUnfused(first_ind_r, last_ind_r,
                    sel_treat_inds_shifted, gram_inv)

            }

            stopifnot(length(psi_r) == length(sel_treat_inds_shifted))

            psi_mat[, r] <- psi_r
            # Get standard errors
            
            cohort_te_ses[r] <- sqrt(sig_eps_sq * as.numeric(t(psi_r) %*%
                gram_inv %*% psi_r) / (N * T))
        } 
    }

    if(fused & calc_ses){
        if(nrow(d_inv_treat_sel) != num_treats){
            err_mes <- paste("nrow(d_inv_treat_sel) == num_treats is not TRUE. ",
                "nrow(d_inv_treat_sel): ", nrow(d_inv_treat_sel), ". num_treats: ",
                num_treats, ". R: ", R, ". first_inds: ",
                paste(first_inds, collapse=", "), ".")
            stop(err_mes)
        }
    }

    stopifnot(length(c_names) == R)
    stopifnot(length(cohort_tes) == R)

    if(calc_ses & all(!is.na(gram_inv))){
        stopifnot(length(cohort_te_ses) == R)

        cohort_te_df <- data.frame(c_names, cohort_tes, cohort_te_ses,
            cohort_tes - stats::qnorm(1 - alpha/2)*cohort_te_ses,
            cohort_tes + stats::qnorm(1 - alpha/2)*cohort_te_ses)

        names(cohort_te_ses) <- c_names
        names(cohort_tes) <- c_names

    } else{
        cohort_te_df <- data.frame(c_names, cohort_tes, rep(NA, R),
            rep(NA, R), rep(NA, R))
    }
    
    colnames(cohort_te_df) <- c("Cohort", "Estimated TE", "SE", "ConfIntLow",
        "ConfIntHigh")

    if(fused){
        stopifnot(is.matrix(d_inv_treat_sel))
        ret <- list(cohort_te_df=cohort_te_df, cohort_tes=cohort_tes,
            cohort_te_ses=cohort_te_ses, psi_mat=psi_mat, gram_inv=gram_inv,
            d_inv_treat_sel=d_inv_treat_sel, calc_ses=calc_ses)
    } else{
        ret <- list(cohort_te_df=cohort_te_df, cohort_tes=cohort_tes,
            cohort_te_ses=cohort_te_ses, psi_mat=psi_mat, gram_inv=gram_inv,
            calc_ses=calc_ses)
    }
    return(ret)
}

getPsiRUnfused <- function(first_ind_r, last_ind_r, sel_treat_inds_shifted,
    gram_inv){

    which_inds_ir <- sel_treat_inds_shifted %in% (first_ind_r:last_ind_r)

    psi_r <- rep(0, length(sel_treat_inds_shifted))

    if(sum(which_inds_ir) > 0){

        inds_r <- which(which_inds_ir)

        stopifnot(is.integer(inds_r) | is.numeric(inds_r))
        stopifnot(identical(inds_r, as.integer(round(inds_r))))
        stopifnot(length(inds_r) >= 1)
        stopifnot(length(inds_r) == length(unique(inds_r)))
        stopifnot(length(inds_r) <= length(sel_treat_inds_shifted))
        stopifnot(all(inds_r %in% 1:length(sel_treat_inds_shifted)))

        stopifnot(max(inds_r) <= nrow(gram_inv))
        stopifnot(max(inds_r) <= ncol(gram_inv))
        stopifnot(min(inds_r) >= 0)

        psi_r[inds_r] <- 1

        stopifnot(sum(psi_r) > 0)

        psi_r <- psi_r/sum(psi_r)
    }

    return(psi_r)
}

genInvTwoWayFusionTransformMat <- function(n_vars, first_inds, R){
    stopifnot(length(n_vars) == 1)
    stopifnot(length(first_inds) == R)
    D <- matrix(0, n_vars, n_vars)

    diag(D) <- 1

    if(R < 2){
        stop("Only one treated cohort detected in data. Currently fetwfe only supports data sets with at least two treated cohorts.")
    }

    for(j in 1:(R - 1)){
        index_j <- first_inds[j]
        next_index <- first_inds[j + 1]

        D[index_j:n_vars, index_j] <- 1
        D[next_index:n_vars, next_index] <- 1

        if(index_j + 1 <= next_index - 1){
            for(k in (index_j + 1):(next_index - 1)){
                D[k, (index_j + 1):k] <- 1
            }
        }
    }

    
    return(D)
}

getBetaBIC <- function(fit, N, T, p, X_mod, y){
    stopifnot(length(y) == N*T)
    n_lambda <- ncol(fit$beta)
    BICs <- rep(as.numeric(NA), n_lambda)
    model_sizes <- rep(as.integer(NA), n_lambda)

    stopifnot(nrow(fit$beta) == p + 1)

    for(k in 1:n_lambda){
        eta_hat_k <- fit$beta[1, k]
        beta_hat_k <- fit$beta[2:(p + 1), k]
        # Residual sum of squares
        mse_hat <- sse_bridge(eta_hat_k, beta_hat_k, y=y, X_mod=X_mod, N=N, T=T)
        # Number of fitted coefficients
        s <- sum(fit$beta[, k]!= 0)
        model_sizes[k] <- s

        stopifnot(is.na(BICs[k]))
        BICs[k] <- N*T*log(mse_hat) + s*log(N*T)
    }

    lambda_star_ind <- which(BICs == min(BICs))
    if(length(lambda_star_ind) == 1){
        lambda_star_final_ind <- lambda_star_ind
        theta_hat <- fit$beta[, lambda_star_final_ind]
    } else{
        # Choose smallest model size among models with equal BIC
        model_sizes_star <- model_sizes[lambda_star_ind]
        min_model_size_ind <- which(model_sizes_star == min(model_sizes_star))
        lambda_star_final_ind <- lambda_star_ind[min_model_size_ind][1]
        stopifnot(length(lambda_star_final_ind) == 1)
        theta_hat <- fit$beta[, lambda_star_final_ind]
    }
    stopifnot(length(lambda_star_final_ind) == 1)
    stopifnot(length(theta_hat) == p + 1)
    stopifnot(all(!is.na(theta_hat)))

    return(list(
        theta_hat=theta_hat,
        lambda_star_ind=lambda_star_final_ind,
        lambda_star_model_size=model_sizes[lambda_star_final_ind]
        )
    )
}

getGramInv <- function(N, T, X_final, sel_feat_inds, treat_inds, num_treats,
    sel_treat_inds_shifted, calc_ses){

    stopifnot(nrow(X_final) == N * T)
    X_sel <- X_final[, sel_feat_inds]
    X_sel_centered <- scale(X_sel, center=TRUE, scale=FALSE)

    gram <- 1/(N*T)*(t(X_sel_centered) %*% X_sel_centered)

    stopifnot(nrow(gram) == length(sel_feat_inds))
    stopifnot(ncol(gram) == length(sel_feat_inds))

    min_gram_eigen <- min(eigen(gram, symmetric = TRUE,
                only.values = TRUE)$values)

    if(min_gram_eigen < 10^(-12)){
        warning("Gram matrix corresponding to selected features is not invertible. Assumptions needed for inference are not satisfied. Standard errors will not be calculated.")
        return(list(gram_inv=NA, calc_ses=FALSE))
    }

    gram_inv <- solve(gram)

    # Get only the parts of gram_inv that have to do with treatment effects
    sel_treat_inds <- sel_feat_inds %in% treat_inds

    stopifnot(is.logical(sel_treat_inds))
    stopifnot(sum(sel_treat_inds) <= length(sel_feat_inds))
    stopifnot(length(sel_feat_inds) == length(sel_feat_inds))

    gram_inv <- gram_inv[sel_treat_inds, sel_treat_inds]

    stopifnot(nrow(gram_inv) <= num_treats)
    stopifnot(nrow(gram_inv) <= length(sel_feat_inds))
    stopifnot(nrow(gram_inv) == ncol(gram_inv))
    stopifnot(nrow(gram_inv) == length(sel_treat_inds_shifted))

    return(list(gram_inv=gram_inv, calc_ses=calc_ses))
}

getPsiRFused <- function(first_ind_r, last_ind_r, sel_treat_inds_shifted,
    d_inv_treat){

    stopifnot(length(sel_treat_inds_shifted) >= 0)
    stopifnot(last_ind_r >= first_ind_r)
    # Get psi vector: the part of D inverse that we need to look at is the
    # block corresponding to the treatment effect estimates, which is the
    # num_treats x num_treats matrix yielded by
    # genInvTwoWayFusionTransformMat(num_treats, first_inds).

    # Correct rows of matrix

    if(last_ind_r > first_ind_r){
        if(length(sel_treat_inds_shifted) > 1){
            psi_r <- colMeans(d_inv_treat[first_ind_r:last_ind_r,
                sel_treat_inds_shifted])
        } else{
            psi_r <- mean(d_inv_treat[first_ind_r:last_ind_r,
                sel_treat_inds_shifted])
        }
        if(length(sel_treat_inds_shifted) == 1){
            # Need to coerce this object to be a matrix with one column so it
            # works smoothly with rbind() later
            d_inv_treat_sel <- matrix(d_inv_treat[first_ind_r:last_ind_r,
                sel_treat_inds_shifted], ncol=1)

        } else{
            d_inv_treat_sel <- d_inv_treat[first_ind_r:last_ind_r,
                sel_treat_inds_shifted]
        }
        
        
    } else{
        psi_r <- d_inv_treat[first_ind_r:last_ind_r, sel_treat_inds_shifted]
        # Since first_ind_r and last_ind_r are the same, need to coerce this
        # object to be a matrix with one row so that it works smoothly with
        # rbind() later
        d_inv_treat_sel <- matrix(d_inv_treat[first_ind_r:last_ind_r,
            sel_treat_inds_shifted], nrow=1)
    }

    stopifnot(is.matrix(d_inv_treat_sel))

    return(list(psi_r=psi_r, d_inv_treat_sel=d_inv_treat_sel))
}


getTeResults2 <- function(
    # model,
    sig_eps_sq,
    N,
    T,
    R,
    num_treats,
    cohort_tes,
    cohort_probs,
    psi_mat,
    gram_inv,
    sel_treat_inds_shifted,
    tes,
    d_inv_treat_sel,
    cohort_probs_overall,
    first_inds,
    theta_hat_treat_sel,
    calc_ses,
    indep_probs=FALSE
    ){

    att_hat <- as.numeric(cohort_tes %*% cohort_probs)

    if(calc_ses){
        # Get ATT standard error
        # first variance term: convergence of theta
        psi_att <- psi_mat %*% cohort_probs

        att_var_1 <- sig_eps_sq * as.numeric(t(psi_att) %*% gram_inv %*%
            psi_att) / (N * T)

        # Second variance term: convergence of cohort membership probabilities
        att_var_2 <- getSecondVarTermDataApp(
            cohort_probs=cohort_probs,
            psi_mat=psi_mat,
            sel_treat_inds_shifted=sel_treat_inds_shifted,
            tes=tes,
            d_inv_treat_sel=d_inv_treat_sel,
            cohort_probs_overall=cohort_probs_overall,
            first_inds=first_inds,
            theta_hat_treat_sel=theta_hat_treat_sel,
            num_treats=num_treats,
            N=N,
            T=T,
            R=R)

        if(indep_probs){
            att_te_se <- sqrt(att_var_1 + att_var_2)
        } else{
            att_te_se <- sqrt(att_var_1 + att_var_2 + 2*sqrt(
                att_var_1 * att_var_2))
        }

        att_te_se_no_prob <- sqrt(att_var_1)

    } else{
        att_te_se <- NA
        att_te_se_no_prob <- NA
    }

    return(list(
        att_hat=att_hat,
        att_te_se=att_te_se,
        att_te_se_no_prob=att_te_se_no_prob
        )
    )
}