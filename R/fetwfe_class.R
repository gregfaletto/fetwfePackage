#' @title Fused Extended Two-Way Fixed Effects Output Class
#' 
#' @description S3 class for the output of \code{fetwfe()}. Contains methods for printing and summarizing the results.
#' 
#' @param x   An object of class \code{fetwfe}.  (Used by \code{print.fetwfe} and \code{print.summary.fetwfe}.)
#' @param object   An object of class \code{fetwfe}.  (Used by \code{coef.fetwfe} and \code{summary.fetwfe}.)
#' @param show_internal Logical; if TRUE, internal outputs will be shown in the summary. Default is FALSE.
#' @param ... Additional arguments passed to methods
#' 
#' @name fetwfe-class
NULL

#' @rdname fetwfe-class
#' @export
coef.fetwfe <- function(object, ...) {
    return(object$beta_hat)
}

#' @rdname fetwfe-class
#' @export
print.fetwfe <- function(x, show_internal = FALSE, ...) {
    cat("Fused Extended Two-Way Fixed Effects Results\n")
    cat("===========================================\n\n")
    
    # Print overall ATT
    cat("Overall Average Treatment Effect (ATT):\n")
    cat(sprintf("  Estimate: %.4f\n", x$att_hat))
    if (!is.na(x$att_se)) {
        cat(sprintf("  Standard Error: %.4f\n", x$att_se))
    }
    cat("\n")
    
    # Print cohort-specific effects
    cat("Cohort Average Treatment Effects (CATT):\n")
    catt_df <- x$catt_df
    for (i in 1:nrow(catt_df)) {
        cat(sprintf("  Cohort %s:\n", catt_df$Cohort[i]))
        cat(sprintf("    Estimate: %.4f\n", catt_df$`Estimated TE`[i]))
        if (!is.na(catt_df$SE[i])) {
            cat(sprintf("    Standard Error: %.4f\n", catt_df$SE[i]))
            cat(sprintf("    %d%% CI: [%.4f, %.4f]\n", 
                       (1 - x$alpha) * 100,
                       catt_df$ConfIntLow[i],
                       catt_df$ConfIntHigh[i]))
        }
        cat("\n")
    }
    
    # Print model details
    cat("Model Details:\n")
    cat(sprintf("  Number of units (N): %d\n", x$N))
    cat(sprintf("  Number of time periods (T): %d\n", x$T))
    cat(sprintf("  Number of treated cohorts (R): %d\n", x$R))
    cat(sprintf("  Number of covariates (d): %d\n", x$d))
    cat(sprintf("  Total number of features (p): %d\n", x$p))
    cat(sprintf("  Selected model size: %d\n", x$lambda_star_model_size))
    cat(sprintf("  Lambda*: %.4f\n", x$lambda_star))
    
    # Print internal details if requested
    if (show_internal) {
        cat("\nInternal Details:\n")
        cat("  Design matrix dimensions:", dim(x$internal$X_ints), "\n")
        cat("  Response vector length:", length(x$internal$y), "\n")
        cat("  Standard errors calculated:", x$internal$calc_ses, "\n")
    }
    
    invisible(x)
}

#' @rdname fetwfe-class
#' @export
summary.fetwfe <- function(object, ..., show_internal = FALSE) {
    # Create a summary list
    summary_list <- list(
        att = list(
            estimate = object$att_hat,
            se       = object$att_se
        ),
        catt = object$catt_df,
        model_info = list(
            N              = object$N,
            T              = object$T,
            R              = object$R,
            d              = object$d,
            p              = object$p,
            lambda_star    = object$lambda_star,
            model_size     = object$lambda_star_model_size,
            sig_eps_sq     = object$sig_eps_sq,
            sig_eps_c_sq   = object$sig_eps_c_sq
        )
    )

    # Add internal details if requested
    if (show_internal) {
        summary_list$internal <- object$internal
    }

    class(summary_list) <- "summary.fetwfe"
    return(summary_list)
}


#' @rdname fetwfe-class
#' @export
print.summary.fetwfe <- function(x, ...) {
    cat("Summary of Fused Extended Two-Way Fixed Effects Results\n")
    cat("=====================================================\n\n")
    
    # Print overall ATT
    cat("Overall Average Treatment Effect (ATT):\n")
    cat(sprintf("  Estimate: %.4f\n", x$att$estimate))
    if (!is.na(x$att$se)) {
        cat(sprintf("  Standard Error: %.4f\n", x$att$se))
    }
    cat("\n")
    
    # Print cohort-specific effects
    cat("Cohort Average Treatment Effects (CATT):\n")
    catt_df <- x$catt
    for (i in 1:nrow(catt_df)) {
        cat(sprintf("  Cohort %s:\n", catt_df$Cohort[i]))
        cat(sprintf("    Estimate: %.4f\n", catt_df$`Estimated TE`[i]))
        if (!is.na(catt_df$SE[i])) {
            cat(sprintf("    Standard Error: %.4f\n", catt_df$SE[i]))
            cat(sprintf("    %d%% CI: [%.4f, %.4f]\n", 
                       (1 - 0.05) * 100,  # Assuming alpha = 0.05
                       catt_df$ConfIntLow[i],
                       catt_df$ConfIntHigh[i]))
        }
        cat("\n")
    }
    
    # Print model details
    cat("Model Details:\n")
    cat(sprintf("  Number of units (N): %d\n", x$model_info$N))
    cat(sprintf("  Number of time periods (T): %d\n", x$model_info$T))
    cat(sprintf("  Number of treated cohorts (R): %d\n", x$model_info$R))
    cat(sprintf("  Number of covariates (d): %d\n", x$model_info$d))
    cat(sprintf("  Total number of features (p): %d\n", x$model_info$p))
    cat(sprintf("  Selected model size: %d\n", x$model_info$model_size))
    cat(sprintf("  Lambda*: %.4f\n", x$model_info$lambda_star))
    cat(sprintf("  Row-level noise variance (sig_eps_sq): %.4f\n", x$model_info$sig_eps_sq))
    cat(sprintf("  Unit-level noise variance (sig_eps_c_sq): %.4f\n", x$model_info$sig_eps_c_sq))
    
    # Print internal details if present
    if (!is.null(x$internal)) {
        cat("\nInternal Details:\n")
        cat("  Design matrix dimensions:", dim(x$internal$X_ints), "\n")
        cat("  Response vector length:", length(x$internal$y), "\n")
        cat("  Standard errors calculated:", x$internal$calc_ses, "\n")
    }
    
    invisible(x)
} 