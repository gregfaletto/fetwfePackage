# predict() methods for fetwfe / etwfe / betwfe (#33).
#
# Implements the closed-form CATT(r, t, x) estimator from
# `paper_arxiv.tex` Eq. (catt.t.r.est) (around line 803):
#
#   tau_hat_CATT(r, t, x) = tau_hat_rt + (x - X_bar_r)^T rho_hat_rt
#
# along with the variance estimator from Theorem
# `te.asym.norm.thm.gen.cond`(a) (Eq. (v.n.r.t.catt.const), around
# paper line 2884). For each scalar (r, t, x) prediction the variance
# decomposes into two pieces:
#
#   (a) Regression-coefficient noise (`v_reg` below):
#         sigma^2 * psi_hat^T (Gram_S)^{-1} psi_hat / (N * T)
#       where `psi_hat` is the length-|S| selected-features vector with
#       the row of D_N^{-1} corresponding to tau_rt at the tau slot, and
#       (x - X_bar_r)^T applied to the d rows of D_N^{-1} at the rho_rt
#       slot. `Gram_S` is the empirical Gram matrix of the selected
#       columns of `X_final` (= Z * D_N^{-1}).
#
#   (b) Cohort-mean sample noise (`v_mean` below):
#         (1 / N_r) * rho_hat_rt^T Sigma_X|r rho_hat_rt
#       where `Sigma_X|r` is the within-cohort sample covariance of the
#       covariates for cohort r, and `N_r` is the number of units in
#       cohort r. This is the contraction of the second term of the
#       theorem's variance estimator
#         (T / (N * T)) * (rho_hat_rt)^T Sigma_X|r rho_hat_rt / (N_r / N)
#       at psi_{rt} = 1 / 0 (single-(r, t) prediction).
#
# The total prediction variance is then `(v_reg + v_mean)`. We floor any
# floating-point-tiny negatives via `.floor_cluster_quad()` (same pattern
# as the cluster-SE machinery; see `R/cluster_floor.R`).

#' Predict cohort-time treatment effects at user-specified covariate values
#'
#' @description
#' For a fitted [fetwfe()] / [etwfe()] / [betwfe()] object, computes the
#' closed-form conditional average treatment effect on the treated at a
#' user-supplied covariate vector `x`, for every cohort `r` and post-
#' treatment time `t` the model estimated:
#'
#' \deqn{\hat\tau_{\mathrm{CATT}}(r, t, x) = \hat\tau_{rt} +
#'     (x - \bar X_r)^\top \hat\rho_{rt}}
#'
#' (paper Eq. `catt.t.r.est`). Asymptotic standard errors and Wald
#' confidence intervals follow Theorem `te.asym.norm.thm.gen.cond`(a)
#' from the methodology paper (Eq. `v.n.r.t.catt.const`).
#'
#' Output is a long-form data frame with one row per
#' `(newdata_row x cohort x post-treatment time)` combination. When
#' `newdata = NULL` (the default), predictions are computed at each
#' cohort's sample-mean covariate vector \eqn{\bar X_r}, so the
#' `(x - \bar X_r)^\top \hat\rho_{rt}` term vanishes and the
#' `estimate` column reduces to \eqn{\hat\tau_{rt}}.
#'
#' @param object A fitted object from [fetwfe()] (class `"fetwfe"`),
#'   [etwfe()] (class `"etwfe"`), or [betwfe()] (class `"betwfe"`).
#' @param newdata Optional data frame of covariate vectors at which to
#'   predict. Must contain numeric columns whose names match every
#'   covariate the fit used (i.e. the `d` covariates encoded in the
#'   fitted design matrix; see `colnames(object$internal$X_ints)`).
#'   If `NULL` (the default), predictions are computed at each cohort's
#'   sample-mean covariate vector \eqn{\bar X_r}.
#' @param conf.int Logical scalar; if `TRUE` (the default), include
#'   `std.error`, `conf.low`, and `conf.high` columns in the returned
#'   data frame. Requires the fit to have computed standard errors
#'   (i.e. `q < 1` for the bridge estimators or `calc_ses = TRUE` for
#'   ETWFE). When `FALSE`, `std.error`/`conf.low`/`conf.high` are still
#'   present but populated with `NA_real_`.
#' @param conf.level Optional numeric in `(0, 1)`; the confidence level
#'   for the Wald intervals. If `NULL` (the default), uses
#'   `1 - object$alpha` (the alpha encoded at fit time).
#' @param ... Currently unused; present for S3 compatibility.
#' @return A data frame with class
#'   `c("predict_fetwfe", "data.frame")` (or `"predict_etwfe"` /
#'   `"predict_betwfe"` for the other two methods) and columns:
#'   \describe{
#'     \item{`x_row`}{Integer; index into `newdata` (or `1:R` when
#'       `newdata = NULL`, in which case `x_row[r]` corresponds to
#'       cohort `r`'s mean covariate vector).}
#'     \item{`cohort`}{Cohort identifier (character or integer, matching
#'       the fit's `catt_df$cohort` column).}
#'     \item{`time`}{Integer; calendar time of the predicted treatment
#'       effect (a post-treatment period for the cohort).}
#'     \item{`estimate`}{Numeric; \eqn{\hat\tau_{\mathrm{CATT}}(r, t, x)}.}
#'     \item{`std.error`}{Numeric; asymptotic standard error from
#'       Theorem `te.asym.norm.thm.gen.cond`(a). `NA_real_` when
#'       the fit didn't compute standard errors or when `conf.int = FALSE`.}
#'     \item{`conf.low`, `conf.high`}{Numeric; Wald confidence interval
#'       bounds at `conf.level`. `NA_real_` when standard errors aren't
#'       available.}
#'   }
#'
#' @details
#' The variance formula (paper Eq. `v.n.r.t.catt.const`) decomposes into
#' two components for a single-`(r, t, x)` prediction:
#'
#' 1. **Regression-coefficient variance.** The joint sampling
#'    uncertainty in \eqn{(\hat\tau_{rt}, \hat\rho_{rt})} carried
#'    through the selected-features Gram inverse
#'    \eqn{\hat\Sigma((Z D_N^{-1})_{\cdot, \hat S})^{-1}} and the
#'    relevant rows of \eqn{D_N^{-1}} (the row corresponding to
#'    \eqn{\tau_{rt}} and the `d` rows corresponding to \eqn{\rho_{rt}}),
#'    scaled by \eqn{(x - \bar X_r)} on the \eqn{\rho} block.
#' 2. **Cohort-mean variance.** Additional \eqn{O_p(N^{-1/2})} noise
#'    from estimating \eqn{\bar X_r} with the sample mean over the
#'    \eqn{N_r} units in cohort `r`. This term equals
#'    \eqn{\hat\rho_{rt}^\top \widehat{\mathrm{Cov}}(X_i \mid W_i = r)
#'        \hat\rho_{rt} / N_r}.
#'
#' Both components are nonnegative in expectation; floating-point-tiny
#' negatives from numerical cancellation are floored to zero via the
#' same diagnostic helper used for the cluster-robust SE machinery.
#'
#' When `newdata = NULL`, the `(x - \bar X_r)^\top \hat\rho_{rt}` term
#' vanishes by construction, and the regression-variance component
#' collapses to the same expression the package already reports for the
#' per-`(r, t)` cohort treatment effect (component (1)). The
#' cohort-mean variance term remains, reflecting the fact that even at
#' \eqn{x = \bar X_r}, the population \eqn{\bar X_r} is estimated.
#'
#' @section Out of scope:
#' `predict.twfeCovs()` is not provided — `twfeCovs` lacks the
#' bridge-penalty selection that grounds the CATT-prediction theorem.
#' See issue #58 for the broader `twfeCovs` class-method gap.
#'
#' Theorem `te.asym.norm.thm.gen.cond`(b), which incorporates estimated
#' generalized propensity scores, is also deferred — see the pluggable-
#' propensity-scores draft (#32).
#'
#' @seealso [fetwfe()] for the estimator that produces the fitted
#' object, [eventStudy()] / [cohortStudy()] for parallel post-fit
#' aggregations, and [tidy()] / [augment()] for broom-style summaries
#' of the same fit.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   coefs <- genCoefs(R = 3, T = 5, density = 0.5, eff_size = 1, d = 2)
#'   sim <- simulateData(
#'     coefs,
#'     N = 60,
#'     sig_eps_sq = 0.5,
#'     sig_eps_c_sq = 0.5
#'   )
#'   res <- fetwfeWithSimulatedData(sim)
#'
#'   # Predict at each cohort's mean covariate vector (recovers tau_rt):
#'   predict(res)
#'
#'   # Predict at a user-specified x vector:
#'   nd <- data.frame(cov1 = c(0, 1), cov2 = c(-0.5, 0.5))
#'   pred <- predict(res, newdata = nd)
#'   head(pred)
#' }
#'
#' @export
predict.fetwfe <- function(
	object,
	newdata = NULL,
	conf.int = TRUE,
	conf.level = NULL,
	...
) {
	.predict_estimator(
		object = object,
		newdata = newdata,
		conf.int = conf.int,
		conf.level = conf.level,
		out_class = "predict_fetwfe"
	)
}

#' @rdname predict.fetwfe
#' @export
predict.etwfe <- function(
	object,
	newdata = NULL,
	conf.int = TRUE,
	conf.level = NULL,
	...
) {
	.predict_estimator(
		object = object,
		newdata = newdata,
		conf.int = conf.int,
		conf.level = conf.level,
		out_class = "predict_etwfe"
	)
}

#' @rdname predict.fetwfe
#' @export
predict.betwfe <- function(
	object,
	newdata = NULL,
	conf.int = TRUE,
	conf.level = NULL,
	...
) {
	.predict_estimator(
		object = object,
		newdata = newdata,
		conf.int = conf.int,
		conf.level = conf.level,
		out_class = "predict_betwfe"
	)
}

# .predict_estimator
#
# Shared body for `predict.fetwfe()` / `.etwfe()` / `.betwfe()`.
# Handles input validation, the closed-form point estimate, the
# variance computation per Theorem `te.asym.norm.thm.gen.cond`(a), and
# assembly of the long-form output.
#
# @param object A fitted `fetwfe` / `etwfe` / `betwfe` object.
# @param newdata Optional data frame of new covariate vectors.
# @param conf.int Logical; whether to attach CI columns.
# @param conf.level Numeric in (0, 1); overrides `object$alpha`.
# @param out_class Character scalar; the leading S3 class to attach.
# @return A long-form data frame as documented for `predict.fetwfe()`.
# @keywords internal
# @noRd
.predict_estimator <- function(
	object,
	newdata,
	conf.int,
	conf.level,
	out_class
) {
	# Class precondition (each public predict.* enforces this via dispatch,
	# but the helper double-checks for direct callers).
	if (!inherits(object, c("fetwfe", "etwfe", "betwfe"))) {
		stop(
			".predict_estimator() requires an object of class ",
			"'fetwfe', 'etwfe', or 'betwfe'.",
			call. = FALSE
		)
	}

	if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
		stop("`conf.int` must be a single non-NA logical value.", call. = FALSE)
	}

	if (!is.null(conf.level)) {
		if (
			!is.numeric(conf.level) ||
				length(conf.level) != 1L ||
				is.na(conf.level) ||
				conf.level <= 0 ||
				conf.level >= 1
		) {
			stop(
				"`conf.level` must be a single number strictly between 0 and 1.",
				call. = FALSE
			)
		}
		alpha <- 1 - conf.level
	} else {
		alpha <- object$alpha
	}

	# Resolve top-level vs internal slots (fetwfe stores X_ints/X_final
	# only under $internal; etwfe/betwfe duplicate them at the top level).
	bundle <- .predict_extract_internals(object)
	X_ints <- bundle$X_ints
	X_final <- bundle$X_final
	N <- bundle$N
	T <- bundle$T
	R <- bundle$R
	d <- bundle$d
	p <- bundle$p
	calc_ses <- bundle$calc_ses
	sig_eps_sq <- bundle$sig_eps_sq

	beta_hat <- object$beta_hat
	treat_inds <- object$treat_inds
	treat_int_inds <- object$treat_int_inds
	c_names <- as.character(object$catt_df$cohort)
	first_inds <- getFirstInds(R = R, T = T)
	num_treats <- length(treat_inds)

	# Validate newdata shape and resolve cohort-mean newdata default.
	model_cov_names <- if (d > 0) {
		colnames(X_ints)[(R + T - 1 + 1):(R + T - 1 + d)]
	} else {
		character(0)
	}

	cohort_means <- .predict_cohort_means(
		X_ints = X_ints,
		R = R,
		N = N,
		T = T,
		d = d,
		model_cov_names = model_cov_names
	)

	if (is.null(newdata)) {
		# Default: predict at each cohort's mean X_bar_r.
		newdata_mat <- cohort_means # R x d
		x_row_for <- function(r_idx) r_idx
	} else {
		if (!is.data.frame(newdata)) {
			stop("`newdata` must be a data frame.", call. = FALSE)
		}
		if (d > 0) {
			missing_cov <- setdiff(model_cov_names, colnames(newdata))
			if (length(missing_cov) > 0L) {
				stop(
					"`newdata` is missing required covariate column(s): ",
					paste(missing_cov, collapse = ", "),
					". The fitted model used covariates: ",
					paste(model_cov_names, collapse = ", "),
					".",
					call. = FALSE
				)
			}
			# Subset & coerce to numeric matrix; error on non-numeric / NA.
			newdata_sub <- newdata[, model_cov_names, drop = FALSE]
			for (cname in model_cov_names) {
				if (!is.numeric(newdata_sub[[cname]])) {
					stop(
						"`newdata$",
						cname,
						"` must be numeric. Got class: ",
						paste(class(newdata_sub[[cname]]), collapse = "/"),
						".",
						call. = FALSE
					)
				}
				if (any(is.na(newdata_sub[[cname]]))) {
					stop(
						"`newdata$",
						cname,
						"` contains NA; cannot predict at a missing covariate value.",
						call. = FALSE
					)
				}
			}
			newdata_mat <- as.matrix(newdata_sub)
		} else {
			# d == 0: no covariate columns to validate. Use nrow(newdata)
			# rows of "empty" features; predictions reduce to tau_hat_rt.
			newdata_mat <- matrix(numeric(0), nrow = nrow(newdata), ncol = 0L)
		}
		x_row_for <- function(r_idx) integer(0) # placeholder, overridden below
	}

	n_rows_pred <- nrow(newdata_mat)
	if (n_rows_pred == 0L) {
		stop(
			"`newdata` must have at least one row to predict for.",
			call. = FALSE
		)
	}

	# Build the cohort-time index table that drives the rest of the loop.
	# cohort_time[i, ] = (r_idx in 1:R, t_calendar in {cohort_start_time:T}).
	# `cohort_start_times` is decoded from the cohort labels in catt_df;
	# the calendar time for the k-th treatment effect in cohort r is
	# cohort_start_times[r] + (k - 1).
	cohort_start_times <- suppressWarnings(as.integer(c_names))
	if (any(is.na(cohort_start_times))) {
		# Fall back to 1-based cohort labeling when the cohort identifier
		# isn't an integer-coercible string (e.g., real-data cohorts that
		# happen to be all-numeric strings are handled by the cast above;
		# non-numeric labels lose their calendar-time interpretation).
		# We still report a `time` column relative to the cohort's first
		# treatment period (1 = first, 2 = second, ...).
		cohort_start_times <- rep(1L, R)
	}

	# tes_rt[k] is tau_hat_rt for the k-th cohort/time entry (k in 1..num_treats),
	# in the same ordering as `treat_inds`.
	tes_rt <- beta_hat[treat_inds]
	# rho_mat[k, ] is rho_hat_rt for the k-th cohort/time entry (length d).
	if (d > 0L) {
		rho_mat <- matrix(
			beta_hat[treat_int_inds],
			nrow = num_treats,
			ncol = d,
			byrow = TRUE
		)
	} else {
		rho_mat <- matrix(numeric(0), nrow = num_treats, ncol = 0L)
	}

	# Per-treatment-effect (r, t) decoding.
	rt_table <- .predict_rt_table(
		R = R,
		num_treats = num_treats,
		first_inds = first_inds,
		cohort_start_times = cohort_start_times
	)
	stopifnot(nrow(rt_table) == num_treats)

	# Within-cohort sample covariance + cohort sizes (N_r), used for the
	# cohort-mean variance term. Both are computed once outside the loop.
	cohort_info <- .predict_cohort_info(
		X_ints = X_ints,
		R = R,
		N = N,
		T = T,
		d = d,
		model_cov_names = model_cov_names
	)
	cov_X_given_r <- cohort_info$cov_X_given_r # list of length R, each d x d
	N_r_vec <- cohort_info$N_r_vec # length R

	# Set up the regression-coefficient-variance machinery (the
	# Gram-inverse of the selected columns of X_final). Skipped when
	# calc_ses is FALSE or conf.int is FALSE.
	want_ses <- isTRUE(conf.int) && isTRUE(calc_ses)

	if (want_ses) {
		# Resolve selected-feature indices.
		sel_info <- .predict_selected_indices(object, p = p)
		sel_feat_inds <- sel_info$sel_feat_inds
		sel_treat_block_inds <- sel_info$sel_treat_block_inds # in 1..num_treats
		sel_treat_int_block_inds <- sel_info$sel_treat_int_block_inds # in 1..(d * num_treats)

		# Centered Gram of selected columns of X_final.
		gram_full_inv <- .predict_gram_inv(
			X_final = X_final,
			sel_feat_inds = sel_feat_inds,
			N = N,
			T = T
		)
		if (is.null(gram_full_inv)) {
			# Singular gram; fall back to NA SEs to avoid hard failure
			# in well-specified-but-poorly-conditioned cases.
			want_ses <- FALSE
		}
	} else {
		gram_full_inv <- NULL
		sel_feat_inds <- integer(0)
		sel_treat_block_inds <- integer(0)
		sel_treat_int_block_inds <- integer(0)
	}

	# Build D_N^{-1} restricted to the treatment-effect block. The full
	# D_N^{-1} is block-diagonal; for the variance formula we need the
	# row of (D^{(2)}(R))^{-1} that gives tau_rt back from the treatment
	# slot of theta, and the same matrix tiled by the covariate index
	# for the rho_rt rows. For unfused (etwfe), D^{(2)}(R)^{-1} is the
	# identity, so the relevant rows are just standard basis vectors.
	if (.predict_is_fused(object)) {
		d_inv_treat_full <- genInvTwoWayFusionTransformMat(
			num_treats,
			first_inds,
			R
		)
	} else {
		d_inv_treat_full <- diag(num_treats)
	}

	# Assemble output, row by (newdata_row x cohort x time).
	out_list <- vector("list", n_rows_pred * num_treats)
	out_idx <- 0L

	for (i in seq_len(n_rows_pred)) {
		x_i <- as.numeric(newdata_mat[i, ])

		# When newdata is NULL, x_i is X_bar_r for cohort i, so the
		# (x - X_bar_r) term vanishes for cohort i but is nonzero for
		# other cohorts r' != i. Loop over all treatment effects (r, t).
		for (k in seq_len(num_treats)) {
			out_idx <- out_idx + 1L
			r_idx <- rt_table$r_idx[k]
			t_cal <- rt_table$t_cal[k]
			x_bar_r <- cohort_means[r_idx, ]
			if (d > 0L) {
				x_diff <- x_i - as.numeric(x_bar_r)
				rho_k <- rho_mat[k, ]
				point <- as.numeric(tes_rt[k] + sum(x_diff * rho_k))
			} else {
				x_diff <- numeric(0)
				rho_k <- numeric(0)
				point <- as.numeric(tes_rt[k])
			}

			if (want_ses) {
				v_reg <- .predict_v_reg(
					k = k,
					d = d,
					num_treats = num_treats,
					sel_treat_block_inds = sel_treat_block_inds,
					sel_treat_int_block_inds = sel_treat_int_block_inds,
					x_diff = x_diff,
					d_inv_treat_full = d_inv_treat_full,
					gram_full_inv = gram_full_inv,
					sig_eps_sq = sig_eps_sq,
					N = N,
					T = T
				)
				v_mean <- .predict_v_mean(
					rho_k = rho_k,
					cov_X_r = cov_X_given_r[[r_idx]],
					N_r = N_r_vec[r_idx]
				)
				var_total <- .floor_cluster_quad(
					v_reg + v_mean,
					"predict.fetwfe/var_total"
				)
				se_k <- sqrt(var_total)
				z <- stats::qnorm(1 - alpha / 2)
				ci_lo <- point - z * se_k
				ci_hi <- point + z * se_k
			} else {
				se_k <- NA_real_
				ci_lo <- NA_real_
				ci_hi <- NA_real_
			}

			out_list[[out_idx]] <- list(
				x_row = i,
				cohort = c_names[r_idx],
				time = t_cal,
				estimate = point,
				std.error = se_k,
				conf.low = ci_lo,
				conf.high = ci_hi
			)
		}
	}

	out_df <- data.frame(
		x_row = vapply(out_list, function(z) z$x_row, integer(1)),
		cohort = vapply(out_list, function(z) z$cohort, character(1)),
		time = vapply(out_list, function(z) z$time, integer(1)),
		estimate = vapply(out_list, function(z) z$estimate, numeric(1)),
		std.error = vapply(out_list, function(z) z$std.error, numeric(1)),
		conf.low = vapply(out_list, function(z) z$conf.low, numeric(1)),
		conf.high = vapply(out_list, function(z) z$conf.high, numeric(1)),
		stringsAsFactors = FALSE
	)

	class(out_df) <- c(out_class, "data.frame")
	out_df
}

# .predict_extract_internals
#
# Resolve where the design-matrix and noise-variance slots live on a
# fitted object. fetwfe stores everything under $internal; etwfe and
# betwfe expose them both at the top level AND under $internal (#154).
# Read $internal preferentially because that's the canonical path per
# the slot-validators in R/{fetwfe,etwfe,betwfe}_class.R.
#
# @keywords internal
# @noRd
.predict_extract_internals <- function(object) {
	internal <- object$internal
	if (is.null(internal)) {
		stop(
			"Fitted object is missing the `$internal` slot; cannot ",
			"compute predictions. Was the object constructed by ",
			"fetwfe() / etwfe() / betwfe()?",
			call. = FALSE
		)
	}
	list(
		X_ints = internal$X_ints,
		X_final = internal$X_final,
		N = object$N,
		T = object$T,
		R = object$R,
		d = object$d,
		p = object$p,
		calc_ses = isTRUE(internal$calc_ses) ||
			(is.null(internal$calc_ses) && isTRUE(object$calc_ses)),
		sig_eps_sq = object$sig_eps_sq
	)
}

# .predict_cohort_means
#
# Compute the per-cohort sample mean X_bar_r of each model covariate.
# Reads from the d "main effect" columns of `X_ints` (which carry the
# original, uncentered covariate values, repeated T times per unit),
# subsetted to the rows belonging to cohort r via the cohort-dummy
# columns 1:R of X_ints. To get one row per unit, we take the rows at
# t = 1 (every T-th row starting from 1).
#
# @return A numeric matrix R x d, rownames set to the cohort labels.
# @keywords internal
# @noRd
.predict_cohort_means <- function(X_ints, R, N, T, d, model_cov_names) {
	if (d == 0L) {
		return(matrix(numeric(0), nrow = R, ncol = 0L))
	}
	cov_cols <- (R + T - 1 + 1):(R + T - 1 + d)
	# Take t = 1 rows for each unit so we have one observation per unit.
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe_t1 <- X_ints[t1_rows, 1:R, drop = FALSE]
	X_t1 <- X_ints[t1_rows, cov_cols, drop = FALSE]
	# Per-unit cohort (1..R, or 0 for never-treated).
	unit_cohort <- apply(cohort_fe_t1, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})
	out <- matrix(NA_real_, nrow = R, ncol = d)
	colnames(out) <- model_cov_names
	for (r in seq_len(R)) {
		rows_r <- which(unit_cohort == r)
		if (length(rows_r) == 0L) {
			stop(
				"Cohort ",
				r,
				" has zero units in the fitted data. ",
				"Cannot compute X_bar_r; refit on a panel with at least ",
				"one unit per cohort.",
				call. = FALSE
			)
		}
		out[r, ] <- colMeans(X_t1[rows_r, , drop = FALSE])
	}
	out
}

# .predict_cohort_info
#
# Compute the within-cohort sample covariance Cov(X_i | W_i = r) and
# cohort size N_r for each cohort r. Used for the cohort-mean variance
# component of the CATT prediction variance.
#
# @return A list with elements:
#   - cov_X_given_r: list of length R; r-th element is a d x d matrix.
#   - N_r_vec: integer vector of length R.
# @keywords internal
# @noRd
.predict_cohort_info <- function(X_ints, R, N, T, d, model_cov_names) {
	cov_list <- vector("list", R)
	N_r_vec <- integer(R)

	if (d == 0L) {
		for (r in seq_len(R)) {
			cov_list[[r]] <- matrix(numeric(0), nrow = 0L, ncol = 0L)
		}
	}

	cov_cols <- if (d > 0L) {
		(R + T - 1 + 1):(R + T - 1 + d)
	} else {
		integer(0)
	}
	t1_rows <- seq(1L, N * T, by = T)
	cohort_fe_t1 <- X_ints[t1_rows, 1:R, drop = FALSE]
	X_t1 <- if (d > 0L) {
		X_ints[t1_rows, cov_cols, drop = FALSE]
	} else {
		NULL
	}
	unit_cohort <- apply(cohort_fe_t1, 1L, function(row) {
		if (all(row == 0)) 0L else which(row == 1)[1L]
	})

	for (r in seq_len(R)) {
		rows_r <- which(unit_cohort == r)
		N_r_vec[r] <- length(rows_r)
		if (d == 0L) {
			next # already filled with 0 x 0 matrix
		}
		X_r <- X_t1[rows_r, , drop = FALSE]
		if (length(rows_r) >= 2L) {
			cov_list[[r]] <- stats::cov(X_r) # uses N_r - 1 in the denominator
		} else {
			# Single-unit cohort: stats::cov() returns NA matrix. The
			# variance contribution from cohort-mean noise is 0 in this
			# degenerate case (no within-cohort variance to propagate);
			# the regression-variance component still fires.
			cov_list[[r]] <- matrix(0, nrow = d, ncol = d)
		}
	}

	list(cov_X_given_r = cov_list, N_r_vec = N_r_vec)
}

# .predict_rt_table
#
# Decode the treatment-effect index k in {1, ..., num_treats} into the
# corresponding (cohort_index_r, calendar_time_t) tuple. The ordering
# matches `treat_inds` / `first_inds` / the columns of `treat_mat_long`
# in R/design_matrix.R.
#
# @return A data frame with `num_treats` rows and columns r_idx (1..R)
#   and t_cal (calendar time, derived from cohort_start_times[r]).
# @keywords internal
# @noRd
.predict_rt_table <- function(R, num_treats, first_inds, cohort_start_times) {
	r_idx <- integer(num_treats)
	t_cal <- integer(num_treats)
	for (r in seq_len(R)) {
		start <- first_inds[r]
		end <- if (r < R) first_inds[r + 1L] - 1L else num_treats
		n_t_r <- end - start + 1L
		r_idx[start:end] <- r
		# k-th treatment effect for cohort r is at calendar time
		# cohort_start_times[r] + (k - 1), where k = 1..n_t_r.
		t_cal[start:end] <- cohort_start_times[r] + (0L:(n_t_r - 1L))
	}
	data.frame(r_idx = r_idx, t_cal = t_cal)
}

# .predict_is_fused
#
# Heuristic for "did the estimator apply the bridge-fusion penalty?"
# True for fetwfe and betwfe (both use the two-way fusion transform);
# False for etwfe (pure OLS, no fusion in the basis).
#
# @keywords internal
# @noRd
.predict_is_fused <- function(object) {
	inherits(object, c("fetwfe", "betwfe"))
}

# .predict_selected_indices
#
# Recover the selected-features indices on the **transformed** scale
# (i.e., among the columns of X_final = Z * D_N^{-1}). For fetwfe, the
# canonical source is the non-zero entries of `internal$theta_hat`.
# For etwfe / betwfe (which store `X_final` and the full beta_hat but
# not theta_hat), we fall back to "all features selected" — etwfe is
# pure OLS, so every column is in the selected support, and betwfe is
# handled by checking which entries of beta_hat are non-zero (the
# bridge penalty applied to the original basis).
#
# @return A list with sel_feat_inds (1..p), sel_treat_block_inds
#   (subset of 1..num_treats — the offsets within the tau block of
#   selected tau_rt's), and sel_treat_int_block_inds (subset of
#   1..(d * num_treats) — the offsets within the rho block of
#   selected rho_rt entries).
# @keywords internal
# @noRd
.predict_selected_indices <- function(object, p) {
	treat_inds <- object$treat_inds
	treat_int_inds <- object$treat_int_inds
	num_treats <- length(treat_inds)
	d <- object$d
	internal <- object$internal

	# Prefer theta_hat (transformed scale) if available.
	if (!is.null(internal$theta_hat)) {
		theta_slopes <- internal$theta_hat[-1L]
		stopifnot(length(theta_slopes) == p)
		sel_feat_inds <- which(theta_slopes != 0)
	} else {
		# Fall back to "all features selected" for OLS-family estimators
		# (etwfe). On betwfe, theta_hat lives in `internal$theta_hat` per
		# the slot validators, so this branch shouldn't fire.
		sel_feat_inds <- seq_len(p)
	}

	# Indices within the tau block (positions in `treat_inds` that were
	# selected).
	in_tau <- sel_feat_inds %in% treat_inds
	# Map each selected tau-feature index back to its 1..num_treats
	# offset within the tau block. `treat_inds` is sequential so we can
	# just subtract `treat_inds[1] - 1`.
	sel_treat_block_inds <- sel_feat_inds[in_tau] - (treat_inds[1L] - 1L)

	# Same for the rho block (treat_int_inds).
	if (d > 0L && length(treat_int_inds) > 0L) {
		in_rho <- sel_feat_inds %in% treat_int_inds
		sel_treat_int_block_inds <- sel_feat_inds[in_rho] -
			(treat_int_inds[1L] - 1L)
	} else {
		sel_treat_int_block_inds <- integer(0)
	}

	list(
		sel_feat_inds = sel_feat_inds,
		sel_treat_block_inds = sel_treat_block_inds,
		sel_treat_int_block_inds = sel_treat_int_block_inds
	)
}

# .predict_gram_inv
#
# Compute the inverse of the centered Gram matrix of the selected
# columns of X_final. Returns NULL if the matrix is too close to
# singular. Otherwise returns a length(sel_feat_inds) x length(...)
# matrix.
#
# Note: Unlike `getGramInv()` in R/core_funcs.R (which subsets the
# result to *just* the treatment-effect rows/cols), we keep the FULL
# selected-feature inverse here because the CATT variance formula
# needs the entries at *both* the tau-block and the rho-block of S.
#
# @keywords internal
# @noRd
.predict_gram_inv <- function(X_final, sel_feat_inds, N, T) {
	if (length(sel_feat_inds) == 0L) {
		# No selected features (no variance contribution from regression).
		# Return an empty 0 x 0 matrix; downstream code multiplies by an
		# empty psi vector and yields 0.
		return(matrix(numeric(0), nrow = 0L, ncol = 0L))
	}
	X_sel <- X_final[, sel_feat_inds, drop = FALSE]
	X_sel_centered <- scale(X_sel, center = TRUE, scale = FALSE)
	gram <- (1 / (N * T)) * crossprod(X_sel_centered)
	min_ev <- min(eigen(gram, symmetric = TRUE, only.values = TRUE)$values)
	if (min_ev < 1e-16) {
		warning(
			"In predict(): Gram matrix on selected features is not ",
			"invertible. Standard errors will be reported as NA.",
			call. = FALSE
		)
		return(NULL)
	}
	solve(gram)
}

# .predict_v_reg
#
# Compute the regression-coefficient variance contribution for a
# single (r, t, x) prediction:
#   v_reg = sigma^2 * psi_hat^T (Gram_S)^{-1} psi_hat / (N * T)
#
# psi_hat is length |S| and has nonzero entries at:
#   - The tau slot, if tau_rt's row of D^{(2)}(R)^{-1} survives the
#     selection mask. For a single (r, t) prediction (psi_{rt} = 1),
#     this is the k-th row of d_inv_treat_full restricted to
#     `sel_treat_block_inds`.
#   - The rho slots, if (some of) the d rho_rt entries survive. The d
#     rows of D_N^{-1} corresponding to rho_rt (block-diagonal layout,
#     I_d kron D^{(2)}(R)^{-1}) get scaled by (x - X_bar_r) and
#     summed.
#
# @param k Integer; the (r, t) index in 1..num_treats.
# @param d Integer; number of covariates.
# @param num_treats Integer; total number of treatment effects.
# @param sel_treat_block_inds Integer vector; selected tau offsets.
# @param sel_treat_int_block_inds Integer vector; selected rho offsets
#   (1..(d * num_treats)).
# @param x_diff Numeric vector of length d; (x - X_bar_r).
# @param d_inv_treat_full Numeric matrix; num_treats x num_treats; the
#   inverse two-way fusion-transform matrix on the treatment block.
# @param gram_full_inv Numeric matrix; |S| x |S|; inverse of the
#   centered Gram of selected X_final columns.
# @param sig_eps_sq Numeric scalar; idiosyncratic noise variance.
# @param N, T Integers.
# @return A nonnegative numeric scalar.
# @keywords internal
# @noRd
.predict_v_reg <- function(
	k,
	d,
	num_treats,
	sel_treat_block_inds,
	sel_treat_int_block_inds,
	x_diff,
	d_inv_treat_full,
	gram_full_inv,
	sig_eps_sq,
	N,
	T
) {
	S_total <- nrow(gram_full_inv)
	if (S_total == 0L) {
		return(0)
	}
	# psi_hat is a length-|S| vector. We need to know the layout of |S|:
	# in `gram_full_inv`, the rows/cols are ordered by `sel_feat_inds`.
	# Within that, the tau-block selected entries come first (those with
	# original index in `treat_inds`), then the rho-block selected
	# entries. Since `sel_feat_inds` from .predict_selected_indices is
	# sorted, the tau-block entries appear before the rho-block entries.
	#
	# psi_hat decomposition:
	#   - At position j in the tau slot (j in 1..length(sel_treat_block_inds))
	#     the entry is d_inv_treat_full[k, sel_treat_block_inds[j]].
	#   - At position j in the rho slot
	#     (j in 1..length(sel_treat_int_block_inds))
	#     the entry is the contraction of x_diff against the d rows of
	#     D_N^{-1} that correspond to rho_rt at index k. Specifically,
	#     the rho block of D_N^{-1} is I_d kron D^{(2)}(R)^{-1}, so the
	#     j-th column in the rho block has index (m - 1) * num_treats + l
	#     where m in 1..d is the covariate slot and l in 1..num_treats is
	#     the tau slot. Mapping: index = (m - 1) * num_treats + l. The
	#     contribution to psi_hat for rho_rt at index k is, for each m,
	#     x_diff[m] * d_inv_treat_full[k, l] applied at the rho-block
	#     position (m - 1) * num_treats + l (if selected).

	psi <- numeric(S_total)

	# Tau-block contribution.
	n_tau_sel <- length(sel_treat_block_inds)
	if (n_tau_sel > 0L) {
		psi[seq_len(n_tau_sel)] <- d_inv_treat_full[k, sel_treat_block_inds]
	}

	# Rho-block contribution.
	n_rho_sel <- length(sel_treat_int_block_inds)
	if (n_rho_sel > 0L && d > 0L) {
		# Decode each selected rho offset into (m, l).
		# offset = (m - 1) * num_treats + l, with l in 1..num_treats and
		# m in 1..d.
		offset <- sel_treat_int_block_inds
		l_vec <- ((offset - 1L) %% num_treats) + 1L
		m_vec <- ((offset - 1L) %/% num_treats) + 1L
		# psi entry = x_diff[m] * d_inv_treat_full[k, l]
		psi[n_tau_sel + seq_along(offset)] <-
			x_diff[m_vec] * d_inv_treat_full[k, l_vec]
	}

	quad <- as.numeric(t(psi) %*% gram_full_inv %*% psi)
	sig_eps_sq * quad / (N * T)
}

# .predict_v_mean
#
# Compute the cohort-mean variance contribution for a single (r, t, x)
# prediction:
#   v_mean = rho_rt^T Cov(X_i | W_i = r) rho_rt / N_r
#
# This is the contraction of the second term of the theorem's variance
# estimator
#   (T / (N * T)) * (psi_r^CATT)^T Cov(X_i | W_i = r) psi_r^CATT / (N_r / N)
# at psi_{r't'} = 1 if (r, t) = (r', t') and 0 otherwise, with the
# resulting v / (N * T) prefactor handled outside (the theorem's
# v_N / (N * T) gives the variance of the test statistic; this helper
# returns the analogous quantity directly).
#
# Algebra: psi_r^CATT = rho_rt; the prefactor T / (N * T) = 1 / N;
# divided by N_r / N gives 1 / N_r.
#
# @param rho_k Numeric vector length d; rho_hat_rt.
# @param cov_X_r Numeric matrix d x d; Cov(X_i | W_i = r).
# @param N_r Integer; number of units in cohort r.
# @return A nonnegative numeric scalar.
# @keywords internal
# @noRd
.predict_v_mean <- function(rho_k, cov_X_r, N_r) {
	if (length(rho_k) == 0L) {
		return(0)
	}
	if (N_r <= 0L) {
		return(0)
	}
	quad <- as.numeric(t(rho_k) %*% cov_X_r %*% rho_k)
	quad / N_r
}
