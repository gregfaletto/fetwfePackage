# Fusion-transform matrix generators used by FETWFE: forward / backward
# transforms applied to the treatment-effect block, their inverses, and
# the full-coefficient analogues used during post-selection
# untransformation. Moved from R/fetwfe_core.R (and one helper from
# R/core_funcs.R) in 1.9.25.

#-------------------------------------------------------------------------------
# Transformation Matrices and Coefficient Transformations
#-------------------------------------------------------------------------------

#' Transform a Full Design Matrix by  \( \boldsymbol D_N^{-1} \)
#'
#' @description
#' Takes the *raw* stacked panel design matrix
#' \eqn{\tilde{\boldsymbol Z}\in\mathbb G^{NT\times p}}
#' and post-multiplies it by the block-diagonal inverse fusion matrix
#' \(\boldsymbol D_N^{-1}\) from Lemma 3:
#' \deqn{
#'   \boldsymbol D_N^{-1}
#'   = \operatorname{diag}\!\bigl(
#'       (D^{(1)}(G))^{-1},\;
#'       (D^{(1)}(T-1))^{-1},\;
#'       I_{d},\;
#'       (D^{(1)}(G))^{-1}\otimes I_{d},\;
#'       (D^{(1)}(T-1))^{-1}\otimes I_{d},\;
#'       (D^{(2)}(\mathcal G))^{-1},\;
#'       (D^{(2)}(\mathcal G))^{-1}\otimes I_{d}
#'     \bigr).
#' }
#' The result is a matrix ready for vanilla bridge regression (Lasso,
#' elastic-net, \eqn{\ell_q} etc.), where the penalty on the transformed
#' coefficients reproduces the complex fusion penalty on the original ones.
#'
#' @details
#' The columns of `X_int` **must** appear in the order:
#' \enumerate{
#'   \item Cohort fixed-effects (length \eqn{G})
#'   \item Time fixed-effects (length \eqn{T-1})
#'   \item Main covariates (length \eqn{d})
#'   \item Covariate \(\times\) cohort interactions (length \eqn{dG})
#'   \item Covariate \(\times\) time interactions (length \eqn{d(T-1)})
#'   \item Base treatment effects (length `num_treats`)
#'   \item Covariate \(\times\) treatment interactions (length `d*num_treats`)
#' }
#' Each block is transformed exactly by the corresponding diagonal block of
#' \(\boldsymbol D_N^{-1}\) using helper functions
#' `genBackwardsInvFusionTransformMat()` and
#' `genInvTwoWayFusionTransformMat()`.
#'
#' Side-effect `stopifnot()` guards verify both the expected column order and
#' the absence of `NA`s after transformation.
#'
#' @param X_int Numeric matrix \eqn{NT\times p}.  Original design matrix.
#' @param N,T,G Integers. Panel dimensions and number of treated cohorts.
#' @param d Integer. Number of time-invariant covariates.
#' @param num_treats Integer. Total number of base treatment-effect dummies
#'   \eqn{\mathfrak W}.
#' @param first_inds Optional integer vector of length \eqn{G}.
#'   Starting column indices of the first treatment dummy of each cohort inside
#'   the treatment block.  If `NA` (default) they are computed by
#'   `getFirstInds()`.
#' @param fusion_structure Character; one of `"cohort"` (default) or
#'   `"event_study"`. Selects which built-in inverse treatment-effect fusion
#'   block is used (ignored when `d_inv_treat` is supplied).
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   user-supplied already-inverted treatment-effect fusion block
#'   (`solve(fusion_matrix)`, #236). When non-`NULL` it overrides
#'   `fusion_structure` for the treatment block. Default `NULL`.
#'
#' @return
#' A numeric matrix `X_mod` with the **same dimensions** as `X_int` but whose
#' columns are the transformed regressors
#' \eqn{[\tilde{\boldsymbol Z}\,\boldsymbol D_N^{-1}]_{NT\times p}}.
#'
#' @seealso
#' * `genBackwardsInvFusionTransformMat()`
#' * `genInvTwoWayFusionTransformMat()`
#'
#' @examples
#' set.seed(1)
#' G <- 2; T <- 5; d <- 1; N <- 10
#' num_treats  <- getNumTreats(G, T)
#' p <- G + (T-1) + d + d*G + d*(T-1) + num_treats + d*num_treats
#' X_int <- matrix(rnorm(N*T*p), N*T, p)
#' X_mod <- transformXintImproved(
#'   X_int, N=N, T=T, G=G, d=d, num_treats=num_treats
#' )
#' # The two matrices have identical dimensions:
#' dim(X_mod)  # NT x p
#' @keywords internal
#' @noRd
transformXintImproved <- function(
	X_int,
	N,
	T,
	G,
	d,
	num_treats,
	first_inds = NA,
	fusion_structure = "cohort",
	d_inv_treat = NULL
) {
	p <- getP(G = G, T = T, d = d, num_treats = num_treats)
	stopifnot(p == ncol(X_int))
	X_mod <- matrix(as.numeric(NA), nrow = N * T, ncol = p)
	stopifnot(nrow(X_int) == N * T)

	# Transform cohort fixed effects
	X_mod[, 1:G] <- X_int[, 1:G] %*% genBackwardsInvFusionTransformMat(G)

	# Transform time fixed effects
	X_mod[, (G + 1):(G + T - 1)] <- X_int[, (G + 1):(G + T - 1)] %*%
		genBackwardsInvFusionTransformMat(T - 1)

	# Copy X (the main covariate block; may be empty when d==0)
	if (d > 0) {
		stopifnot(all(is.na(X_mod[, (G + T - 1 + 1):(G + T - 1 + d)])))
		X_mod[, (G + T - 1 + 1):(G + T - 1 + d)] <- X_int[,
			(G + T - 1 + 1):(G + T - 1 + d)
		]

		stopifnot(all(!is.na(X_mod[, 1:(G + T - 1 + d)])))
		stopifnot(all(is.na(X_mod[, (G + T - 1 + d + 1):p])))
	}

	# For cohort effects interacted with X: we have d*G columns to deal with.
	# For each individual feature, this will be handled using
	# genTransformedMatFusion.
	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(G = G, T = T)
	}

	if (d > 0) {
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on G times
			feat_1 <- G + T - 1 + d + j
			feat_R <- G + T - 1 + d + (G - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == G)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))

			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(G)
		}
		stopifnot(all(!is.na(X_mod[, 1:(G + T - 1 + d + G * d)])))
		stopifnot(all(is.na(X_mod[, (G + T - 1 + d + G * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- G + T - 1 + d + G * d + j
			feat_T_minus_1 <- G + T - 1 + d + G * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(X_mod[, feat_inds_j])))
			X_mod[, feat_inds_j] <- X_int[, feat_inds_j] %*%
				genBackwardsInvFusionTransformMat(T - 1)
		}
		stopifnot(all(!is.na(X_mod[, 1:(G + T - 1 + d + G * d + (T - 1) * d)])))
		stopifnot(all(is.na(X_mod[,
			(G + T - 1 + d + G * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects. For each cohort, will penalize base term, then
	# fuse remaining terms toward it. Also, for each cohort, will penalize base
	# treatment effect of this cohort to base of previous cohort. New function
	# genTransformedMatTwoWayFusion does this.

	feat_inds <- (G + T - 1 + d + G * d + (T - 1) * d + 1):(G +
		T -
		1 +
		d +
		G * d +
		(T - 1) * d +
		num_treats)

	# Now ready to generate the appropriate transformed matrix
	stopifnot(all(is.na(X_mod[, feat_inds])))

	X_mod[, feat_inds] <- X_int[, feat_inds] %*%
		.gen_inv_treat_block(
			num_treats,
			first_inds,
			G,
			fusion_structure,
			d_inv_treat = d_inv_treat
		)

	if (d > 0) {
		# Lastly, penalize interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with genTransformedMatTwoWayFusion, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + G + T - 1 + d + G * d + (T - 1) * d + num_treats

			# Now ready to generate the appropriate transformed matrix
			stopifnot(all(is.na(X_mod[, inds_j])))

			X_mod[, inds_j] <- X_int[, inds_j] %*%
				.gen_inv_treat_block(
					num_treats,
					first_inds,
					G,
					fusion_structure,
					d_inv_treat = d_inv_treat
				)

			stopifnot(all(!is.na(X_mod[, inds_j])))
		}
	}

	stopifnot(all(!is.na(X_mod)))
	stopifnot(ncol(X_mod) == p)
	stopifnot(nrow(X_mod) == N * T)

	return(X_mod)
}


#' Back-transform Bridge-Regression Coefficients
#' \( \widehat{\boldsymbol\beta}
#'   = \boldsymbol D_N^{-1}\,\widehat{\boldsymbol\theta}\)
#'
#' @description
#' After fitting a bridge (Lasso, Elastic-Net, \eqn{\ell_q}) regression on the
#' transformed design matrix
#' \eqn{\widetilde{\boldsymbol Z}\,\boldsymbol D_N^{-1}}
#' (see `transformXintImproved()`), the solver returns parameter estimates
#' \eqn{\widehat{\boldsymbol\theta}}.
#' This helper multiplies that vector by the block-diagonal matrix
#' \eqn{\boldsymbol D_N^{-1}} from the paper,
#' thereby recovering the original-scale coefficients
#' \eqn{\widehat{\boldsymbol\beta}} for the FETWFE model.
#'
#' @details
#' The matrix
#' \deqn{
#' \boldsymbol D_N^{-1}
#'  = \operatorname{diag}\Bigl(
#'        (D^{(1)}(G))^{-1},\;
#'        (D^{(1)}(T\!-\!1))^{-1},\;
#'        I_d,\;
#'        (D^{(1)}(G))^{-1}\!\otimes\!I_d,\;
#'        (D^{(1)}(T\!-\!1))^{-1}\!\otimes\!I_d,\;
#'        (D^{(2)}(\mathcal G))^{-1},\;
#'        (D^{(2)}(\mathcal G))^{-1}\!\otimes\!I_d
#'     \Bigr)
#' }
#' corresponds to seven consecutive blocks in the column ordering used by
#' `transformXintImproved()`:
#' \enumerate{
#'   \item Cohort fixed-effect coefficients (length \eqn{G})
#'   \item Time fixed-effects (length \eqn{T-1})
#'   \item Main covariates (length \eqn{d})
#'   \item Covariate x cohort interactions (\eqn{dG})
#'   \item Covariate x time interactions \eqn{d(T-1)}
#'   \item Base treatment effects (\eqn{\mathfrak W =} `num_treats`)
#'   \item Covariate x treatment interactions (\eqn{d\mathfrak W})
#' }
#'
#' For each block the function premultiplies the slice of
#' \code{beta_hat_mod} with the *same* inverse-fusion matrix that was used as a
#' post-multiplier when building the transformed design matrix:
#'
#' | block | helper used | mathematical symbol |
#' |-------|-------------|---------------------|
#' | cohort FE | `genBackwardsInvFusionTransformMat(G)` | \((D^{(1)}(G))^{-1}\) |
#' | time FE | `genBackwardsInvFusionTransformMat(T-1)` | \((D^{(1)}(T-1))^{-1}\) |
#' | covariate blocks | identity copy | \(I_d\) |
#' | covariate x cohort | same helper, one copy per feature (interleaved, level-major) | \((D^{(1)}(G))^{-1}\otimes I_d\) |
#' | covariate x time | idem with \(T-1\) | \((D^{(1)}(T-1))^{-1}\otimes I_d\) |
#' | base treatment | `genInvTwoWayFusionTransformMat()` | \((D^{(2)}(\mathcal G))^{-1}\) |
#' | covariate x treatment | same two-way helper per feature (interleaved, level-major) | \((D^{(2)}(\mathcal G))^{-1}\otimes I_d\) |
#'
#' @param beta_hat_mod Numeric vector (length \eqn{p}).
#'   Estimated coefficients returned by a penalised fit on the transformed
#'   design matrix.
#' @param T Integer. Total time periods.
#' @param G Integer. Number of treated cohorts.
#' @param p Integer. Total number of columns in the **original** design matrix
#'   \eqn{p = p_N}.
#' @param d Integer. Number of time-invariant covariates.
#' @param num_treats Integer. Count of base treatment-effect dummies
#'   \eqn{\mathfrak W}.
#' @param first_inds Optional integer vector of length \eqn{G}.
#'   Starting indices (1-based, inside the treatment-dummy block) of the first
#'   effect for each cohort.  If `NA` they are reconstructed with
#'   **`getFirstInds()`**.
#' @param fusion_structure Character; one of `"cohort"` (default) or
#'   `"event_study"`. Selects which built-in inverse treatment-effect fusion
#'   block is used (ignored when `d_inv_treat` is supplied).
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   user-supplied already-inverted treatment-effect fusion block
#'   (`solve(fusion_matrix)`, #236). When non-`NULL` it overrides
#'   `fusion_structure` for the treatment block. Default `NULL`.
#'
#' @return Numeric vector \code{beta_hat} (length \eqn{p}) equal to
#'   \eqn{\boldsymbol D_N^{-1}\,\widehat{\boldsymbol\theta}}.
#'
#' @references
#' Wooldridge, J. M. (2021). Two-way fixed effects, the two-way mundlak
#' regression, and difference-in-differences estimators.
#' \emph{Available at SSRN 3906345}.
#' \doi{10.2139/ssrn.3906345}.
#' Tibshirani & Taylor (2011), "The Solution Path of the Generalized
#' Lasso".
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @seealso
#' * `transformXintImproved()` - forward transformation of the design matrix.
#' * `genBackwardsInvFusionTransformMat()`,
#'   `genInvTwoWayFusionTransformMat()` - block-wise inverse matrices.
#'
#' @examples
#' ## toy example: one covariate, two treated cohorts, T = 5
#' G <- 2; T <- 5; d <- 1
#' num_treats <- getNumTreats(G, T)
#' p <- G + (T-1) + d + d*G + d*(T-1) + num_treats + d*num_treats
#'
#' ## pretend we ran a penalised regression in the transformed space
#' beta_hat_mod <- rnorm(p)
#'
#' ## back-transform:
#' beta_hat <- untransformCoefImproved(
#'   beta_hat_mod, T=T, G=G, p=p, d=d,
#'   num_treats=num_treats
#' )
#' length(beta_hat)  # == p
#' @keywords internal
#' @noRd
untransformCoefImproved <- function(
	beta_hat_mod,
	T,
	G,
	p,
	d,
	num_treats,
	first_inds = NA,
	fusion_structure = "cohort",
	d_inv_treat = NULL
) {
	stopifnot(length(beta_hat_mod) == p)
	beta_hat <- rep(as.numeric(NA), p)

	if (any(is.na(first_inds))) {
		first_inds <- getFirstInds(G = G, T = T)
	}

	# First handle G cohort fixed effects effects
	beta_hat[1:G] <- genBackwardsInvFusionTransformMat(G) %*% beta_hat_mod[1:G]

	stopifnot(all(!is.na(beta_hat[1:G])))
	stopifnot(all(is.na(beta_hat[(G + 1):p])))

	# Next, T - 1 time fixed effects
	beta_hat[(G + 1):(G + T - 1)] <- genBackwardsInvFusionTransformMat(
		T - 1
	) %*%
		beta_hat_mod[(G + 1):(G + T - 1)]

	stopifnot(all(!is.na(beta_hat[1:(G + T - 1)])))
	stopifnot(all(is.na(beta_hat[(G + T):p])))

	# Coefficients for X (if any)
	if (d > 0) {
		beta_hat[(G + T):(G + T - 1 + d)] <- beta_hat_mod[
			(G + T):(G + T - 1 + d)
		]

		stopifnot(all(!is.na(beta_hat[1:(G + T - 1 + d)])))
		stopifnot(all(is.na(beta_hat[(G + T + d):p])))

		# Next, coefficients for cohort effects interacted with X. For each individual
		# feature, this will be handled using untransformVecFusion.
		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and cohort
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on G times
			feat_1 <- G + T - 1 + d + j
			feat_R <- G + T - 1 + d + (G - 1) * d + j
			feat_inds_j <- seq(feat_1, feat_R, by = d)
			stopifnot(length(feat_inds_j) == G)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(G) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(!is.na(beta_hat[1:(G + T - 1 + d + G * d)])))
		stopifnot(all(is.na(beta_hat[(G + T - 1 + d + G * d + 1):p])))

		# Similar for time effects interacted with X

		for (j in 1:d) {
			# Get indices corresponding to interactions between feature j and time
			# fixed effects--these are the first feature, the (1 + d)th feature, and
			# so on T - 1 times
			feat_1 <- G + T - 1 + d + G * d + j
			feat_T_minus_1 <- G + T - 1 + d + G * d + (T - 2) * d + j
			feat_inds_j <- seq(feat_1, feat_T_minus_1, by = d)
			stopifnot(length(feat_inds_j) == T - 1)

			stopifnot(all(is.na(beta_hat[feat_inds_j])))

			beta_hat[feat_inds_j] <- genBackwardsInvFusionTransformMat(
				T - 1
			) %*%
				beta_hat_mod[feat_inds_j]
			stopifnot(all(!is.na(beta_hat[feat_inds_j])))
		}
		stopifnot(all(
			!is.na(beta_hat[1:(G + T - 1 + d + G * d + (T - 1) * d)])
		))
		stopifnot(all(is.na(beta_hat[
			(G + T - 1 + d + G * d + (T - 1) * d + 1):p
		])))
	}

	# Now base treatment effects.

	feat_inds <- (G + T - 1 + d + G * d + (T - 1) * d + 1):(G +
		T -
		1 +
		d +
		G * d +
		(T - 1) * d +
		num_treats)

	stopifnot(all(is.na(beta_hat[feat_inds])))

	beta_hat[feat_inds] <- .gen_inv_treat_block(
		num_treats,
		first_inds,
		G,
		fusion_structure,
		d_inv_treat = d_inv_treat
	) %*%
		beta_hat_mod[feat_inds]

	if (d > 0) {
		stopifnot(all(
			!is.na(beta_hat[
				1:(G + T - 1 + d + G * d + (T - 1) * d + num_treats)
			])
		))
		stopifnot(all(is.na(beta_hat[
			(G + T - 1 + d + G * d + (T - 1) * d + num_treats + 1):p
		])))
		# Lastly, interactions between each treatment effect and each feature.
		# Feature-wise, we can do this with untransformTwoWayFusionCoefs, in the same
		# way that we did for previous interactions with X.
		for (j in 1:d) {
			# Recall that we have arranged the last d*num_Feats features in X_int
			# as follows: the first d are the first column of treat_mat_long interacted
			# with all of the columns of X, and so on. So, the columns that interact
			# the jth feature with all of the treatment effects are columns j, j + 1*d,
			# j + 2*d, ..., j + (num_treats - 1)*d.
			inds_j <- seq(j, j + (num_treats - 1) * d, by = d)
			stopifnot(length(inds_j) == num_treats)
			inds_j <- inds_j + G + T - 1 + d + G * d + (T - 1) * d + num_treats

			# Now ready to untransform the estimated coefficients
			stopifnot(all(is.na(beta_hat[inds_j])))

			beta_hat[inds_j] <- .gen_inv_treat_block(
				num_treats,
				first_inds,
				G,
				fusion_structure,
				d_inv_treat = d_inv_treat
			) %*%
				beta_hat_mod[inds_j]

			stopifnot(all(!is.na(beta_hat[inds_j])))
		}
	}

	stopifnot(all(!is.na(beta_hat)))

	return(beta_hat)
}


# genBackwardsFusionTransformMat
#' @title Generate Backward Fusion Transformation Matrix
#' @description Creates a square transformation matrix `D` of size `n_vars` x
#'   `n_vars`. When pre-multiplied by a coefficient vector `beta`, `D %*% beta`
#'   yields a transformed vector `theta` where `theta_i = beta_i - beta_{i+1}`
#'   for `i < n_vars`, and `theta_{n_vars} = beta_{n_vars}`. This is used to
#'   penalize coefficients towards the *next* coefficient in sequence, and the
#'   last coefficient directly.
#' @param n_vars Integer; the number of variables (coefficients) in the block
#'   to be transformed. This will be the dimension of the output matrix.
#' @return A numeric matrix of dimension `n_vars` x `n_vars`.
#' @details The resulting matrix `D` has 1s on the main diagonal. For each row
#'   `i` (from 1 to `n_vars - 1`), it has a -1 at column `i+1`. All other
#'   elements are 0.
#' @examples
#'   genBackwardsFusionTransformMat(3)
#'   # Output:
#'   #      [,1] [,2] [,3]
#'   # [1,]    1   -1    0
#'   # [2,]    0    1   -1
#'   # [3,]    0    0    1
#' @keywords internal
#' @noRd
genBackwardsFusionTransformMat <- function(n_vars) {
	# Generates D matrix in relation theta = D beta, where D beta is what
	# we want to penalize (for a single set of coefficients where we want to
	# penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D <- matrix(0, n_vars, n_vars)

	for (i in 1:n_vars) {
		for (j in 1:n_vars) {
			if (i == j) {
				D[i, j] <- 1
			}
			if (j == i + 1) {
				D[i, j] <- -1
			}
		}
	}

	return(D)
}


#' Inverse "Backwards-Difference" Transformation Matrix  \( \bigl(D^{(1)}(t)\bigr)^{-1} \)
#'
#' @description
#' Generates the \eqn{t\times t} upper-triangular matrix of 1s whose inverse is
#' the first-difference operator
#' \eqn{D^{(1)}(t)}
#' defined in Eq. (14) of the paper:
#' \deqn{
#'   D^{(1)}(t) =
#'   \begin{bmatrix}
#'     1 & -1 &        &        & 0\\
#'       &  1 & -1     &        &  \\
#'       &    & \ddots & \ddots &  \\
#'       &    &        & 1 & -1 \\
#'     0 &    &        &   &  1
#'   \end{bmatrix}.
#' }
#' Multiplying a block of coefficients by this inverse converts finite
#' differences back to cumulative sums, which is what the bridge-penalty
#' expects.
#'
#' @details
#' The matrix has ones on and strictly above the main diagonal
#' (\eqn{D^{-1}_{ij}=1_{\{i\le j\}}}).
#' It is its own Cholesky factor, \eqn{D^{-1} = U = U^\top}.
#'
#' After constructing `D_inv`, the routine calls
#' `genBackwardsFusionTransformMat(n_vars)`, forms both
#' `D %*% D_inv` and `D_inv %*% D`, and checks that their
#' maximum element-wise deviation from the identity matrix is below
#' `tol = 1e-12`.  Failing the check raises an error, protecting downstream
#' computations from a silent algebraic bug.
#'
#' @param n_vars Integer. Dimension \eqn{t}.
#'
#' @return A \eqn{n\_vars \times n\_vars} numeric matrix equal to
#'   \(\bigl(D^{(1)}(n\_vars)\bigr)^{-1}\).
#' @keywords internal
#' @noRd
genBackwardsInvFusionTransformMat <- function(n_vars) {
	# Generates inverse of D matrix in relation theta = D beta, where D beta is
	# what we want to penalize (for a single set of coefficients where we want
	# to penalize last coefficient directly and penalize remaining coefficints
	# towards the next coefficient)
	D_inv <- matrix(0, n_vars, n_vars)

	diag(D_inv) <- 1

	D_inv[upper.tri(D_inv)] <- 1

	stopifnot(nrow(D_inv) == n_vars)
	stopifnot(ncol(D_inv) == n_vars)

	## -- self-test: D_inv %*% D  ==  I
	D <- genBackwardsFusionTransformMat(n_vars)
	tol <- 1e-12
	if (!isTRUE(all.equal(D_inv %*% D, diag(n_vars), tolerance = tol))) {
		stop(
			"genBackwardsInvFusionTransformMat(): self-test failed - ",
			"result is not the matrix inverse of D."
		)
	}

	return(D_inv)
}


#' Inverse Two-Way-Fusion Transformation Matrix  \( \bigl(D^{(2)}(\mathcal G)\bigr)^{-1} \)
#'
#' @description
#' Constructs the **inverse** of the block-lower-triangular matrix
#' \eqn{D^{(2)}(\mathcal G)} that appears in Lemma \eqn{3} of the paper.
#' Each treated cohort \(g\) has a "first" post-treatment coefficient
#' \(\tau_{g,0}\) and a string of subsequent coefficients
#' \(\tau_{g,1},\dots,\tau_{g,T-g}\).
#' The fusion penalty (1) pulls every \(\tau_{g,k}\;(k>0)\) toward its
#' predecessor \(\tau_{g,k-1}\) *within* the same cohort **and**
#' (2) pulls the first effect of cohort \(g\) toward the first effect of cohort
#' \(g-1\).
#' Multiplying the raw design sub-matrix by the output of
#' `genInvTwoWayFusionTransformMat()` therefore changes coordinates from
#' \eqn{\boldsymbol\beta} to
#' \eqn{\boldsymbol\theta = D^{(2)}(\mathcal G)\,\boldsymbol\beta},
#' so that an \eqn{\ell_q} penalty on \eqn{\theta} is exactly the desired
#' two-level fusion penalty on \eqn{\beta}.
#'
#' @details
#' * The returned matrix is **lower-triangular with 1s** on and below the main
#'   diagonal and zeros elsewhere, except for extra 1s that implement the
#'   cross-cohort link on the first effect of each cohort (see paper, Eq. (18)).
#' * Its inverse contains only \{-1,0,1\} and recreates Eq. (17) of the paper.
#' * Determinant is 1, so the transformation is volume-preserving.
#'
#' @param n_vars Integer. Total number of base treatment-effect coefficients
#'   ( \eqn{\mathfrak W}  in the paper).
#' @param first_inds Integer vector of length \eqn{G}.
#'   `first_inds[g]` is the **1-based** column index of \(\tau_{g,0}\)
#'   inside the block of the \eqn{n\_vars} treatment columns.
#' @param G Integer. Number of treated cohorts.
#'
#' @return A numeric matrix of size \eqn{n\_vars \times n\_vars} that is
#'   exactly \(\bigl(D^{(2)}(\mathcal G)\bigr)^{-1}\).
#'
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#'
#' @examples
#' G  <- 3;  T <- 6
#' num_treats <- getNumTreats(G, T)
#' first <- getFirstInds(G, T)
#' Dinv <- genInvTwoWayFusionTransformMat(num_treats, first, G)
#' # verify Dinv %*% solve(Dinv) = I
#' all.equal(Dinv %*% solve(Dinv), diag(num_treats))
#' @keywords internal
#' @noRd
genInvTwoWayFusionTransformMat <- function(n_vars, first_inds, G) {
	stopifnot(length(n_vars) == 1, n_vars >= 0)
	stopifnot(is.numeric(G), length(G) == 1, G >= 0)
	if (G > 0) {
		stopifnot(is.numeric(first_inds), length(first_inds) == G)
		# Add checks for consistency of first_inds if n_vars > 0
		if (n_vars > 0) {
			stopifnot(first_inds[1] == 1)
			if (G > 1) {
				stopifnot(all(diff(first_inds) > 0)) # Must be strictly increasing
			}
			# Check that the last effect of the last cohort aligns with n_vars
			# M_R = n_vars - first_inds[G] + 1. This M_R must be > 0 if first_inds[G] <= n_vars.
			stopifnot(first_inds[G] <= n_vars)
		} else {
			# n_vars == 0
			stopifnot(G == 0) # If n_vars is 0, G must be 0. first_inds should be empty.
		}
	} else {
		# G == 0
		stopifnot(length(first_inds) == 0, n_vars == 0)
	}

	D_inv <- matrix(0, nrow = n_vars, ncol = n_vars)

	if (n_vars == 0) {
		# Handles G=0 correctly
		return(D_inv)
	}

	# Part 1: Set up the column structure for \tilde{U} blocks and
	# the first column of each diagonal block.
	# For each cohort `c` (from 1 to G), the column in D_inv corresponding to its
	# first treatment effect (i.e., absolute column index `first_inds[c]`)
	# should have 1s from its own row (`first_inds[c]`) down to `n_vars`.
	if (G > 0) {
		for (cohort_idx in 1:G) {
			col_to_fill <- first_inds[cohort_idx]
			D_inv[col_to_fill:n_vars, col_to_fill] <- 1
		}
	}

	# Part 2: Form the correct diagonal blocks.
	# Each diagonal block is (D^{(1)}(M_c)^{-1})^T, which is a fully
	# lower triangular matrix of 1s (all elements on and below diagonal are 1).
	# This will correctly overwrite the 1s placed in Part 1 for these diagonal parts.
	if (G > 0) {
		for (cohort_idx in 1:G) {
			block_start_idx <- first_inds[cohort_idx]

			if (cohort_idx < G) {
				block_end_idx <- first_inds[cohort_idx + 1] - 1
			} else {
				# This is the last cohort
				block_end_idx <- n_vars
			}

			# Ensure block indices are valid (e.g. cohort has at least one effect)
			if (block_start_idx > block_end_idx) {
				# This case implies cohort_idx has zero effects.
				# This should ideally not happen if first_inds and n_vars are consistent
				# with cohorts having at least one effect.
				# If it can happen, 'continue' or 'next' might be appropriate.
				# For now, assume M_x >= 1 for all cohorts included in G.
				next
			}

			# Fill this diagonal block D_inv[block_start_idx:block_end_idx, block_start_idx:block_end_idx]
			for (row_abs in block_start_idx:block_end_idx) {
				# Absolute row index in D_inv
				# For the current row_abs within its block, set columns from
				# the start of the block (block_start_idx) up to the current row_abs to 1.
				D_inv[row_abs, block_start_idx:row_abs] <- 1
			}
		}
	}

	# The original code had a stop for G < 2.
	# This revised logic handles G=0 and G=1 correctly without a special stop.
	# If G=1, Part 1 sets D_inv[1:n_vars, 1] <- 1.
	# Part 2 (with cohort_idx=1, block_start_idx=1, block_end_idx=n_vars) then correctly
	# forms the full (D^{(1)}(n_vars)^{-1})^T matrix.

	return(D_inv)
}


#' Build the inverse event-study fusion transform matrix \((D^{(2)}_{ES})^{-1}\)
#'
#' Constructs the inverse of the event-study treatment-effect differences matrix
#' \eqn{\boldsymbol{D}^{(2)}_{\mathrm{ES}}(\mathcal G)} of Faletto (2025) Lemma
#' \code{event.study.sing.val.lem} (closed form \code{es.2.inv.exp}). The
#' event-study penalty fuses treatment effects at the same time since treatment
#' (event time \eqn{e = t - g}) across cohorts, in place of the default
#' within-cohort / between-cohort two-way fusion built by
#' \code{genInvTwoWayFusionTransformMat()}.
#'
#' With treatment effects ordered by cohort --- cohort \eqn{g_k} occupies
#' \code{first_inds[k]:(first_inds[k + 1] - 1)}, of size \eqn{n_k = T - g_k + 1}
#' --- the inverse is block lower-triangular: its first column-block is the
#' first \eqn{n_k} rows of the all-ones lower-triangular matrix (the "spine" of
#' the earliest cohort \eqn{g_1}); its strictly-lower off-diagonal blocks are
#' same-event-time selectors \eqn{S(g_j, g_k) = [\boldsymbol I_{n_k}\ \boldsymbol
#' 0]}; and its diagonal blocks are \eqn{\boldsymbol I_{n_k}} for \eqn{k \ge 2}
#' (or the all-ones lower-triangular spine for \eqn{k = 1}).
#'
#' @param n_vars Integer; the number of treatment effects \eqn{\mathfrak W}.
#' @param first_inds Integer vector of length \code{G}; the index, within the
#'   treatment-effect block, of each cohort's first treatment effect.
#' @param G Integer; the number of treated cohorts.
#' @return The \code{n_vars x n_vars} matrix \eqn{(D^{(2)}_{ES})^{-1}}.
#' @keywords internal
#' @noRd
genInvEventStudyFusionTransformMat <- function(n_vars, first_inds, G) {
	stopifnot(length(n_vars) == 1, n_vars >= 0)
	stopifnot(is.numeric(G), length(G) == 1, G >= 0)
	if (G > 0) {
		stopifnot(is.numeric(first_inds), length(first_inds) == G)
		if (n_vars > 0) {
			stopifnot(first_inds[1] == 1)
			if (G > 1) {
				stopifnot(all(diff(first_inds) > 0))
			}
			stopifnot(first_inds[G] <= n_vars)
		} else {
			stopifnot(G == 0)
		}
	} else {
		stopifnot(length(first_inds) == 0, n_vars == 0)
	}

	D_inv <- matrix(0, nrow = n_vars, ncol = n_vars)
	if (n_vars == 0) {
		return(D_inv)
	}

	# Cohort block boundaries and sizes n_k = T - g_k + 1.
	block_starts <- first_inds
	block_ends <- c(first_inds[-1] - 1L, n_vars)
	n_k <- block_ends - block_starts + 1L
	n1 <- n_k[1]
	cols1 <- block_starts[1]:block_ends[1]

	# (tilde D^{(1)}(1)^{-1})^T is the n1 x n1 all-ones lower-triangular matrix.
	LT1 <- matrix(0, n1, n1)
	LT1[lower.tri(LT1, diag = TRUE)] <- 1

	for (k in 1:G) {
		rows_k <- block_starts[k]:block_ends[k]
		# Block (k, 1): first n_k rows of LT1, at the earliest cohort's columns.
		# For k = 1 this is the full spine block (the only directly penalized
		# base term plus its adjacent-event-time differences).
		D_inv[rows_k, cols1] <- LT1[seq_len(n_k[k]), , drop = FALSE]
		if (k >= 2) {
			# Strictly-lower off-diagonal same-event-time selectors
			# S(g_j, g_k) = [I_{n_k} | 0].
			if (k >= 3) {
				for (j in 2:(k - 1)) {
					cols_j <- block_starts[j]:block_ends[j]
					D_inv[rows_k, cols_j[seq_len(n_k[k])]] <- diag(n_k[k])
				}
			}
			# Diagonal block I_{n_k}.
			D_inv[rows_k, rows_k] <- diag(n_k[k])
		}
	}

	D_inv
}


# .gen_inv_treat_block
#' @title Dispatch the inverse treatment-effect fusion block by structure
#' @description Returns the inverse of the treatment-effect fusion difference
#'   matrix for the requested `fusion_structure`: the default cohort two-way
#'   fusion (`genInvTwoWayFusionTransformMat()`) or the event-study fusion
#'   (`genInvEventStudyFusionTransformMat()`, #40). A single dispatch point so
#'   every consumer (transform / untransform / variance / accessors) stays in
#'   sync, and the `"cohort"` (default) branch is the exact pre-existing call.
#'
#'   When `d_inv_treat` is supplied (a user `fusion_matrix`'s already-inverted
#'   `num_treats x num_treats` block; #236), it is returned verbatim, overriding
#'   `fusion_structure`. This is the single override choke point: every consumer
#'   reaches it, so threading `d_inv_treat` here makes the user's custom block
#'   reach the transform / untransform / variance / ridge / accessor paths at
#'   once. With `d_inv_treat = NULL` (the default) the call is byte-identical to
#'   the prior `fusion_structure` dispatch.
#' @param num_treats,first_inds,G As in `genInvTwoWayFusionTransformMat()`.
#' @param fusion_structure One of `"cohort"` or `"event_study"`.
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   already-inverted user-supplied treatment-effect fusion block
#'   (`solve(fusion_matrix)`). If non-`NULL`, overrides `fusion_structure` and is
#'   returned as-is. Default `NULL`.
#' @return The `num_treats x num_treats` inverse treatment-effect fusion block.
#' @keywords internal
#' @noRd
.gen_inv_treat_block <- function(
	num_treats,
	first_inds,
	G,
	fusion_structure = "cohort",
	d_inv_treat = NULL
) {
	if (!is.null(d_inv_treat)) {
		stopifnot(identical(
			dim(d_inv_treat),
			c(as.integer(num_treats), as.integer(num_treats))
		))
		return(d_inv_treat)
	}
	if (identical(fusion_structure, "event_study")) {
		genInvEventStudyFusionTransformMat(num_treats, first_inds, G)
	} else {
		genInvTwoWayFusionTransformMat(num_treats, first_inds, G)
	}
}


# .validate_fusion_matrix
#' @title Validate a user-supplied `fusion_matrix` and return its inverse block
#' @description Validates the optional user `fusion_matrix` argument of
#'   `fetwfe()` (#236) and returns the already-inverted treatment-effect block
#'   `solve(fusion_matrix)` (a `num_treats x num_treats` matrix) that the
#'   estimator threads through the `.gen_inv_treat_block()` choke point. When
#'   `fusion_matrix` is `NULL` (the default) it returns `NULL` unchanged, so the
#'   built-in `fusion_structure` dispatch is used and the fit is byte-identical
#'   to omitting the argument.
#'
#'   Checks (in order), all referencing the cohort-major `(g, t)` row order of
#'   `getFirstInds()` / `getTreatInds()`:
#'   \enumerate{
#'     \item type/shape: a finite numeric matrix with
#'       `nrow == ncol == num_treats` (else `stop()`, naming the expected
#'       `num_treats` derived from `G` and `T`);
#'     \item invertibility: `solve(fusion_matrix)` must succeed (else `stop()`);
#'     \item conditioning: a `warning()` (not error) when `D_N` is numerically
#'       near-singular (reciprocal condition number `< 1e-10`) -- it stays
#'       invertible and yields a valid point estimator, but
#'       `solve(fusion_matrix)` may be unreliable. Under the paper's
#'       fixed-dimension scoping any finite invertible `D_N` inherits the
#'       guarantees (Assumption (D)), so the singular-value magnitudes
#'       themselves are not policed.
#'   }
#'   When a non-default `fusion_structure` was also supplied, a `message()`
#'   (gated on `verbose`) notes that `fusion_matrix` overrides it for the
#'   treatment block.
#' @param fusion_matrix The user-supplied forward difference matrix `D_N`, or
#'   `NULL`.
#' @param num_treats Integer; the number of base treatment effects (the required
#'   matrix dimension).
#' @param G Integer; number of treated cohorts (for the error message).
#' @param T Integer; number of time periods (for the error message).
#' @param fusion_structure_supplied Logical; `TRUE` if the user explicitly
#'   passed `fusion_structure` (drives the override-notice message). Default
#'   `FALSE`.
#' @param verbose Logical; gates the override-notice `message()`. Default
#'   `FALSE`.
#' @return `NULL` when `fusion_matrix` is `NULL`; otherwise the
#'   `num_treats x num_treats` numeric matrix `solve(fusion_matrix)`.
#' @references
#' Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
#' Difference-in-Differences with Staggered Adoptions.
#' \emph{arXiv preprint arXiv:2312.05985}.
#' \url{https://arxiv.org/abs/2312.05985}.
#' @keywords internal
#' @noRd
.validate_fusion_matrix <- function(
	fusion_matrix,
	num_treats,
	G,
	T,
	fusion_structure_supplied = FALSE,
	verbose = FALSE
) {
	if (is.null(fusion_matrix)) {
		return(NULL)
	}

	# --- 1. type / shape -------------------------------------------------
	if (
		!is.matrix(fusion_matrix) ||
			!is.numeric(fusion_matrix) ||
			!all(is.finite(fusion_matrix))
	) {
		stop(
			"fusion_matrix must be a finite numeric matrix (it is the ",
			"num_treats x num_treats forward-differences matrix D_N). Got an ",
			"object of class ",
			paste(class(fusion_matrix), collapse = "/"),
			if (is.matrix(fusion_matrix) && !all(is.finite(fusion_matrix))) {
				" containing non-finite (NA/NaN/Inf) entries."
			} else {
				"."
			},
			call. = FALSE
		)
	}
	if (
		nrow(fusion_matrix) != num_treats || ncol(fusion_matrix) != num_treats
	) {
		stop(
			"fusion_matrix must be a num_treats x num_treats numeric matrix; ",
			"got ",
			nrow(fusion_matrix),
			" x ",
			ncol(fusion_matrix),
			", expected num_treats = ",
			num_treats,
			" (from G = ",
			G,
			", T = ",
			T,
			"). The rows/columns follow the cohort-major (g, t) order of ",
			"getFirstInds() / getTreatInds().",
			call. = FALSE
		)
	}

	# --- 2. invertibility ------------------------------------------------
	d_inv_treat <- tryCatch(
		solve(fusion_matrix),
		error = function(e) {
			stop(
				"fusion_matrix is singular (or numerically non-invertible) and ",
				"cannot be inverted: ",
				conditionMessage(e),
				". It must be an invertible forward-differences matrix D_N.",
				call. = FALSE
			)
		}
	)

	# --- 3. conditioning heads-up (warning, not error) ------------------
	# Under the paper's fixed-dimension scoping, any finite invertible D_N
	# satisfies Assumption (D) and inherits the inferential guarantees -- its
	# singular-value magnitudes affect only constant factors -- so we do not
	# police them. We do warn when D_N is numerically near-singular: solve()
	# above succeeded, but its inverse may be unreliable. The threshold sits well
	# above R's hard "computationally singular" error (~machine epsilon) and far
	# below any legitimate D_N (the built-in penalties have reciprocal condition
	# number ~1 / (T * sqrt(2 * T)), orders of magnitude larger). rcond() is an
	# O(num_treats^2) LAPACK estimate -- much cheaper than the former full svd().
	rcond_floor <- 1e-10
	rc <- rcond(fusion_matrix)
	if (is.finite(rc) && rc < rcond_floor) {
		warning(
			"fusion_matrix is numerically near-singular (reciprocal condition ",
			"number = ",
			format(rc, digits = 6),
			" < ",
			format(rcond_floor, digits = 6),
			"). It is still invertible and yields a valid point estimator, but ",
			"solve(fusion_matrix) may be numerically unreliable; consider a ",
			"better-conditioned D_N.",
			call. = FALSE
		)
	}

	# --- 4. override notice ----------------------------------------------
	if (isTRUE(fusion_structure_supplied) && isTRUE(verbose)) {
		message(
			"fusion_matrix was supplied, so it overrides fusion_structure for ",
			"the treatment-effect block."
		)
	}

	d_inv_treat
}


#' Build the **entire** inverse-fusion matrix \(D_N^{-1}\)
#'
#' Constructs the \eqn{p\times p} block-diagonal matrix that sends the sparse
#' coefficient vector \(\theta\) used in the bridge loss back to the original
#' coefficient scale \(\beta\), and is also used to multiply \(Z\) on the right
#' to transform the features into a transformed feature space:
#' \deqn{\beta \;=\;D_N^{-1}\,\theta.}
#' The block layout follows the factorisation proved in the
#' paper -- repeated here using the helper generators already available in the
#' package.
#'
#' \preformatted{
#'       D_N^{-1} = diag(
#'         (D^{(1)}(G))^{-1},                          # 1. cohort FEs
#'         (D^{(1)}(T-1))^{-1},                        # 2. time  FEs
#'         I_d,                                        # 3. X main effects
#'         (D^{(1)}(G))^{-1} ⊗ I_d,                   # 4. cohort × X
#'         (D^{(1)}(T-1))^{-1} ⊗ I_d,                 # 5. time   × X
#'         (D^{(2)}(G))^{-1},                         # 6. base τ_{g,t}
#'         (D^{(2)}(G))^{-1} ⊗ I_d )                  # 7. τ_{g,t} × X
#' }
#'
#' @section Block dimensions:
#' \itemize{
#'   \item Cohort FEs: \(G\times G\)
#'   \item Time-period FEs: \((T-1)\times(T-1)\)
#'   \item Identitites: \(d\), \(dG\), \(d(T-1)\)
#'   \item Treatment blocks: \( \mathfrak W \times \mathfrak W\) with
#'         \(\mathfrak W = \texttt{num_treats}\)
#' }
#' All non-identity blocks contain only 0/1 entries, so the determinant of
#' the whole matrix is 1 (volume-preserving transform).
#'
#' @param first_inds Integer vector (length \code{G}).
#'   `first_inds[g]` is the 1-based column index of the first base
#'   treatment-effect parameter \(\tau_{g,0}\) for cohort \code{g}.
#' @param T Integer. Total number of time periods \(\ge 2\).
#' @param G Integer. Number of treated cohorts (\(\ge 1\)).
#'   The function stops if you accidentally pass \code{G = 0}.
#' @param d Integer. Number of time-invariant covariates (can be 0).
#' @param num_treats Integer.  Total number of base treatment-effect
#'   coefficients \(\mathfrak W = T G - G(G+1)/2\).
#' @param fusion_structure Character; one of `"cohort"` (default) or
#'   `"event_study"`. Selects which built-in inverse treatment-effect fusion
#'   block is embedded (ignored when `d_inv_treat` is supplied).
#' @param d_inv_treat Optional `num_treats x num_treats` numeric matrix; the
#'   user-supplied already-inverted treatment-effect fusion block
#'   (`solve(fusion_matrix)`, #236). When non-`NULL` it overrides
#'   `fusion_structure` for the treatment block embedded in the full
#'   block-diagonal. Default `NULL`.
#'
#' @return A dense base-R matrix of size
#'   \eqn{p \times p} with
#'   \eqn{p = G + (T-1) + d + dG + d(T-1) + \mathfrak W + d\mathfrak W}.
#'
#' @examples
#' G  <- 3; T <- 6; d <- 2
#' nt <- getNumTreats(G, T)
#' Dinv <- genFullInvFusionTransformMat(getFirstInds(G,T), T, G, d, nt)
#' dim(Dinv)   # should be p x p
#'
#' @keywords internal
#' @noRd
genFullInvFusionTransformMat <- function(
	first_inds,
	T,
	G,
	d,
	num_treats,
	fusion_structure = "cohort",
	d_inv_treat = NULL
) {
	##———— Safety checks ————————————————————————————————————————————
	stopifnot(is.numeric(G), length(G) == 1L, G >= 1L)
	stopifnot(is.numeric(T), length(T) == 1L, T >= 2L, G <= T - 1)
	stopifnot(is.numeric(d), length(d) == 1L, d >= 0L)
	stopifnot(length(first_inds) == G)

	##———— 1. Cohort fixed-effects block:  (D^{(1)}(G))^{-1} ————————
	block1 <- genBackwardsInvFusionTransformMat(G)

	##———— 2. Time fixed-effects block:    (D^{(1)}(T-1))^{-1} ———————
	block2 <- genBackwardsInvFusionTransformMat(T - 1)

	##———— 3. Covariate main effects:      I_d ————————————————
	block3 <- if (d > 0) diag(d) else NULL

	##———— 4. Cohort × X interactions:     (D^{(1)}(G))^{-1} ⊗ I_d ——
	block4 <- if (d > 0) {
		kronecker(genBackwardsInvFusionTransformMat(G), diag(d))
	} else {
		NULL
	}

	##———— 5. Time × X interactions:       (D^{(1)}(T-1))^{-1} ⊗ I_d —
	block5 <- if (d > 0) {
		kronecker(genBackwardsInvFusionTransformMat(T - 1), diag(d))
	} else {
		NULL
	}

	##———— 6. Base treatment effects:      (D^{(2)}(G))^{-1} ————————
	block6 <- .gen_inv_treat_block(
		num_treats = num_treats,
		first_inds = first_inds,
		G = G,
		fusion_structure = fusion_structure,
		d_inv_treat = d_inv_treat
	)

	##———— 7. Treatment × X interactions:  (D^{(2)}(G))^{-1} ⊗ I_d ——
	block7 <- if (d > 0) {
		kronecker(
			.gen_inv_treat_block(
				num_treats = num_treats,
				first_inds = first_inds,
				G = G,
				fusion_structure = fusion_structure,
				d_inv_treat = d_inv_treat
			),
			diag(d)
		)
	} else {
		NULL
	}

	## Gather present blocks in the same order as the theoretical expression
	blocks <- list(block1, block2, block3, block4, block5, block6, block7)
	blocks <- Filter(Negate(is.null), blocks) # drop NULLs for d = 0

	##———— Assemble block-diagonal matrix ————————————————
	## Matrix::bdiag() returns a sparse dgCMatrix.  We convert to base-R matrix
	## here because downstream (ridge-row augmentation) works with dense objects
	## via rbind().
	full_D_inv <- as.matrix(Matrix::bdiag(blocks))

	##———— Dimension cross-check ————————————————————————————
	p <- getP(G = G, T = T, d = d, num_treats = num_treats)
	stopifnot(nrow(full_D_inv) == p, ncol(full_D_inv) == p)

	return(full_D_inv)
}
