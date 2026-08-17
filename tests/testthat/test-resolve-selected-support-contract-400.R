library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #400 follow-up: a DIRECT contract test for .resolve_selected_support() (the
# selected-support class dispatch single-sourced in PR #414, Phase 2b). The
# cross-accessor guardrail (test-cross-accessor-scaffold-guardrail-400.R) only
# reaches these conventions INDIRECTLY, through full accessor output on ordinary
# fits. This pins them directly:
#   - etwfe/twfeCovs: sel_feat_inds is the scalar NA "all-features" sentinel, all
#     cells selected, d_inv_treat_sel = full diag(num_treats).
#   - fetwfe/betwfe: d_inv_treat_sel keeps drop = FALSE (a matrix even with a
#     single selected treatment cell) and is NULL when nothing is selected.
#
# #429 item 8: the k = 1 drop = FALSE case is now covered for BOTH branches. It
# used to be exercised for betwfe only, while this header claimed it for
# "fetwfe/betwfe" -- so removing drop = FALSE from the fetwfe branch of
# .resolve_selected_support() left this file green.
# ------------------------------------------------------------------------------

.resolve <- get(".resolve_selected_support", envir = asNamespace("fetwfe"))
.getFirstInds <- get("getFirstInds", envir = asNamespace("fetwfe"))

.args_for <- function(fit, beta_hat = fit$beta_hat) {
	treat_inds <- fit$treat_inds
	list(
		x = fit,
		treat_inds = treat_inds,
		num_treats = length(treat_inds),
		first_inds = .getFirstInds(G = fit$G, T = fit$T),
		G = fit$G,
		beta_hat = beta_hat,
		tes = beta_hat[treat_inds]
	)
}

.make_fits <- function() {
	cf <- genCoefs(G = 4, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 270)
	sim <- simulateData(
		cf,
		N = 160,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 270
	)
	base <- list(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		verbose = FALSE
	)
	list(
		fetwfe = do.call(fetwfe, base),
		betwfe = do.call(betwfe, base),
		etwfe = do.call(etwfe, base),
		twfeCovs = do.call(twfeCovs, base)
	)
}

fits <- .make_fits()

test_that("etwfe/twfeCovs use the NA all-features sentinel + full diag (#400)", {
	for (cls in c("etwfe", "twfeCovs")) {
		r <- do.call(.resolve, .args_for(fits[[cls]]))
		nt <- length(fits[[cls]]$treat_inds)
		expect_identical(r$sel_feat_inds, NA) # scalar NA sentinel
		expect_identical(r$sel_treat_inds_shifted, seq_len(nt))
		expect_true(is.matrix(r$d_inv_treat_sel))
		expect_identical(dim(r$d_inv_treat_sel), c(nt, nt))
		expect_equal(r$d_inv_treat_sel, diag(nt))
	}
})

test_that("fetwfe/betwfe d_inv_treat_sel keeps drop=FALSE on ordinary fits (#400)", {
	for (cls in c("fetwfe", "betwfe")) {
		r <- do.call(.resolve, .args_for(fits[[cls]]))
		k <- length(r$sel_treat_inds_shifted)
		expect_gte(k, 1L) # the fit selected something
		expect_true(is.matrix(r$d_inv_treat_sel))
		expect_identical(ncol(r$d_inv_treat_sel), k)
		expect_false(is.na(r$sel_feat_inds[1])) # NOT the sentinel for these classes
	}
})

test_that("betwfe d_inv_treat_sel stays a 1-col matrix with a single selection (#400 drop=FALSE)", {
	a <- .args_for(fits$betwfe)
	a$beta_hat[a$treat_inds] <- 0
	a$beta_hat[a$treat_inds[1]] <- 1 # exactly one selected treatment cell
	a$tes <- a$beta_hat[a$treat_inds]
	r <- do.call(.resolve, a)
	expect_length(r$sel_treat_inds_shifted, 1L)
	expect_true(is.matrix(r$d_inv_treat_sel)) # drop=FALSE: NOT dropped to a vector
	expect_identical(ncol(r$d_inv_treat_sel), 1L)
})

test_that("fetwfe d_inv_treat_sel stays a 1-col matrix with a single selection (#429 item 8)", {
	# The fetwfe branch has to be forced DIFFERENTLY from betwfe's above: betwfe
	# reads its selection from the `beta_hat` ARGUMENT, while fetwfe ignores
	# beta_hat entirely and reads x$internal$theta_hat. So the mutation goes on
	# the fit object, with the `1L +` offset the resolver's intercept-dropping
	# slice implies. The ordinary fit selects k = 6 here (a 10x6 subset is not
	# dropped), which is why the k >= 1 block above cannot see the mutation.
	#
	# This deliberately leaves `f` OFF-CONTRACT: its theta_hat no longer agrees
	# with its beta_hat, att_selected, catt_hats or catt_df, so
	# .validate_fetwfe() would reject it. That is safe only because
	# .resolve_selected_support() is a pure internal reader with no precondition
	# call. DO NOT extend this block to call any of the seven
	# .check_for_*-guarded accessors (eventStudy, augment, tidy, glance, plot,
	# coef, simultaneousCIs) on `f` -- each would fail on a contract unrelated
	# to drop = FALSE.
	f <- fits$fetwfe
	num_treats <- length(f$treat_inds)
	f$internal$theta_hat[1L + f$treat_inds] <- 0
	f$internal$theta_hat[1L + f$treat_inds[1]] <- 1 # exactly one selected cell
	r <- do.call(.resolve, .args_for(f))

	expect_length(r$sel_treat_inds_shifted, 1L)
	# The two load-bearing assertions: without drop = FALSE, d_inv_treat[, 1] is
	# a bare length-num_treats vector, so is.matrix() is FALSE and dim() is NULL.
	expect_true(is.matrix(r$d_inv_treat_sel))
	expect_identical(dim(r$d_inv_treat_sel), c(as.integer(num_treats), 1L))

	# Branch identity: the fetwfe branch must not be silently producing what the
	# betwfe branch would from the same selection. Measured here, the fetwfe
	# column is the all-ones vector (cohort 1's fused base column lifts every
	# CATT equally), which diag(num_treats)[, 1] is not. Do NOT "strengthen"
	# this into a pin on the column's exact contents -- that is
	# genInvTwoWayFusionTransformMat()'s contract and is tested elsewhere.
	expect_false(identical(r$d_inv_treat_sel[, 1], diag(num_treats)[, 1]))
})

test_that("fetwfe/betwfe d_inv_treat_sel is NULL on empty selection (#400)", {
	# betwfe: zero the beta-space treatment cells
	a <- .args_for(fits$betwfe)
	a$beta_hat[a$treat_inds] <- 0
	a$tes <- a$beta_hat[a$treat_inds]
	rb <- do.call(.resolve, a)
	expect_length(rb$sel_treat_inds_shifted, 0L)
	expect_null(rb$d_inv_treat_sel)
	# fetwfe: zero the theta-space treatment slopes (theta_hat[1 + treat_inds])
	f <- fits$fetwfe
	f$internal$theta_hat[1L + f$treat_inds] <- 0
	rf <- do.call(.resolve, .args_for(f))
	expect_length(rf$sel_treat_inds_shifted, 0L)
	expect_null(rf$d_inv_treat_sel)
})
