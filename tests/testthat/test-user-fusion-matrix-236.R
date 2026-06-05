library(testthat)
library(fetwfe)

# Issue #236: user-supplied `fusion_matrix` (a num_treats x num_treats forward
# differences matrix D_N) that overrides `fusion_structure` for the
# treatment-effect block. The estimator inverts it once
# (`d_inv_treat <- solve(fusion_matrix)`) and threads the inverse through the
# single `.gen_inv_treat_block()` choke point, storing the inverse on the fit so
# the access-time accessors (eventStudy() / simultaneousCIs()) reuse it.
#
# Headline guarantee: feeding the event-study built-in's forward D_N as a custom
# matrix reproduces the `fusion_structure = "event_study"` fit byte-identically.
#
# The fixture below is DISCRIMINATING (per the plan-review finding): its selected
# treatment support hits treatment columns where the cohort and event-study
# inverse blocks differ, so the accessor-time stored-slot read is genuinely
# exercised (a fixture selecting only coinciding columns would make the
# mutation-check vacuous). The test asserts that intersection at RUNTIME rather
# than hard-coding the support, so it stays valid if selection shifts.

# ------------------------------------------------------------------------------
# Shared fixture (pinned): genCoefs seed 1, density 0.9, eff_size 8, N 120,
# data seed 101, lambda_selection = "bic", q = 0.5.
# ------------------------------------------------------------------------------

.fm236_sim <- function() {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.9,
		eff_size = 8,
		seed = 1
	)
	set.seed(101)
	simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5)
}

.fm236_fit <- function(sim, ...) {
	suppressWarnings(fetwfeWithSimulatedData(
		sim,
		q = 0.5,
		lambda_selection = "bic",
		...
	))
}

# The event-study built-in's forward D_N for this fixture's dimensions.
.fm236_event_study_DN <- function(fit) {
	num_treats <- length(fit$treat_inds)
	first_inds <- getFirstInds(fit$G, fit$T)
	solve(fetwfe:::genInvEventStudyFusionTransformMat(
		num_treats,
		first_inds,
		fit$G
	))
}

# ------------------------------------------------------------------------------
# 1. NULL no-op: fusion_matrix = NULL is byte-identical to omitting it, under
#    both fusion_structure settings. (Choke-point else branch unchanged; every
#    threaded d_inv_treat = NULL default leaves existing calls byte-identical.)
# ------------------------------------------------------------------------------

test_that("fusion_matrix = NULL is byte-identical to omitting it (cohort + event_study)", {
	skip_on_cran()
	sim <- .fm236_sim()

	fit_co_omit <- .fm236_fit(sim, fusion_structure = "cohort")
	fit_co_null <- .fm236_fit(
		sim,
		fusion_structure = "cohort",
		fusion_matrix = NULL
	)
	expect_identical(fit_co_omit, fit_co_null)

	fit_es_omit <- .fm236_fit(sim, fusion_structure = "event_study")
	fit_es_null <- .fm236_fit(
		sim,
		fusion_structure = "event_study",
		fusion_matrix = NULL
	)
	expect_identical(fit_es_omit, fit_es_null)

	# The NULL-default still materializes the (NULL) slots, so the slot names
	# survive for doc-slot parity.
	expect_true("fusion_matrix" %in% names(fit_co_omit))
	expect_null(fit_co_omit$fusion_matrix)
	expect_true("d_inv_treat" %in% names(fit_co_omit$internal))
	expect_null(fit_co_omit$internal$d_inv_treat)
})

# ------------------------------------------------------------------------------
# 2. Headline cross-check: a custom D_N equal to the event-study built-in's
#    forward matrix, supplied with fusion_structure = "cohort", reproduces the
#    fusion_structure = "event_study" fit byte-identically -- including the
#    SE/CI columns of eventStudy() / simultaneousCIs(), which exercise the
#    accessor-time stored-slot read.
# ------------------------------------------------------------------------------

test_that("custom event-study D_N reproduces the event_study fit byte-identically", {
	skip_on_cran()
	sim <- .fm236_sim()

	fit_es <- .fm236_fit(sim, fusion_structure = "event_study")
	D_N <- .fm236_event_study_DN(fit_es)
	fit_custom <- .fm236_fit(
		sim,
		fusion_structure = "cohort",
		fusion_matrix = D_N
	)

	# Fit-time quantities.
	expect_identical(fit_custom$beta_hat, fit_es$beta_hat)
	expect_identical(fit_custom$att_hat, fit_es$att_hat)
	expect_identical(fit_custom$att_se, fit_es$att_se)
	expect_identical(fit_custom$catt_df, fit_es$catt_df)
	expect_identical(fit_custom$catt_ses, fit_es$catt_ses)

	# Confirm the fixture genuinely discriminates: the cohort fit differs, and
	# the selected support hits columns where cohort and event-study inverse
	# blocks differ.
	fit_co <- .fm236_fit(sim, fusion_structure = "cohort")
	expect_false(isTRUE(all.equal(fit_co$beta_hat, fit_es$beta_hat)))

	num_treats <- length(fit_es$treat_inds)
	first_inds <- getFirstInds(fit_es$G, fit_es$T)
	theta_slopes <- fit_es$internal$theta_hat[2:(fit_es$p + 1)]
	sel_treat <- which(theta_slopes[fit_es$treat_inds] != 0)
	d_co <- fetwfe:::.gen_inv_treat_block(
		num_treats,
		first_inds,
		fit_es$G,
		"cohort"
	)
	d_es <- fetwfe:::.gen_inv_treat_block(
		num_treats,
		first_inds,
		fit_es$G,
		"event_study"
	)
	diff_cols <- which(colSums(abs(d_co - d_es)) > 0)
	expect_gt(length(intersect(sel_treat, diff_cols)), 0)

	# Accessor-time SE/CI columns (these exercise x$internal$d_inv_treat).
	es_es <- eventStudy(fit_es)
	es_custom <- eventStudy(fit_custom)
	expect_identical(es_custom$se, es_es$se)
	expect_identical(es_custom$ci_low, es_es$ci_low)
	expect_identical(es_custom$ci_high, es_es$ci_high)
	expect_identical(es_custom$p_value, es_es$p_value)

	sc_es <- simultaneousCIs(fit_es)
	sc_custom <- simultaneousCIs(fit_custom)
	expect_identical(sc_custom$ci, sc_es$ci)
	expect_identical(sc_custom$adjusted_p_values, sc_es$adjusted_p_values)
})

# ------------------------------------------------------------------------------
# 2b. Mutation check for the accessor stored-slot read. If the eventStudy()
#     accessor ignored x$internal$d_inv_treat, the custom fit's event-study SEs
#     would fall back to the cohort built-in and DISAGREE with the event_study
#     fit. We prove the assertion in test 2 bites by simulating that mutation
#     here (recomputing the SE/CI the accessor would produce if it dropped the
#     stored block) and confirming it differs.
# ------------------------------------------------------------------------------

test_that("the accessor stored-slot read is load-bearing (mutation bites)", {
	skip_on_cran()
	sim <- .fm236_sim()
	fit_es <- .fm236_fit(sim, fusion_structure = "event_study")
	D_N <- .fm236_event_study_DN(fit_es)
	fit_custom <- .fm236_fit(
		sim,
		fusion_structure = "cohort",
		fusion_matrix = D_N
	)

	es_custom <- eventStudy(fit_custom)

	# Mutated accessor: a fit identical to fit_custom but with the stored slot
	# emptied to NULL (kept present so the validator's slot-presence contract
	# still holds -- `[<-` with list(NULL) preserves the name). With the value
	# NULL the choke point falls back to fit_custom$fusion_structure ("cohort"),
	# i.e. exactly what an accessor that ignored x$internal$d_inv_treat would do.
	fit_mut <- fit_custom
	fit_mut$internal["d_inv_treat"] <- list(NULL)
	es_mut <- suppressWarnings(eventStudy(fit_mut))

	# At least one SE must differ -- otherwise test 2's accessor assertions
	# would pass even with the read removed (a vacuous test).
	expect_false(isTRUE(all.equal(es_custom$se, es_mut$se)))
})

# ------------------------------------------------------------------------------
# 3. Validation errors. Each malformed fusion_matrix stops with the documented
#    message. Tested both directly against the helper (cheap, runs on CRAN) and
#    -- for the non-square case -- end-to-end through fetwfe().
# ------------------------------------------------------------------------------

test_that(".validate_fusion_matrix() rejects malformed matrices", {
	nt <- 9L
	G <- 3L
	T <- 5L
	D <- diag(nt)

	# non-square
	expect_error(
		fetwfe:::.validate_fusion_matrix(D[, -1], nt, G, T),
		"num_treats x num_treats"
	)
	# wrong num_treats (square but wrong dimension)
	expect_error(
		fetwfe:::.validate_fusion_matrix(diag(nt + 1L), nt, G, T),
		"expected num_treats = 9"
	)
	# non-numeric
	expect_error(
		fetwfe:::.validate_fusion_matrix(
			matrix("a", nt, nt),
			nt,
			G,
			T
		),
		"finite numeric matrix"
	)
	# non-finite
	bad <- D
	bad[1, 1] <- NA_real_
	expect_error(
		fetwfe:::.validate_fusion_matrix(bad, nt, G, T),
		"finite numeric matrix"
	)
	# singular
	sing <- matrix(0, nt, nt)
	sing[1, ] <- 1
	expect_error(
		fetwfe:::.validate_fusion_matrix(sing, nt, G, T),
		"singular"
	)
	# NULL passes through unchanged
	expect_null(fetwfe:::.validate_fusion_matrix(NULL, nt, G, T))
	# valid identity matrix returns its inverse (identity)
	expect_identical(
		fetwfe:::.validate_fusion_matrix(D, nt, G, T),
		solve(D)
	)
})

test_that("fetwfe() stops end-to-end on a malformed fusion_matrix", {
	skip_on_cran()
	sim <- .fm236_sim()
	fit <- .fm236_fit(sim)
	nt <- length(fit$treat_inds)

	expect_error(
		.fm236_fit(sim, fusion_matrix = diag(nt)[, -1]),
		"num_treats x num_treats"
	)
	sing <- matrix(0, nt, nt)
	sing[1, ] <- 1
	expect_error(
		.fm236_fit(sim, fusion_matrix = sing),
		"singular"
	)
})

# ------------------------------------------------------------------------------
# 4. Conditioning. The svd singular-value-lemma policing is gone: under fixed-dim
#    Assumption (D) any finite invertible D_N inherits the guarantees, so a
#    well-conditioned D_N never warns however large its singular values, and only
#    a numerically near-singular (but invertible) D_N warns -- a heads-up, not an
#    error -- while still returning its inverse / a finite fit. (nt = 9, G = 3,
#    T = 5 is the consistent triple num_treats = T*G - G*(G + 1)/2 = 9.)
# ------------------------------------------------------------------------------

test_that("a well-conditioned fusion_matrix never warns, whatever its scale", {
	# sigma_max = 10 (>> sqrt(6)) but perfectly conditioned (rcond = 1). Under the
	# removed svd bound this warned; it must now be silent. Direct, CRAN-cheap.
	nt <- 9L
	G <- 3L
	T <- 5L
	expect_warning(
		fetwfe:::.validate_fusion_matrix(diag(nt) * 10, nt, G, T),
		regexp = NA
	)
})

test_that("a near-singular (but invertible) fusion_matrix warns, not stops", {
	# Invertible but ill-conditioned: solve() succeeds (rcond ~ 1e-12, above
	# machine epsilon) yet the inverse is numerically fragile, so the validator
	# warns rather than erroring -- and still returns solve(near_singular).
	nt <- 9L
	G <- 3L
	T <- 5L
	near_singular <- diag(nt)
	near_singular[1, 1] <- 1e-12
	expect_warning(
		fetwfe:::.validate_fusion_matrix(near_singular, nt, G, T),
		"near-singular"
	)
	expect_equal(
		suppressWarnings(
			fetwfe:::.validate_fusion_matrix(near_singular, nt, G, T)
		),
		solve(near_singular)
	)
})

test_that("a non-default but well-conditioned fusion_matrix still fits end-to-end", {
	skip_on_cran()
	sim <- .fm236_sim()
	fit <- .fm236_fit(sim)
	nt <- length(fit$treat_inds)

	# 10 * I is a valid (if unusual) restriction matrix: a finite fit, no warning.
	fit_w <- .fm236_fit(sim, fusion_matrix = diag(nt) * 10)
	expect_true(is.finite(fit_w$att_hat))
	expect_equal(dim(fit_w$fusion_matrix), c(nt, nt))
})

# ------------------------------------------------------------------------------
# 5. add_ridge EXACT reproduction. add_ridge = TRUE is the only path that builds
#    the FULL inverse-fusion matrix (genFullInvFusionTransformMat), so this
#    proves the custom block reaches the ridge-penalty rows: a custom
#    event-study D_N with add_ridge = TRUE reproduces the add_ridge = TRUE
#    event_study fit byte-identically.
# ------------------------------------------------------------------------------

test_that("add_ridge = TRUE with a custom event-study D_N reproduces the add_ridge event_study fit", {
	skip_on_cran()
	sim <- .fm236_sim()

	fit_es <- .fm236_fit(
		sim,
		fusion_structure = "event_study",
		add_ridge = TRUE
	)
	D_N <- .fm236_event_study_DN(fit_es)
	fit_custom <- .fm236_fit(
		sim,
		fusion_structure = "cohort",
		fusion_matrix = D_N,
		add_ridge = TRUE
	)

	expect_identical(fit_custom$beta_hat, fit_es$beta_hat)
	expect_identical(fit_custom$att_hat, fit_es$att_hat)
	expect_identical(fit_custom$att_se, fit_es$att_se)
	expect_identical(fit_custom$catt_df, fit_es$catt_df)
})

# ------------------------------------------------------------------------------
# 6. Composition with se_type and lambda_selection. A well-formed custom D_N
#    composes with each se_type and lambda_selection (finite, well-formed fit).
# ------------------------------------------------------------------------------

test_that("a custom fusion_matrix composes with each se_type", {
	skip_on_cran()
	sim <- .fm236_sim()
	fit <- .fm236_fit(sim)
	D_N <- .fm236_event_study_DN(fit)

	for (st in c("default", "conservative", "cluster")) {
		f <- suppressWarnings(fetwfeWithSimulatedData(
			sim,
			q = 0.5,
			lambda_selection = "bic",
			se_type = st,
			fusion_matrix = D_N
		))
		expect_true(is.finite(f$att_hat), info = st)
		expect_equal(dim(f$fusion_matrix), dim(D_N), info = st)
		expect_identical(f$internal$d_inv_treat, solve(D_N), info = st)
	}
})

test_that("a custom fusion_matrix composes with each lambda_selection", {
	skip_on_cran()
	sim <- .fm236_sim()
	fit <- .fm236_fit(sim)
	D_N <- .fm236_event_study_DN(fit)

	for (ls in c("cv", "bic")) {
		f <- suppressWarnings(fetwfeWithSimulatedData(
			sim,
			q = 0.5,
			lambda_selection = ls,
			fusion_matrix = D_N
		))
		expect_true(is.finite(f$att_hat), info = ls)
		expect_identical(f$internal$d_inv_treat, solve(D_N), info = ls)
	}
})
