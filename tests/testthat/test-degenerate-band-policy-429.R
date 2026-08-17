# Issue #429 items 6 and 7: the #304 degenerate-band silence policy, and the
# schema of the degenerate early return, in the HIGH-DIMENSIONAL (p >= N*T) arm.
#
# Item 6. .fit_band_for_family() (R/simultaneous_cis.R:1173-1198) hardcodes
# `warn_degenerate_highdim = FALSE` at :1182, which is what keeps a plain
# fetwfe() / betwfe() call SILENT when the bridge penalty zeroes every treatment
# effect in the high-dimensional regime, while a direct simultaneousCIs() call
# WARNS. Nothing asserted the fit-time half, so flipping that FALSE to TRUE --
# which would make every high-dim degenerate fit start warning -- was invisible.
#
# Item 7. The degenerate early return's data.frame() (R/simultaneous_cis.R:
# 641-649) had no absolute column-name pin on this arm. Adding, removing or
# reordering a column there was invisible to the high-dim path.
# (test-degenerate-jacobian-225.R now carries the same pin on the fixed-p arm.)
#
# THE TWO PINS ARE REDUNDANT TODAY, and the honest reason to keep both is the
# same one test-c6-toplevel-routing-429.R:20-23 gives for its nine redundant
# cells. R/simultaneous_cis.R:623 is a SINGLE outer guard; the inner if/else
# below it chooses only between warning() and message(), and both branches fall
# through to one `ses <- rep(0, K)` and one data.frame() at :641-649. So either
# pin alone catches an add, a remove or a reorder -- measured, 4 failures from
# one and 3 from the other. They stop being redundant the moment that branch is
# split per regime, which is exactly the shape #429 exists because of. (An
# earlier draft of this comment claimed the two arms "enter the same
# data.frame() through different guards, so both are needed" -- there is one
# guard, and that reasoning was backwards.)
#
# FIXTURE FAILURE MODE, stated because CI is gating over four Linux/OpenBLAS
# jobs: `p >= N*T` is deterministic, but `all(catt_hats == 0)` is the outcome of
# a cross-validated penalized fit. If this file goes red on, say,
# `ubuntu-latest (oldrel-2)` and green everywhere else, the FIXTURE drifted out
# of the degenerate high-dim regime -- the policy did not change. The four
# invariants in the first block are what makes that read as a loud failure
# rather than as silent vacuity. (simulateData(..., seed = 55) calls set.seed()
# without restoring, so the CV folds are deterministic within a platform; the
# residual risk is BLAS, not RNG.)

# ------------------------------------------------------------------------------
# Fixture, at file scope: every block below needs it, and test_that()
# environments do not leak.
# ------------------------------------------------------------------------------
.dbp429_cf <- genCoefs(
	G = 4,
	T = 6,
	d = 12,
	density = 0.5,
	eff_size = 2,
	seed = 55
)
.dbp429_sim <- simulateData(
	.dbp429_cf,
	N = 30,
	sig_eps_sq = 1,
	sig_eps_c_sq = 0.5,
	seed = 55
)

# NO suppressWarnings() here, deliberately. Block (a) below asserts that these
# exact calls emit no warning, so suppressing one here would contradict the
# file's own claim and would hide a fixture that started warning -- the failure
# this file exists to detect. Measured: both fits are warning-free.
.dbp429_fits <- list(
	fetwfe = fetwfeWithSimulatedData(.dbp429_sim, q = 0.5),
	betwfe = betwfeWithSimulatedData(.dbp429_sim, q = 0.5)
)

.dbp429_ci_names <- c(
	"effect",
	"estimate",
	"simultaneous_ci_low",
	"simultaneous_ci_high",
	"pointwise_ci_low",
	"pointwise_ci_high"
)

# ------------------------------------------------------------------------------
# The fixture invariants get their own block, FIRST in the file, so that both of
# the blocks below run on an asserted regime. Tucking them inside block (a)
# would leave block (b)'s column pins running on an unasserted fixture, which is
# the same vacuity these assertions exist to prevent, moved one block over.
# ------------------------------------------------------------------------------
test_that("the degenerate high-dim fixture really is degenerate and high-dim (#429 items 6-7)", {
	for (nm in names(.dbp429_fits)) {
		fit <- .dbp429_fits[[nm]]
		# Measured on this fixture: p = 311, N*T = 180.
		expect_true(fit$p >= fit$N * fit$T, info = nm)
		expect_true(all(fit$catt_hats == 0), info = nm)
		# calc_ses and ci_type are the two the silence assertions actually
		# depend on: .apply_simultaneous_catt_band() returns NULL unless
		# has_valid_ses, and .finalize_ci_type() early-returns unless ci_type is
		# "simultaneous". Without these, a fixture that lost either would make
		# every silence assertion below pass for the wrong reason.
		expect_true(fit$calc_ses, info = nm)
		expect_identical(fit$ci_type, "simultaneous", info = nm)
	}
})

# ------------------------------------------------------------------------------
# (a) Item 6 -- both halves of the #304 policy: the fit is silent, the direct
# accessor warns.
# ------------------------------------------------------------------------------
test_that("the #304 policy keeps the fit silent while simultaneousCIs() warns (#429 item 6)", {
	# Fit-time silence. These are deliberately SECOND fits, not the file-scope
	# ones: expect_no_warning() has to wrap the call itself, so the fit-time
	# half of the policy necessarily costs one extra fit per estimator. Do not
	# "optimize" these into the file-scope fits -- that would silently delete
	# the fit-time half.
	expect_no_warning(fetwfeWithSimulatedData(.dbp429_sim, q = 0.5))
	expect_no_warning(betwfeWithSimulatedData(.dbp429_sim, q = 0.5))

	for (nm in names(.dbp429_fits)) {
		fit <- .dbp429_fits[[nm]]

		for (fam in c("cohort", "event_study")) {
			# The anchor is deliberately NOT "zeroed out every treatment", which
			# is a shared prefix of the warning() arm (:625-633) and the
			# message() arm (:634-639). "bridge selection is NOT " appears only
			# in the high-dimensional warning, which is the arm this file exists
			# to reach.
			expect_warning(
				simultaneousCIs(fit, family = fam),
				"bridge selection is NOT ",
				fixed = TRUE,
				info = paste(nm, fam)
			)
		}

		# The other route into the same mutated line. Do NOT write
		# expect_no_warning(eventStudy(fit)) here: measured, eventStudy() calls
		# .fit_band_for_family() ZERO times on this fixture and can never reach
		# the degenerate arm. An all-zeroed bridge makes .event_study_fetwfe()
		# force its local calc_ses <- FALSE (R/event_study.R:550), which closes
		# the :820 gate -- and that is the SAME condition as the degenerate
		# branch at R/simultaneous_cis.R:623. A fixture that reaches one cannot
		# reach the other. The fit itself calls the helper exactly once, for the
		# COHORT family only, so a direct call is the only coverage the
		# event-study family can have at all.
		#
		# expected_rows is derived, not hardcoded: nrow(catt_df) is K for the
		# cohort family, T - 1L is K for the event-study family. A wrong value
		# makes the helper return NULL rather than error, which is why the
		# non-vacuity control below matters.
		#
		# One standing caution for whoever reuses these calls:
		# .fit_band_for_family() hardcodes has_valid_ses = TRUE. That is sound
		# here only because the invariants block above pins fit$calc_ses, so
		# these fits genuinely have valid standard errors. Pointed at a
		# calc_ses = FALSE fit, the helper would assert an invariant the object
		# does not satisfy, and whatever it returned would be meaningless
		# rather than wrong-and-loud.
		expected_rows <- c(
			cohort = nrow(fit$catt_df),
			event_study = fit$T - 1L
		)
		for (fam in c("cohort", "event_study")) {
			expect_no_warning(
				fetwfe:::.fit_band_for_family(
					fit,
					fam,
					0.05,
					expected_rows[[fam]]
				)
			)
			# NON-VACUITY CONTROL, not coverage. It adds nothing against the
			# `warn_degenerate_highdim = TRUE` mutation (the return stays
			# non-NULL under it). Its job is that .fit_band_for_family() returns
			# NULL on both a tryCatch()ed error and a row-count mismatch, and
			# expect_no_warning() passes perfectly happily on a NULL.
			expect_false(
				is.null(fetwfe:::.fit_band_for_family(
					fit,
					fam,
					0.05,
					expected_rows[[fam]]
				)),
				info = paste(nm, fam)
			)
		}
	}
})

# ------------------------------------------------------------------------------
# (b) Item 7 -- the degenerate $ci schema, high-dim arm.
#
# THIS BLOCK IS NOT SELF-PROTECTING, and the dependency is declared rather than
# assumed: a NON-degenerate fit returns these same six column names from the
# ordinary path at R/simultaneous_cis.R:1014-1022, so nothing asserted below
# would notice if the fixture stopped being degenerate. What rules that out is
# the invariants block at the top of this file (`all(catt_hats == 0)` plus
# `p >= N*T`), which fails loudly on the same run. Do not move these assertions
# into a file that lacks those invariants.
# ------------------------------------------------------------------------------
test_that("the degenerate band's $ci carries exactly six columns in order (#429 item 7)", {
	for (nm in names(.dbp429_fits)) {
		for (fam in c("cohort", "event_study")) {
			# Block (a) owns the warning assertions; here the warning is noise.
			sci <- suppressWarnings(simultaneousCIs(
				.dbp429_fits[[nm]],
				family = fam
			))
			expect_identical(
				names(sci$ci),
				.dbp429_ci_names,
				info = paste(nm, fam)
			)
		}
	}
})
