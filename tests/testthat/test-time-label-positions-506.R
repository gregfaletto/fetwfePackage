library(testthat)
library(fetwfe)

# Issue #506: the post-fit accessors map time labels to period positions
# through the sorted time values the fit saw (`.fit_time_values()`), so a fit
# gives the same eventStudy(), cohortTimeATTs(), simultaneousCIs() and catt_df
# however its time variable is coded, and cohortTimeATTs() and the cohort family
# of simultaneousCIs() label with the panel's own time values.
#
# One synthetic panel is fitted with its time column coded 1..8 (the reference)
# and recoded three ways: in steps of 2, with a gap, and as YYYYMM across a year
# boundary. Every recoded fit is compared with the reference fit.

.coefs506 <- genCoefs(
	G = 3,
	T = 8,
	d = 1,
	density = 0.5,
	eff_size = 2,
	seed = 7
)
.sim506 <- simulateData(
	.coefs506,
	N = 150,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 8
)
.codes506 <- list(
	steps2 = as.integer(2000 + 2 * (0:7)),
	gap = as.integer(c(2000:2002, 2004:2008)),
	yyyymm = as.integer(c(199006:199012, 199101))
)
.est506 <- list(
	fetwfe = fetwfe,
	etwfe = etwfe,
	betwfe = betwfe,
	twfeCovs = twfeCovs
)

# The panel with its k-th smallest time value replaced by `codes[k]`.
.recode506 <- function(pd, time_var, codes) {
	u <- sort(unique(pd[[time_var]]))
	stopifnot(length(codes) == length(u))
	pd[[time_var]] <- codes[match(pd[[time_var]], u)]
	pd
}

# Each fit is computed once and cached; a fit that errored re-raises its error
# in every block that asks for it.
.fits506 <- new.env()
.cached506 <- function(key, expr) {
	if (!exists(key, envir = .fits506, inherits = FALSE)) {
		assign(key, tryCatch(expr, error = function(e) e), envir = .fits506)
	}
	fit <- get(key, envir = .fits506)
	if (inherits(fit, "error")) {
		stop(fit)
	}
	fit
}

# `est_name` fitted on the synthetic panel under `coding` ("ref" is 1..8).
.fit506 <- function(
	est_name,
	coding,
	tag = "",
	pd = .sim506$pdata,
	covs = .sim506$covs,
	...
) {
	.cached506(paste(est_name, coding, tag, sep = "|"), {
		if (coding != "ref") {
			pd <- .recode506(pd, "time", .codes506[[coding]])
		}
		.est506[[est_name]](
			pdata = pd,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = covs,
			...
		)
	})
}

# A data frame's columns, without its class, attributes or the `drop` columns.
.plain506 <- function(d, drop = character(0)) {
	data.frame(
		unclass(d)[setdiff(names(d), drop)],
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}

# Columns that hold a simultaneous band or its p-values, which pass through
# mvtnorm's quasi-Monte Carlo integration: another BLAS can move them past the
# default tolerance (#541), so they compare at 1e-4, every other column at the
# default.
.band_cols506 <- c(
	"ci_low",
	"ci_high",
	"p_value",
	"simultaneous_ci_low",
	"simultaneous_ci_high"
)
.expect_frame_agrees506 <- function(alt, ref, drop = character(0), ...) {
	expect_identical(names(alt), names(ref), ...)
	band <- intersect(names(ref), .band_cols506)
	not_band <- setdiff(names(ref), band)
	expect_equal(
		.plain506(alt, c(drop, band)),
		.plain506(ref, c(drop, band)),
		...
	)
	expect_equal(
		.plain506(alt, not_band),
		.plain506(ref, not_band),
		...,
		tolerance = 1e-4
	)
}

# The checks shared by every fixture, each comparing a recoded fit `alt` with
# its reference fit `ref`.
.expect_event_study_agrees506 <- function(alt, ref) {
	es_alt <- eventStudy(alt)
	es_ref <- eventStudy(ref)
	.expect_frame_agrees506(es_alt, es_ref)
	expect_identical(
		attr(es_alt, "band_applied"),
		attr(es_ref, "band_applied")
	)
}
.expect_catt_df_agrees506 <- function(alt, ref) {
	expect_identical(names(alt$catt_df), names(ref$catt_df))
	for (col in setdiff(names(ref$catt_df), "cohort")) {
		expect_equal(
			alt$catt_df[[col]],
			ref$catt_df[[col]],
			info = col,
			tolerance = if (col %in% .band_cols506) {
				1e-4
			} else {
				testthat_tolerance()
			}
		)
	}
	expect_identical(alt$catt_band_applied, ref$catt_band_applied)
}
# The cohort family's `effect` labels differ across codings by design and have
# blocks of their own, so they are left out.
.expect_sci_agrees506 <- function(alt, ref, families) {
	for (fam in families) {
		sci_alt <- simultaneousCIs(alt, family = fam)
		sci_ref <- simultaneousCIs(ref, family = fam)
		drop <- if (fam == "cohort") "effect" else character(0)
		.expect_frame_agrees506(sci_alt$ci, sci_ref$ci, drop = drop, info = fam)
		expect_equal(
			sci_alt$adjusted_p_values,
			sci_ref$adjusted_p_values,
			info = fam,
			tolerance = 1e-4
		)
		expect_equal(
			sci_alt$critical_value,
			sci_ref$critical_value,
			info = fam,
			tolerance = 1e-4
		)
	}
}
.expect_cta_agrees506 <- function(alt, ref) {
	cta_alt <- cohortTimeATTs(alt)
	cta_ref <- cohortTimeATTs(ref)
	expect_identical(names(cta_alt), names(cta_ref))
	for (col in setdiff(names(cta_ref), c("cohort", "time"))) {
		expect_equal(cta_alt[[col]], cta_ref[[col]], info = col)
	}
}
# cohortTimeATTs() labels each cell with the time value at the cell's period.
# `codes` holds the recoded panel's sorted time values and `ref_values` the
# reference panel's, in the same order.
.expect_cta_labels506 <- function(alt, ref, codes, ref_values) {
	time_alt <- cohortTimeATTs(alt)$time
	time_ref <- cohortTimeATTs(ref)$time
	expect_equal(time_alt, codes[match(time_ref, ref_values)])
	expect_true(all(time_alt %in% codes))
}
# The cohort family labels each cohort with the label `catt_df` carries.
.expect_cohort_family_labels506 <- function(fit) {
	expect_identical(
		simultaneousCIs(fit, family = "cohort")$ci$effect,
		paste0("Cohort ", fit$catt_df$cohort)
	)
}

.all_families506 <- c("cohort", "event_study", "all_post_treatment")

for (est_name in c("fetwfe", "etwfe", "betwfe")) {
	for (coding in names(.codes506)) {
		lab <- paste0(est_name, " ", coding, ": ")

		test_that(paste0(lab, "the fit succeeds"), {
			expect_s3_class(.fit506(est_name, coding), est_name)
		})

		test_that(paste0(lab, "eventStudy() agrees with the 1..8 fit"), {
			.expect_event_study_agrees506(
				.fit506(est_name, coding),
				.fit506(est_name, "ref")
			)
		})

		test_that(paste0(lab, "catt_df agrees with the 1..8 fit"), {
			.expect_catt_df_agrees506(
				.fit506(est_name, coding),
				.fit506(est_name, "ref")
			)
		})

		test_that(paste0(lab, "simultaneousCIs() agrees with the 1..8 fit"), {
			# For etwfe() and betwfe() the all_post_treatment family reads no
			# cohort offsets, so it is compared for fetwfe() only.
			.expect_sci_agrees506(
				.fit506(est_name, coding),
				.fit506(est_name, "ref"),
				if (est_name == "fetwfe") {
					.all_families506
				} else {
					c("cohort", "event_study")
				}
			)
		})

		test_that(paste0(lab, "cohortTimeATTs() agrees with the 1..8 fit"), {
			.expect_cta_agrees506(
				.fit506(est_name, coding),
				.fit506(est_name, "ref")
			)
		})

		test_that(
			paste0(lab, "cohortTimeATTs() labels cells with time values"),
			{
				.expect_cta_labels506(
					.fit506(est_name, coding),
					.fit506(est_name, "ref"),
					codes = .codes506[[coding]],
					ref_values = 1:8
				)
			}
		)

		test_that(
			paste0(lab, "the cohort family labels cohorts by time value"),
			{
				.expect_cohort_family_labels506(.fit506(est_name, coding))
			}
		)
	}
}

for (coding in names(.codes506)) {
	test_that(
		paste0(
			"twfeCovs ",
			coding,
			": the cohort family labels cohorts by time value"
		),
		{
			.expect_cohort_family_labels506(.fit506("twfeCovs", coding))
		}
	)
}

test_that("event-study fusion, gap coding: the accessors agree with the 1..8 fit", {
	ref <- .fit506(
		"fetwfe",
		"ref",
		tag = "es_fusion",
		fusion_structure = "event_study"
	)
	alt <- .fit506(
		"fetwfe",
		"gap",
		tag = "es_fusion",
		fusion_structure = "event_study"
	)
	.expect_event_study_agrees506(alt, ref)
	.expect_catt_df_agrees506(alt, ref)
	.expect_sci_agrees506(alt, ref, .all_families506)
	.expect_cta_agrees506(alt, ref)
})

test_that("event-study fusion, gap coding: cohortTimeATTs() and the cohort family label by time value", {
	ref <- .fit506(
		"fetwfe",
		"ref",
		tag = "es_fusion",
		fusion_structure = "event_study"
	)
	alt <- .fit506(
		"fetwfe",
		"gap",
		tag = "es_fusion",
		fusion_structure = "event_study"
	)
	.expect_cta_labels506(alt, ref, codes = .codes506$gap, ref_values = 1:8)
	.expect_cohort_family_labels506(alt)
})

# ----------------------------------------------------------------------------
# A panel whose cohorts adopt at periods 3, 5 and 6 of 7, plus never-treated
# units, built as `test-event-study-present-in-print-summary-174.R` builds it.
# Its cohorts' positions are not the consecutive 2, 3, 4.
# ----------------------------------------------------------------------------

# 30 units per entry of `adopt_periods` (`NA` is never-treated), over
# `T_panel` periods; the defaults give the panel above.
.scattered_panel506 <- function(
	adopt_periods = c(3L, 5L, 6L, NA),
	T_panel = 7L
) {
	set.seed(174)
	make_cohort <- function(ids, adopt) {
		do.call(
			rbind,
			lapply(ids, function(uid) {
				x1 <- rnorm(1)
				unit_fe <- rnorm(1, sd = 0.5)
				treated <- if (is.na(adopt)) {
					integer(T_panel)
				} else {
					as.integer(seq_len(T_panel) >= adopt)
				}
				y <- unit_fe +
					0.1 * seq_len(T_panel) +
					0.6 * x1 +
					1.5 * treated +
					rnorm(T_panel, sd = 0.3)
				data.frame(
					unit = as.character(uid),
					year = seq_len(T_panel),
					x1 = x1,
					treated = treated,
					y = y,
					stringsAsFactors = FALSE
				)
			})
		)
	}
	do.call(
		rbind,
		Map(
			function(k, adopt) make_cohort((k - 1L) * 30L + 1:30, adopt),
			seq_along(adopt_periods),
			adopt_periods
		)
	)
}
.scattered_codes506 <- as.integer(2000 + 2 * (0:6))

# `etwfe()` on the scattered panel, with its years coded 1..7 ("ref") or in
# steps of 2 ("steps2").
.fit_scattered506 <- function(coding) {
	.cached506(paste("scattered", coding, sep = "|"), {
		pd <- .scattered_panel506()
		if (coding != "ref") {
			pd <- .recode506(pd, "year", .scattered_codes506)
		}
		etwfe(
			pdata = pd,
			time_var = "year",
			unit_var = "unit",
			treatment = "treated",
			response = "y",
			covs = "x1"
		)
	})
}

test_that("scattered cohorts in steps of 2: eventStudy() agrees with the 1..7 fit", {
	.expect_event_study_agrees506(
		.fit_scattered506("steps2"),
		.fit_scattered506("ref")
	)
})

test_that("scattered cohorts in steps of 2: simultaneousCIs() agrees with the 1..7 fit", {
	.expect_sci_agrees506(
		.fit_scattered506("steps2"),
		.fit_scattered506("ref"),
		.all_families506
	)
})

test_that("scattered cohorts in steps of 2: catt_df agrees with the 1..7 fit, with the band applied", {
	ref <- .fit_scattered506("ref")
	alt <- .fit_scattered506("steps2")
	.expect_catt_df_agrees506(alt, ref)
	expect_true(ref$catt_band_applied)
	expect_true(alt$catt_band_applied)
})

# ----------------------------------------------------------------------------
# The same builder with cohorts adopting at periods 3, 5, 6 and 8 of 8 and no
# never-treated units, so the fit truncates the panel
# (`allow_no_never_treated`).
# ----------------------------------------------------------------------------

.truncated_codes506 <- as.integer(2000 + 2 * (0:7))

.fit_truncated506 <- function(coding) {
	.cached506(paste("truncated", coding, sep = "|"), {
		pd <- .scattered_panel506(c(3L, 5L, 6L, 8L), T_panel = 8L)
		if (coding != "ref") {
			pd <- .recode506(pd, "year", .truncated_codes506)
		}
		# The fits warn that they truncate the panel.
		suppressWarnings(etwfe(
			pdata = pd,
			time_var = "year",
			unit_var = "unit",
			treatment = "treated",
			response = "y",
			covs = "x1"
		))
	})
}

test_that("no never-treated units, in steps of 2: eventStudy() and cohortTimeATTs() agree with the 1..8 fit, with the band applied", {
	ref <- .fit_truncated506("ref")
	alt <- .fit_truncated506("steps2")
	.expect_event_study_agrees506(alt, ref)
	.expect_cta_agrees506(alt, ref)
	expect_true(ref$catt_band_applied)
	expect_true(alt$catt_band_applied)
})

# ----------------------------------------------------------------------------
# bacondecomp::divorce, the paper's data application, with its years as given
# ("ref") and in steps of 2 ("steps2"). Its fits are slow, so its blocks skip
# on CRAN.
# ----------------------------------------------------------------------------

.divorce_panel506 <- function() {
	data("divorce", package = "bacondecomp", envir = environment())
	divorce[divorce$sex == 2, ]
}
.divorce_years506 <- function() {
	sort(unique(.divorce_panel506()$year))
}
.divorce_codes506 <- function() {
	yrs <- .divorce_years506()
	as.integer(yrs[1] + 2L * (seq_along(yrs) - 1L))
}

# `fetwfe()` with the paper's call (`q = 0.5` and three covariates), or
# `etwfe()` without covariates.
.fit_divorce506 <- function(est_name, coding) {
	.cached506(paste("divorce", est_name, coding, sep = "|"), {
		pd <- .divorce_panel506()
		if (coding != "ref") {
			pd <- .recode506(pd, "year", .divorce_codes506())
		}
		# The divorce fits warn about the units and covariates they drop.
		suppressWarnings(
			if (est_name == "fetwfe") {
				fetwfe(
					pdata = pd,
					time_var = "year",
					unit_var = "st",
					treatment = "changed",
					covs = c("murderrate", "lnpersinc", "afdcrolls"),
					response = "suiciderate_elast_jag",
					q = 0.5
				)
			} else {
				etwfe(
					pdata = pd,
					time_var = "year",
					unit_var = "st",
					treatment = "changed",
					response = "suiciderate_elast_jag"
				)
			}
		)
	})
}

test_that("divorce in steps of 2, fetwfe() paper call: the accessors agree with the yearly fit", {
	skip_on_cran()
	skip_if_not_installed("bacondecomp")
	ref <- .fit_divorce506("fetwfe", "ref")
	alt <- .fit_divorce506("fetwfe", "steps2")
	.expect_event_study_agrees506(alt, ref)
	.expect_cta_agrees506(alt, ref)
	.expect_catt_df_agrees506(alt, ref)
	.expect_sci_agrees506(alt, ref, c("cohort", "event_study"))
	expect_error(capture.output(print(alt)), NA)
	expect_error(capture.output(print(summary(alt))), NA)
})

test_that("divorce, fetwfe() paper call: cohortTimeATTs() and the cohort family label by time value", {
	skip_on_cran()
	skip_if_not_installed("bacondecomp")
	ref <- .fit_divorce506("fetwfe", "ref")
	alt <- .fit_divorce506("fetwfe", "steps2")
	.expect_cta_labels506(
		alt,
		ref,
		codes = .divorce_codes506(),
		ref_values = .divorce_years506()
	)
	.expect_cohort_family_labels506(ref)
	.expect_cohort_family_labels506(alt)
	expect_identical(
		simultaneousCIs(ref, family = "cohort")$ci$effect[1],
		"Cohort 1969"
	)
})

test_that("divorce in steps of 2, etwfe() without covariates: eventStudy() and cohortTimeATTs() agree with the yearly fit", {
	skip_on_cran()
	skip_if_not_installed("bacondecomp")
	ref <- .fit_divorce506("etwfe", "ref")
	alt <- .fit_divorce506("etwfe", "steps2")
	.expect_event_study_agrees506(alt, ref)
	.expect_cta_agrees506(alt, ref)
})

# ----------------------------------------------------------------------------
# `.fit_time_values()` directly.
# ----------------------------------------------------------------------------

test_that("a covariate named like a time dummy leaves the map intact (steps of 2)", {
	pd <- .sim506$pdata
	names(pd)[names(pd) == "cov1"] <- "t_2004"
	for (est_name in c("fetwfe", "etwfe")) {
		ref <- .fit506(
			est_name,
			"ref",
			tag = "t_2004",
			pd = pd,
			covs = "t_2004"
		)
		alt <- .fit506(
			est_name,
			"steps2",
			tag = "t_2004",
			pd = pd,
			covs = "t_2004"
		)
		.expect_frame_agrees506(
			eventStudy(alt),
			eventStudy(ref),
			info = est_name
		)
		expect_equal(
			cohortTimeATTs(alt)$time,
			.codes506$steps2[cohortTimeATTs(ref)$time],
			info = est_name
		)
	}
})

test_that(".fit_time_values() returns the sorted time values the fit saw, or NULL", {
	fit <- .fit506("etwfe", "steps2")
	expect_identical(.fit_time_values(fit), .codes506$steps2)

	# `fit` with `f()` applied to the names of its time-dummy block, in both
	# copies of the design an etwfe() fit carries.
	rename_block <- function(f) {
		x <- fit
		block <- x$G + seq_len(x$T - 1L)
		cn <- colnames(x$X_ints)
		cn[block] <- f(cn[block])
		colnames(x$X_ints) <- cn
		colnames(x$internal$X_ints) <- cn
		x
	}

	unnamed <- fit
	colnames(unnamed$X_ints) <- NULL
	colnames(unnamed$internal$X_ints) <- NULL
	expect_null(.fit_time_values(unnamed))
	# Two block names swapped, so the values are not increasing.
	expect_null(.fit_time_values(rename_block(function(nm) {
		replace(nm, 1:2, nm[2:1])
	})))
	# One block name without its `t_` prefix.
	expect_null(.fit_time_values(rename_block(function(nm) {
		replace(nm, 1L, sub("^t_", "", nm[1L]))
	})))
	# One block name whose suffix is not an integer.
	expect_null(.fit_time_values(rename_block(function(nm) {
		replace(nm, 1L, "t_x")
	})))
})

test_that("without time-dummy names, a fit coded 2001..2008 takes the arithmetic fallback", {
	fit <- .fit506(
		"etwfe",
		"ref",
		tag = "2001",
		pd = .recode506(.sim506$pdata, "time", 2001:2008)
	)
	unnamed <- fit
	colnames(unnamed$X_ints) <- NULL
	colnames(unnamed$internal$X_ints) <- NULL
	expect_null(.fit_time_values(unnamed))
	expect_identical(.derive_cohort_offsets_from_fit(unnamed), 2:4)
	expect_equal(cohortTimeATTs(unnamed)$time, cohortTimeATTs(fit)$time)
})
