# etwfe() returns class "fetwfe_etwfe" (#523), so that fetwfe leaves the etwfe
# package's fits, of class c("etwfe", "fixest") or c("etwfe", "fixest_multi"),
# to their own methods. The stand-ins below carry those two class vectors.
#
# Each dispatch block defines its stubs beside the generics it calls; see #537.

# The refit sentence. Blocks that expect it pin it whole; blocks that expect
# none look for its last words.
.refit523 <- "must be refit"
.refit_sentence523 <- paste0(
	"An object of class \"etwfe\" is either a fit from the etwfe package or an ",
	"etwfe() fit saved by an earlier version of fetwfe, which ",
	.refit523,
	"."
)

.sim523 <- function() {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 1,
		seed = 523
	)
	simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 523
	)
}

test_that("an etwfe-package fit of class (etwfe, fixest) reaches its own methods (#523)", {
	m <- structure(list(), class = c("etwfe", "fixest"))
	coef.fixest <- function(object, ...) "stub:coef.fixest"
	print.fixest <- function(x, ...) "stub:print.fixest"
	summary.fixest <- function(object, ...) "stub:summary.fixest"
	tidy.fixest <- function(x, ...) "stub:tidy.fixest"
	glance.fixest <- function(x, ...) "stub:glance.fixest"
	augment.fixest <- function(x, ...) "stub:augment.fixest"
	plot.fixest <- function(x, ...) "stub:plot.fixest"
	expect_identical(
		tryCatch(coef(m), error = conditionMessage),
		"stub:coef.fixest"
	)
	expect_identical(
		tryCatch(print(m), error = conditionMessage),
		"stub:print.fixest"
	)
	expect_identical(
		tryCatch(summary(m), error = conditionMessage),
		"stub:summary.fixest"
	)
	expect_identical(
		tryCatch(generics::tidy(m), error = conditionMessage),
		"stub:tidy.fixest"
	)
	expect_identical(
		tryCatch(generics::glance(m), error = conditionMessage),
		"stub:glance.fixest"
	)
	expect_identical(
		tryCatch(generics::augment(m), error = conditionMessage),
		"stub:augment.fixest"
	)
	expect_identical(
		tryCatch(plot(m), error = conditionMessage),
		"stub:plot.fixest"
	)
})

test_that("an etwfe-package fit of class (etwfe, fixest_multi) reaches its own methods (#523)", {
	m <- structure(list(), class = c("etwfe", "fixest_multi"))
	coef.fixest_multi <- function(object, ...) "stub:coef.fixest_multi"
	print.fixest_multi <- function(x, ...) "stub:print.fixest_multi"
	summary.fixest_multi <- function(object, ...) "stub:summary.fixest_multi"
	tidy.fixest_multi <- function(x, ...) "stub:tidy.fixest_multi"
	glance.fixest_multi <- function(x, ...) "stub:glance.fixest_multi"
	augment.fixest_multi <- function(x, ...) "stub:augment.fixest_multi"
	plot.fixest_multi <- function(x, ...) "stub:plot.fixest_multi"
	expect_identical(
		tryCatch(coef(m), error = conditionMessage),
		"stub:coef.fixest_multi"
	)
	expect_identical(
		tryCatch(print(m), error = conditionMessage),
		"stub:print.fixest_multi"
	)
	expect_identical(
		tryCatch(summary(m), error = conditionMessage),
		"stub:summary.fixest_multi"
	)
	expect_identical(
		tryCatch(generics::tidy(m), error = conditionMessage),
		"stub:tidy.fixest_multi"
	)
	expect_identical(
		tryCatch(generics::glance(m), error = conditionMessage),
		"stub:glance.fixest_multi"
	)
	expect_identical(
		tryCatch(generics::augment(m), error = conditionMessage),
		"stub:augment.fixest_multi"
	)
	expect_identical(
		tryCatch(plot(m), error = conditionMessage),
		"stub:plot.fixest_multi"
	)
})

test_that("fetwfe registers etwfe()'s methods for fetwfe_etwfe and none for etwfe (#523)", {
	reg <- getNamespaceInfo("fetwfe", "S3methods")
	expect_false(any(reg[, 2] %in% c("etwfe", "summary.etwfe")))
	moved <- reg[, 2] %in% c("fetwfe_etwfe", "summary.fetwfe_etwfe")
	expect_setequal(
		paste(reg[moved, 1], reg[moved, 2], sep = "."),
		c(
			"augment.fetwfe_etwfe",
			"coef.fetwfe_etwfe",
			"glance.fetwfe_etwfe",
			"plot.fetwfe_etwfe",
			"print.fetwfe_etwfe",
			"print.summary.fetwfe_etwfe",
			"simultaneousCIs.fetwfe_etwfe",
			"summary.fetwfe_etwfe",
			"tidy.fetwfe_etwfe"
		)
	)
})

test_that("etwfe() returns class fetwfe_etwfe, with fetwfe's methods and accessors (#523)", {
	sim <- .sim523()
	fit <- etwfeWithSimulatedData(sim)
	expect_true(inherits(fit, "fetwfe_etwfe"))
	expect_false(inherits(fit, "etwfe"))
	expect_true(
		"Extended Two-Way Fixed Effects Results" %in% capture.output(print(fit))
	)
	expect_s3_class(summary(fit), "summary.fetwfe_etwfe")
	expect_true(
		"Summary of Extended Two-Way Fixed Effects" %in%
			capture.output(print(summary(fit)))
	)
	expect_true(is.numeric(coef(fit)))
	expect_s3_class(generics::tidy(fit), "data.frame")
	expect_s3_class(generics::glance(fit), "data.frame")
	expect_s3_class(generics::augment(fit, data = sim$pdata), "data.frame")
	expect_s3_class(eventStudy(fit), "data.frame")
	expect_s3_class(cohortStudy(fit), "data.frame")
	expect_s3_class(cohortTimeATTs(fit), "data.frame")
	expect_no_error(simultaneousCIs(fit))
	expect_false(any(grepl(
		"more event times",
		capture.output(print(fit)),
		fixed = TRUE
	)))
	old <- options(etwfe.max_event_times = 1)
	on.exit(options(old), add = TRUE)
	expect_true(any(grepl(
		"more event times",
		capture.output(print(fit)),
		fixed = TRUE
	)))
	skip_if_not_installed("ggplot2")
	expect_s3_class(plot(fit), "ggplot")
})

.accessors523 <- list(
	eventStudy = eventStudy,
	cohortStudy = cohortStudy,
	cohortTimeATTs = cohortTimeATTs,
	simultaneousCIs = simultaneousCIs
)

# Each accessor's first sentence, through "Got class: ".
.first523 <- c(
	eventStudy = "eventStudy() requires a fitted object from fetwfe(), etwfe(), or betwfe(). Got class: ",
	cohortStudy = "cohortStudy() requires a fitted object from fetwfe(), etwfe(), betwfe(), or twfeCovs(). Got class: ",
	cohortTimeATTs = "cohortTimeATTs() requires a fitted object from fetwfe(), etwfe(), or betwfe(). Got class: ",
	simultaneousCIs = "simultaneousCIs() requires a fitted object from fetwfe(), etwfe(), betwfe(), or twfeCovs(). Got class: "
)

test_that("the accessors reject an etwfe-package fit with the refit sentence (#523)", {
	m <- structure(list(), class = c("etwfe", "fixest"))
	for (nm in names(.accessors523)) {
		e <- tryCatch(.accessors523[[nm]](m), error = function(e) e)
		expect_true(grepl(.first523[[nm]], msg_of(e), fixed = TRUE), info = nm)
		expect_true(
			grepl("Got class: etwfe, fixest.", msg_of(e), fixed = TRUE),
			info = nm
		)
		expect_true(
			grepl(.refit_sentence523, msg_of(e), fixed = TRUE),
			info = nm
		)
		expect_null(call_of(e), info = nm)
	}
})

test_that("the accessors reject an unrelated object without the refit sentence (#523)", {
	m <- stats::lm(mpg ~ wt, data = mtcars)
	for (nm in names(.accessors523)) {
		e <- tryCatch(.accessors523[[nm]](m), error = function(e) e)
		expect_true(grepl(.first523[[nm]], msg_of(e), fixed = TRUE), info = nm)
		expect_true(grepl("Got class: lm.", msg_of(e), fixed = TRUE), info = nm)
		expect_false(grepl(.refit523, msg_of(e), fixed = TRUE), info = nm)
		expect_null(call_of(e), info = nm)
	}
})

test_that("the accessors reject a saved etwfe() fit of class etwfe with the refit sentence (#523)", {
	m <- structure(list(), class = "etwfe")
	for (nm in names(.accessors523)) {
		e <- tryCatch(.accessors523[[nm]](m), error = function(e) e)
		expect_true(grepl(.first523[[nm]], msg_of(e), fixed = TRUE), info = nm)
		expect_true(
			grepl("Got class: etwfe.", msg_of(e), fixed = TRUE),
			info = nm
		)
		expect_true(
			grepl(.refit_sentence523, msg_of(e), fixed = TRUE),
			info = nm
		)
		expect_null(call_of(e), info = nm)
	}
})

test_that("the accessors reject a plain fixest fit of class fixest without the refit sentence (#523)", {
	m <- structure(list(), class = "fixest")
	for (nm in names(.accessors523)) {
		e <- tryCatch(.accessors523[[nm]](m), error = function(e) e)
		expect_true(grepl(.first523[[nm]], msg_of(e), fixed = TRUE), info = nm)
		expect_true(
			grepl("Got class: fixest.", msg_of(e), fixed = TRUE),
			info = nm
		)
		expect_false(grepl(.refit523, msg_of(e), fixed = TRUE), info = nm)
		expect_null(call_of(e), info = nm)
	}
})

test_that(".stop_unsupported_class() lists one constructor alone and joins two with or (#523)", {
	one <- tryCatch(
		fetwfe:::.stop_unsupported_class(1L, "f", "fetwfe"),
		error = function(e) e
	)
	two <- tryCatch(
		fetwfe:::.stop_unsupported_class(1L, "f", c("fetwfe", "etwfe")),
		error = function(e) e
	)
	expect_true(grepl(
		"f() requires a fitted object from fetwfe(). Got class: ",
		msg_of(one),
		fixed = TRUE
	))
	expect_true(grepl(
		"f() requires a fitted object from fetwfe() or etwfe(). Got class: ",
		msg_of(two),
		fixed = TRUE
	))
})
