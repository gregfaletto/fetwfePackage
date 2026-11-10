# eventStudy() and simultaneousCIs(family = "event_study") report only the event
# times at least one cohort reaches (#515), on the issue's panel: the simulation
# of test-plot-methods-29.R's `.plot_setup()`, with the units first treated in
# period 2 dropped.

.esu515_cache <- new.env(parent = emptyenv())

.esu515_memo <- function(key, build) {
	if (!exists(key, envir = .esu515_cache, inherits = FALSE)) {
		assign(key, build(), envir = .esu515_cache)
	}
	get(key, envir = .esu515_cache, inherits = FALSE)
}

# The panel, its number of periods `T`, and the event times its cohorts reach,
# computed from the adoption periods rather than from a fit: a cohort adopting
# in period `a` reaches event times 0 to `T - a`.
.esu515_panel <- function() {
	.esu515_memo("panel", function() {
		set.seed(2026)
		coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
		sim <- simulateData(
			coefs,
			N = 60,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = NA
		)
		pd <- sim$pdata
		treated <- pd$treatment == 1
		first_treated <- tapply(pd$time[treated], pd$unit[treated], min)
		pd <- pd[!pd$unit %in% names(first_treated)[first_treated == 2], ]

		times <- sort(unique(pd$time))
		T <- length(times)
		treated <- pd$treatment == 1
		adopt <- match(tapply(pd$time[treated], pd$unit[treated], min), times)
		list(pdata = pd, T = T, reached = 0:(T - min(adopt)))
	})
}

.esu515_fit <- function(est, se_type = "default") {
	.esu515_memo(paste(est, se_type), function() {
		fn <- switch(est, fetwfe = fetwfe, etwfe = etwfe, betwfe = betwfe)
		fn(
			pdata = .esu515_panel()$pdata,
			time_var = "time",
			unit_var = "unit",
			treatment = "treatment",
			response = "y",
			covs = c("cov1", "cov2"),
			se_type = se_type
		)
	})
}

test_that("the issue's panel has an event time no cohort reaches (#515)", {
	panel <- .esu515_panel()
	expect_lt(length(panel$reached), panel$T - 1L)
})

test_that(".reached_event_times() returns the event times some cohort reaches", {
	expect_identical(fetwfe:::.reached_event_times(2:6, 6L), 0:4)
	expect_identical(fetwfe:::.reached_event_times(6L, 6L), 0L)
	expect_identical(fetwfe:::.reached_event_times(c(3L, 4L), 6L), 0:3)
})

for (est in c("fetwfe", "etwfe", "betwfe")) {
	# Pins: values computed from the panel.

	test_that(
		paste0("eventStudy() reports the reached event times (", est, ")"),
		{
			es <- eventStudy(.esu515_fit(est))
			expect_identical(es$event_time, .esu515_panel()$reached)
			expect_true(all(es$n_cohorts >= 1L))
		}
	)

	test_that(
		paste0(
			"eventStudy(ci_type = 'pointwise') reports the reached event times (",
			est,
			")"
		),
		{
			es <- eventStudy(.esu515_fit(est), ci_type = "pointwise")
			expect_identical(es$event_time, .esu515_panel()$reached)
		}
	)

	test_that(
		paste0("the event-study family's K is the reached count (", est, ")"),
		{
			sc <- simultaneousCIs(.esu515_fit(est), family = "event_study")
			expect_identical(sc$K, length(.esu515_panel()$reached))
		}
	)

	test_that(
		paste0(
			"the event-study family labels the reached event times (",
			est,
			")"
		),
		{
			sc <- simultaneousCIs(.esu515_fit(est), family = "event_study")
			expect_identical(sc$ci$effect, paste0("e", .esu515_panel()$reached))
		}
	)

	test_that(
		paste0(
			"the conservative critical value is Bonferroni over the reached ",
			"event times (",
			est,
			")"
		),
		{
			sc <- suppressMessages(simultaneousCIs(
				.esu515_fit(est, "conservative"),
				family = "event_study"
			))
			expect_equal(
				sc$critical_value,
				stats::qnorm(1 - 0.05 / (2 * length(.esu515_panel()$reached)))
			)
		}
	)

	# Parity checks: the two surfaces agree row for row.

	test_that(
		paste0("eventStudy() applies the simultaneous band (", est, ")"),
		{
			es <- eventStudy(.esu515_fit(est))
			expect_identical(attr(es, "band_applied"), TRUE)
		}
	)

	test_that(
		paste0(
			"the event-study family's K is eventStudy()'s row count (",
			est,
			")"
		),
		{
			fit <- .esu515_fit(est)
			sc <- simultaneousCIs(fit, family = "event_study")
			expect_identical(sc$K, nrow(eventStudy(fit)))
		}
	)

	test_that(
		paste0(
			"the event-study family's labels are eventStudy()'s event times (",
			est,
			")"
		),
		{
			fit <- .esu515_fit(est)
			sc <- simultaneousCIs(fit, family = "event_study")
			expect_identical(
				sc$ci$effect,
				paste0("e", eventStudy(fit)$event_time)
			)
		}
	)

	test_that(
		paste0(
			"the event-study family's estimates are eventStudy()'s (",
			est,
			")"
		),
		{
			fit <- .esu515_fit(est)
			sc <- simultaneousCIs(fit, family = "event_study")
			expect_equal(sc$ci$estimate, eventStudy(fit)$estimate)
		}
	)

	test_that(
		paste0(
			"the conservative critical value is Bonferroni over the family's K (",
			est,
			")"
		),
		{
			sc <- suppressMessages(simultaneousCIs(
				.esu515_fit(est, "conservative"),
				family = "event_study"
			))
			expect_equal(sc$critical_value, stats::qnorm(1 - 0.05 / (2 * sc$K)))
		}
	)

	test_that(
		paste0(
			"plot() draws a point at each reached event time only (",
			est,
			")"
		),
		{
			skip_if_not_installed("ggplot2")
			p <- plot(.esu515_fit(est))
			is_point <- vapply(
				p$layers,
				function(l) inherits(l$geom, "GeomPoint"),
				logical(1)
			)
			ld <- ggplot2::layer_data(p, which(is_point))
			expect_setequal(ld$x, .esu515_panel()$reached)
		}
	)
}

test_that("the bootstrap event-study family's K is the reached count (fetwfe)", {
	sc <- simultaneousCIs(
		.esu515_fit("fetwfe"),
		family = "event_study",
		method = "bootstrap",
		B = 200,
		seed = 1
	)
	expect_identical(sc$K, length(.esu515_panel()$reached))
})

test_that("the bootstrap event-study family labels the reached event times (fetwfe)", {
	sc <- simultaneousCIs(
		.esu515_fit("fetwfe"),
		family = "event_study",
		method = "bootstrap",
		B = 200,
		seed = 1
	)
	expect_identical(sc$ci$effect, paste0("e", .esu515_panel()$reached))
})
