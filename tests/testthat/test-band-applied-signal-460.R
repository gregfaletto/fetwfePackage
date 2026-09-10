library(testthat)
library(fetwfe)

# Issue #460: `print()` and `summary()` labelled BOTH preview blocks from
# `x$ci_type`, which records what the user ASKED FOR rather than what was
# APPLIED -- so a fit whose band came back `NULL` rendered
# `[simultaneous 95% CI]` over bounds that are not a simultaneous band.
#
# The fix adds two positive signals, one per family, because the two bands are
# built from DIFFERENT contrast matrices over the same sandwich and fail
# independently (`.workflow/PROFILE.md` section 12 gotcha 8):
#
#   fit$catt_band_applied              -- the fit-time COHORT band, written by
#                                         `.finalize_ci_type()`
#   attr(eventStudy(fit), "band_applied") -- the render-time EVENT-STUDY band,
#                                         written by `.finish_event_study()`
#
# `ci_type` is deliberately unchanged: it is documented as the user's argument.
#
# ---------------------------------------------------------------------------
# ASSERTION FORM. Two rules, both measured, both mandatory here.
#
# 1. Applied-signal assertions are written LITERALLY --
#    `expect_true(fit$catt_band_applied)`, never
#    `expect_true(isTRUE(fit$catt_band_applied))`. `isTRUE()` is the right
#    PRODUCTION idiom (a missing signal must land on "pointwise", the
#    conservative side) and exactly the wrong TEST idiom, because it collapses
#    "absent" and `FALSE` -- the two states this file exists to distinguish.
#    Measured under this suite's edition (2): `expect_false(NULL)` FAILS
#    ("Expected NULL to be FALSE"), which is what makes the literal form red on
#    a tree where the slot and the attribute do not exist.
#
# 2. `fixed = TRUE` ON EVERY RENDERED-LITERAL ASSERTION, without exception.
#    `[pointwise 95% CI]` and `[simultaneous 95% CI]` are CHARACTER CLASSES, so
#    each matches the OTHER header's line. Measured:
#
#      grepl("[pointwise 95% CI]", "...(CATT) [simultaneous 95% CI]:")   TRUE (!)
#      grepl("[pointwise 95% CI]", "...(CATT) [simultaneous 95% CI]:",
#            fixed = TRUE)                                              FALSE
#
#    and the hazard is symmetric. Without `fixed = TRUE` every header assertion
#    below -- blocks 6 and 9 included, the two that carry the change -- passes
#    on the unfixed tree. The three local helpers below bake the argument in so
#    it cannot be dropped by inattention, which is the convention
#    `tests/testthat/test-print-summary-single-source-439.R` and
#    `tests/testthat/test-class-helpers.R` already use.
#
# 3. HEADERS ARE ASSERTED PER LINE, never on the whole blob. After this change
#    the two headers routinely differ, so a blob-level presence assertion is
#    satisfied by a mutant that SWAPS them.
# ---------------------------------------------------------------------------
#
# The fixtures below are rebuilt here rather than borrowed. `.hd_fit_326` and
# `.fp_fit_326` live in `tests/testthat/test-print-highdim-diagnostics-326.R`,
# whose own assertions read them, and a fitted object cannot go in a
# `helper-*.R` (those are sourced for every test file, so the fit would be paid
# for on every run). Duplicating the construction is the cheaper trade.

.render_print <- function(x, ...) {
	capture.output(print(x, ...))
}

.render_summary <- function(x, ...) {
	capture.output(print(summary(x, ...)))
}

.expect_absent <- function(lines, literal) {
	expect_false(any(grepl(literal, lines, fixed = TRUE)))
}

# Select the one header line starting with `prefix`, then assert `literal` is on
# THAT line. `startsWith()` rather than `grep()`: it is fixed-string and
# anchored at once, so neither `fixed = TRUE` nor a `^` can go missing.
.expect_header <- function(lines, prefix, literal) {
	hits <- which(startsWith(lines, prefix))
	expect_length(hits, 1L)
	expect_match(lines[hits[1]], literal, fixed = TRUE)
}

.CATT_HEADER <- "Cohort Average Treatment Effects (CATT)"
.ES_HEADER <- "Event-Study Average Treatment Effects"
.CATT_PREVIEW <- "CATT (preview)"
.ES_PREVIEW <- "Event Study (preview)"
.POINTWISE <- "[pointwise 95% CI]"
.SIMULTANEOUS <- "[simultaneous 95% CI]"

# ------------------------------------------------------------------------------
# Fixtures.
# ------------------------------------------------------------------------------

# Cell 1, the demonstrated reproducer: a plain public `fetwfe(gls = FALSE)` call
# on a 24-unit, 8-period, 10-covariate scattered panel. Same construction as the
# `.hd_fit_326` fixture. p = 274 against N*T = 192, `calc_ses` FALSE, `catt_df`
# bounds all NA -- so NEITHER band is applied, and before #460 both preview
# headers claimed otherwise.
.bas460_hd_fit <- local({
	set.seed(5)
	d <- 10L
	N <- 24L
	T <- 8L
	adopt <- c(2L, 4L, 7L)
	covs <- matrix(stats::rnorm(N * d), N, d)
	cou <- c(
		rep(0L, N - 18L),
		rep(adopt[1], 6),
		rep(adopt[2], 6),
		rep(adopt[3], 6)
	)
	eff <- stats::setNames(c(0.5, 2, 3.5), as.character(adopt))
	panel <- do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cou[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T,
				treat = as.integer(g > 0 & (1:T) >= g)
			)
			for (j in 1:d) {
				df[[paste0("x", j)]] <- covs[i, j]
			}
			df$y <- 0.3 *
				(1:T) /
				T +
				(if (g > 0) eff[[as.character(g)]] else 0) *
					df$treat +
				0.2 * covs[i, 1] +
				stats::rnorm(T, 0, 0.4)
			df
		})
	)
	fetwfe(
		pdata = panel,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:10),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
})

# The well-conditioned fixed-p simulation behind blocks 4, 5 and 6. Same
# construction as the `.fp_fit_326` fixture: p = 50 against N*T = 600, every
# cohort SE finite and positive, so BOTH bands apply.
.bas460_fp_sim <- local({
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	simulateData(cf, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
})

.bas460_fp_fit <- function(...) {
	fetwfe(
		pdata = .bas460_fp_sim$pdata,
		time_var = .bas460_fp_sim$time_var,
		unit_var = .bas460_fp_sim$unit_var,
		treatment = .bas460_fp_sim$treatment,
		response = .bas460_fp_sim$response,
		covs = .bas460_fp_sim$covs,
		q = 0.5,
		verbose = FALSE,
		...
	)
}

# ------------------------------------------------------------------------------
# 1. (red) Cell 1, the demonstrated reproducer, on the print path. Before #460
#    BOTH headers read `[simultaneous 95% CI]` over all-NA bounds; the issue's
#    acceptance criterion is stated over this output.
# ------------------------------------------------------------------------------
test_that("a fit with neither band applied prints two pointwise headers (#460)", {
	fit <- .bas460_hd_fit

	expect_false(fit$catt_band_applied)
	expect_false(attr(eventStudy(fit), "band_applied"))
	# `ci_type` still records the REQUEST. This is what distinguishes #460's
	# resolution from downgrading `ci_type`, and it is a contract in its own
	# right: the slot is documented as the user's argument.
	expect_identical(fit$ci_type, "simultaneous")
	# Cell 1: no standard errors, so the stored bounds are all NA.
	expect_true(all(is.na(fit$catt_df$ci_low)))

	out <- .render_print(fit)
	# The issue's criterion, over the WHOLE output: not one line says
	# "simultaneous" anywhere.
	.expect_absent(out, "simultaneous")
	# ...and each header is asserted separately, by its own leading text, so a
	# fix to one is not scored as a fix to both.
	.expect_header(out, .CATT_HEADER, .POINTWISE)
	.expect_header(out, .ES_HEADER, .POINTWISE)
})

# ------------------------------------------------------------------------------
# 2. (red) The same fit through `summary()` / `print(summary())`. The summary
#    path resolves its labels from a `summary.<class>` list, not from the fit,
#    so it is a second code path and not a second view of the first.
# ------------------------------------------------------------------------------
test_that("the summary path renders the same two pointwise headers (#460)", {
	out <- .render_summary(.bas460_hd_fit)

	.expect_absent(out, "simultaneous")
	.expect_header(out, .CATT_PREVIEW, .POINTWISE)
	.expect_header(out, .ES_PREVIEW, .POINTWISE)
})

# ------------------------------------------------------------------------------
# 3. (guardrail) The summary object CARRIES the field. `.summary_estimator_output()`
#    builds `out <- list(...)` then `out <- out[keep]`, and `keep` is a
#    whitelist that drops any unlisted name silently. Omitting the field there
#    would relabel every summary as pointwise with no error anywhere, so this
#    asserts the field's presence directly rather than only through a rendering.
# ------------------------------------------------------------------------------
test_that("summary() carries catt_band_applied through the keep whitelist (#460)", {
	s <- summary(.bas460_hd_fit)
	expect_true("catt_band_applied" %in% names(s))
	expect_false(s$catt_band_applied)
	# `ci_type` stays on the summary object too -- it is part of that object's
	# shape and still reports, truthfully, what the user asked for.
	expect_identical(s$ci_type, "simultaneous")
})

# ------------------------------------------------------------------------------
# 4. (guardrail, positive control) A well-conditioned fixed-p fit: both bands
#    apply, both headers say so. WITHOUT THIS BLOCK every assertion above is
#    satisfied by a `.band_label()` hardcoded to "pointwise".
# ------------------------------------------------------------------------------
test_that("a fit with both bands applied prints two simultaneous headers (#460)", {
	fit <- .bas460_fp_fit()

	expect_true(fit$catt_band_applied)
	expect_true(attr(eventStudy(fit), "band_applied"))

	out <- .render_print(fit)
	.expect_header(out, .CATT_HEADER, .SIMULTANEOUS)
	.expect_header(out, .ES_HEADER, .SIMULTANEOUS)

	out_s <- .render_summary(fit)
	.expect_header(out_s, .CATT_PREVIEW, .SIMULTANEOUS)
	.expect_header(out_s, .ES_PREVIEW, .SIMULTANEOUS)
})

# ------------------------------------------------------------------------------
# 5. (guardrail) `ci_type = "pointwise"` is unaffected. The same fit refitted:
#    neither signal is set, both headers stay pointwise. This is pre-existing
#    behavior, pinned here so #460 cannot quietly alter it.
# ------------------------------------------------------------------------------
test_that("ci_type = 'pointwise' sets neither signal and labels both headers (#460)", {
	fit <- .bas460_fp_fit(ci_type = "pointwise")

	expect_false(fit$catt_band_applied)
	expect_false(attr(eventStudy(fit), "band_applied"))
	expect_identical(fit$ci_type, "pointwise")

	out <- .render_print(fit)
	.expect_absent(out, "simultaneous")
	.expect_header(out, .CATT_HEADER, .POINTWISE)
	.expect_header(out, .ES_HEADER, .POINTWISE)
})

# ------------------------------------------------------------------------------
# 6. (red) CELL 3, mocked: `calc_ses = TRUE`, `ci_type = "simultaneous"`, and NO
#    applied cohort band -- so `catt_df` holds the fit-time POINTWISE Wald
#    intervals, which are FINITE. That finiteness is what makes this cell 3
#    rather than cell 1, and it is why a bounds-aware label (`all(is.na(bounds))`)
#    could not have fixed #460.
#
#    The mock gates on `.apply_simultaneous_catt_band`, the FIT-TIME route
#    (PROFILE section 12 gotcha 8: a fixture for the renderer route must gate on
#    `.event_study_simultaneous_bounds` instead). That is deliberate here: the
#    event-study band is left alone precisely so the two families DISAGREE.
#
#    The disagreement assertion is load-bearing. Labelling the event-study header
#    from `x$catt_band_applied` instead of the frame's own attribute renders
#    byte-identical output on every other block in this file.
# ------------------------------------------------------------------------------
test_that("a calc_ses fit whose cohort band failed labels only that header pointwise (#460)", {
	fit <- with_mocked_bindings(
		.bas460_fp_fit(),
		.apply_simultaneous_catt_band = function(x, alpha, has_valid_ses) NULL,
		.package = "fetwfe"
	)

	expect_identical(fit$ci_type, "simultaneous")
	expect_true(fit$internal$calc_ses)
	expect_false(fit$catt_band_applied)
	# Cell 3, not cell 1: the displayed bounds exist.
	expect_true(all(is.finite(fit$catt_df$ci_low)))
	expect_true(all(is.finite(fit$catt_df$ci_high)))
	# The other family is untouched, so its band still applies.
	expect_true(attr(eventStudy(fit), "band_applied"))

	out <- .render_print(fit)
	.expect_header(out, .CATT_HEADER, .POINTWISE)
	.expect_header(out, .ES_HEADER, .SIMULTANEOUS)

	out_s <- .render_summary(fit)
	.expect_header(out_s, .CATT_PREVIEW, .POINTWISE)
	.expect_header(out_s, .ES_PREVIEW, .SIMULTANEOUS)
})

# ------------------------------------------------------------------------------
# 7. (red) Cell 3's CAVEAT TEXT. #433's rendered caveat is kept on `p >= N*T`
#    fits whose band was not applied -- the displayed intervals are the fit-time
#    pointwise Wald bounds on the selected support, which still under-cover
#    (#308) -- but its noun ("this band") is wrong when there is no band.
#
#    `skip_on_cran()` because this block builds the `p = 356` bridge fit TWICE.
#    It follows the file that owns the fixture
#    (`tests/testthat/test-highdim-postselection-band-warning-433.R`), whose
#    fit-building blocks all skip for the same reason. Per
#    `.workflow/PROFILE.md` section 3 this skips almost nowhere it will actually
#    run: `setup-r` exports `NOT_CRAN=true`, so all six CI jobs exercise it on
#    every PR.
# ------------------------------------------------------------------------------
test_that("the #433 caveat renames its subject when no band was applied (#460)", {
	skip_on_cran()

	# The `.hpb433_build("hd_betwfe")` construction: p = 356 against N*T = 300,
	# calc_ses = TRUE, high-dimensional and NOT degenerate.
	sim <- simulateData(
		genCoefs(G = 3, T = 5, d = 20, density = 0.08, eff_size = 6, seed = 11),
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 1001
	)

	# Cell 3: the band is forced to fail, the fit still has SEs.
	no_band <- with_mocked_bindings(
		suppressWarnings(betwfeWithSimulatedData(sim)),
		.apply_simultaneous_catt_band = function(x, alpha, has_valid_ses) NULL,
		.package = "fetwfe"
	)
	expect_true(no_band$calc_ses)
	expect_false(no_band$catt_band_applied)
	expect_true(no_band$p >= no_band$N * no_band$T)

	out_no_band <- suppressWarnings(.render_print(no_band))
	expect_true(any(grepl(
		"these are post-selection intervals",
		out_no_band,
		fixed = TRUE
	)))
	.expect_absent(out_no_band, "this band is the post-selection fallback")

	# Cell 2, unmocked: the band IS applied, so the caveat keeps its original
	# subject. Asserted so the two cells are DISCRIMINATED rather than one of
	# them merely being present.
	band <- suppressWarnings(betwfeWithSimulatedData(sim))
	expect_true(band$catt_band_applied)

	out_band <- suppressWarnings(.render_print(band))
	expect_true(any(grepl(
		"this band is the post-selection fallback",
		out_band,
		fixed = TRUE
	)))
	.expect_absent(out_band, "these are post-selection intervals")
})

# ------------------------------------------------------------------------------
# 8. (red) The event-study attribute survives truncation. A PIN, not a test of a
#    carry: `[.data.frame` copies the whole attribute list on ROW subsetting
#    (only the COLUMN branch rebuilds the frame), so `.truncate_event_study()`
#    needs no carry today. What this locks is that a future version which starts
#    doing column work cannot silently relabel every long panel as pointwise.
#
#    The fixture must be a long panel whose event-study band WAS applied, or the
#    assertion is false on the correct tree. Assert the VALUE, never equality
#    with the input: both sides move together under every mutation, so an
#    equality would be green on base, green patched, and green with the
#    attribute removed.
# ------------------------------------------------------------------------------
test_that("the event-study band_applied attribute survives truncation (#460)", {
	fit <- local({
		cf <- genCoefs(
			G = 3,
			T = 12,
			d = 2,
			density = 0.5,
			eff_size = 2,
			seed = 1
		)
		sim <- simulateData(
			cf,
			N = 200,
			sig_eps_sq = 1,
			sig_eps_c_sq = 0.5,
			seed = 1
		)
		fetwfe(
			pdata = sim$pdata,
			time_var = sim$time_var,
			unit_var = sim$unit_var,
			treatment = sim$treatment,
			response = sim$response,
			covs = sim$covs,
			q = 0.5,
			verbose = FALSE
		)
	})

	es <- eventStudy(fit)
	# More event times than `print()`'s `max_event_times` default of 10, so the
	# truncation below is a real one.
	expect_gt(nrow(es), 10L)
	expect_true(attr(es, "band_applied"))

	truncated <- fetwfe:::.truncate_event_study(es, 10L)
	expect_identical(nrow(truncated), 10L)
	expect_true(attr(truncated, "band_applied"))

	# ...and the rendered header on the truncated preview agrees.
	out <- .render_print(fit)
	.expect_header(out, .ES_HEADER, .SIMULTANEOUS)
})

# ------------------------------------------------------------------------------
# 9. (red) THE TWO HEADERS DISAGREE ON A PLAIN PUBLIC CALL, with no mocking.
#    This is what makes the event-study half of the change real rather than
#    mock-only, and it is the fixture whose print goldens move.
#
#    `generate_panel_data(N = 30, T = 5, R = 2, seed = 123)` with
#    `lambda_selection = "bic"` is the panel behind `_snaps/print-method-snapshot.md`
#    and `test-print-summary-single-source-439.R`'s `.fit_439()`. It is
#    degenerate: the bridge selects nothing, so `catt_df$se` is all zero -- the
#    COHORT band still applies (over an all-zero band, which is why
#    `catt_band_applied = TRUE` means "computed and written", not "informative")
#    while `.finish_event_study()`'s LOCAL `calc_ses` is FALSE and the
#    event-study band is never even requested.
# ------------------------------------------------------------------------------
test_that("the two families disagree on the standard print fixture (#460)", {
	fit <- fetwfe(
		pdata = generate_panel_data(N = 30, T = 5, R = 2, seed = 123),
		time_var = "time",
		unit_var = "unit",
		treatment = "treatment",
		response = "y",
		covs = c("cov1", "cov2"),
		verbose = FALSE,
		lambda_selection = "bic"
	)

	expect_identical(fit$ci_type, "simultaneous")
	expect_true(fit$catt_band_applied)
	expect_false(attr(eventStudy(fit), "band_applied"))

	out <- .render_print(fit)
	.expect_header(out, .CATT_HEADER, .SIMULTANEOUS)
	.expect_header(out, .ES_HEADER, .POINTWISE)

	out_s <- .render_summary(fit)
	.expect_header(out_s, .CATT_PREVIEW, .SIMULTANEOUS)
	.expect_header(out_s, .ES_PREVIEW, .POINTWISE)
})

# ------------------------------------------------------------------------------
# 10. (red under the gate's own reversion) C10 gates on the applied signal, not
#     on `ci_type`. This is the one hunk of #460 that no other assertion
#     observes: reverting `.check_ci_band_width()`'s gate to
#     `identical(x$ci_type, "simultaneous")` leaves the whole suite green,
#     measured. C10 is a wrong-number guardrail, so a gate that can un-narrow
#     itself silently is worth one block.
#
#     The block discriminates in BOTH directions, which is what makes it a pin
#     rather than a restatement: on a cell-3 object (band forced NULL, so
#     `catt_band_applied` FALSE while `ci_type` is still "simultaneous") the new
#     gate skips C10 and the old one would not; flipping ONLY the signal to TRUE
#     re-arms C10 on the identical `catt_df`, so the assertion cannot be
#     satisfied by a validator that has simply stopped checking.
# ------------------------------------------------------------------------------
test_that("C10 gates on the applied signal, not on ci_type (#460)", {
	fit <- with_mocked_bindings(
		.bas460_fp_fit(),
		.apply_simultaneous_catt_band = function(x, alpha, has_valid_ses) NULL,
		.package = "fetwfe"
	)
	expect_identical(fit$ci_type, "simultaneous")
	expect_false(fit$catt_band_applied)
	expect_true(all(is.finite(fit$catt_df$se) & fit$catt_df$se > 0))

	# Narrow every interval to half a pointwise width -- a C10 violation by
	# construction, on an object whose `ci_type` still reads "simultaneous".
	z <- stats::qnorm(1 - fit$alpha / 2)
	narrow <- fit
	narrow$catt_df$ci_low <- narrow$catt_df$estimate -
		0.5 * z * narrow$catt_df$se
	narrow$catt_df$ci_high <- narrow$catt_df$estimate +
		0.5 * z * narrow$catt_df$se

	# The gate this PR installs skips C10 here, because no band was applied.
	expect_silent(fetwfe:::.validate_fetwfe(narrow))

	# ...and C10 is still armed on the identical catt_df: flip only the signal.
	flipped <- narrow
	flipped$catt_band_applied <- TRUE
	expect_error(fetwfe:::.validate_fetwfe(flipped), "C10", fixed = TRUE)
})
