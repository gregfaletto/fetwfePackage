# Add S3 class FETWFE_tes with print/summary methods for getTes() output (issue #16)

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds. This document must be maintained in accordance with `.workflow/PLANS.md`.

## Purpose / Big Picture

Resolves [#16](https://github.com/gregfaletto/fetwfePackage/issues/16). After this change, the output of `getTes()` will be an S3 object of class `FETWFE_tes` instead of a plain `list`, with `print.FETWFE_tes` and `summary.FETWFE_tes` methods that render it nicely. The underlying data layout (`att_true`, `actual_cohort_tes`) is unchanged, so any code accessing fields by name continues to work; the only behavior change is that `getTes(coefs)` at the REPL now produces a self-describing display showing the overall true ATT, per-cohort effects labeled `Cohort 1` … `Cohort R`, and the generating parameters (`R`, `T`, `d`, `seed`) carried over from the source `FETWFE_coefs` object. `summary.FETWFE_tes` returns a structured list with the same content plus simple dispersion stats over the cohort effects (min, max, median, sd) — useful for quickly judging treatment-effect heterogeneity in a simulation; `print.summary.FETWFE_tes` renders that. New code lives in a new file `R/getTes_class.R` matching the `fetwfe_class.R` / `etwfe_class.R` convention; `getTes()` itself in `R/gen_funcs.R` is modified only to add a `class()` assignment to its return value. Existing tests in `tests/testthat/test-getTes.R` pass unchanged (the `expect_type(res, "list")` assertion still succeeds for a list-with-class). New tests cover that the class is set, that the print output contains the expected lines, and that summary returns the documented fields. Version bumps to 1.5.1 in `DESCRIPTION` / `NEWS.md` / `inst/CITATION`.

**Bundled scope (added in round 1 of plan review).** The PR also adds five lines to `.Rbuildignore` to silence two pre-existing CRAN NOTEs surfaced by the plan-review subagent's empirical `devtools::check()` runs against `origin/main` (one NOTE for hidden top-level dirs `.workflow` / `.claude` / `.plans`, one NOTE for top-level files `AGENTS.md` and `paper_arxiv.tex`). With those entries in `.Rbuildignore`, both NOTEs disappear from `R CMD check`. The change does not affect the installed package surface — `.Rbuildignore` only controls what is included in the source tarball — and it is bundled with #16 because the alternative is a redundant follow-up PR cycle for ~5 lines of edit.

The user-visible win, demonstrated:

    > coefs <- genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1, seed = 1)
    > getTes(coefs)
    True Treatment Effects (from FETWFE_coefs)
    ==========================================
    Overall true ATT: <number>

    Cohort effects:
      Cohort 1: <number>
      Cohort 2: <number>
      Cohort 3: <number>

    Generated from:
      Cohorts (R)      : 3
      Time periods (T) : 5
      Covariates (d)   : 2
      Seed             : 1

Before this change, the same call printed the raw list:

    > getTes(coefs)
    $att_true
    [1] <number>

    $actual_cohort_tes
    [1] <number> <number> <number>

## Progress

- [x] (2026-05-10) Issue clarification with Greg complete; scope locked.
- [x] (2026-05-10) Feature branch `feat/s3-class-getTes-16` created off `origin/main`.
- [x] (2026-05-10) `.plans/feat-s3-class-getTes-16/PLAN.md` drafted.
- [x] (2026-05-10) Plan-review subagent invoked per `.workflow/PLANS.md` § "Plan review and iteration"; feedback at `PLAN_feedback.md`.
- [x] (2026-05-10) Round-1 feedback applied to plan; response at `PLAN_feedback_response.md`. No round-2 needed (no substantive items). Bundled `.Rbuildignore` cleanup added to scope.
- [x] (2026-05-10) Greg reviewed the updated plan; sign-off received in chat.
- [x] (2026-05-10) Implementation milestone 1: wrote `R/getTes_class.R`, added `class()` assignment in `R/gen_funcs.R::getTes()`.
- [x] (2026-05-10) Implementation milestone 2: `devtools::document()` regenerated `man/getTes.Rd`, created `man/FETWFE_tes-class.Rd`, updated `NAMESPACE`. `\dontrun{}` example updated.
- [x] (2026-05-10) Implementation milestone 3: appended Tests 5/6/7 to `tests/testthat/test-getTes.R`.
- [x] (2026-05-10) Implementation milestone 4: version bumped 1.5.0 → 1.5.1 in `DESCRIPTION`, `NEWS.md`, `inst/CITATION`.
- [x] (2026-05-10) Implementation milestone 5: added 5 lines to `.Rbuildignore`. Bootstrap of `inst/WORDLIST` (36 entries) added in same commit.
- [x] (2026-05-10) Per-PR CRAN gate clean: `air format .` no-op; `devtools::document()` idempotent; `devtools::check(args = "--as-cran")` `0/0/0`; `devtools::test()` 1065/1065 PASS; `devtools::spell_check()` empty; `urlchecker::url_check()` only the 2 pre-existing redirects in `vignettes/vignette.Rmd`. (`build_vignettes()` not required; vignettes untouched.)
- [x] (2026-05-10) Implementation committed as `3688aff`.
- [x] (2026-05-10) Post-execution review subagent invoked per `.workflow/REVIEW_AGENT.md`; feedback at `PLAN_implementation_feedback.md`. Verdict: LGTM, no blockers, no action items.
- [ ] Push branch to `origin/feat/s3-class-getTes-16`, open PR targeting `main`, link issue #16. (Awaiting Greg's final checkpoint.)

## Surprises & Discoveries

- Observation: `expect_type(structure(list(), class = "FETWFE_tes"), "list")` returns TRUE.
  Evidence: `typeof()` is unaffected by S3 class assignment, and testthat's `expect_type` uses `typeof()`. So existing Test 1 (`expect_type(res, "list")`) continues to pass after `class()` is assigned. Empirically verified by the plan-review subagent.

- Observation: `coefs_obj$seed` is `NULL`, not `NA`, when `genCoefs()` is called without a `seed` argument.
  Evidence: Plan-review subagent ran `coefs <- genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1)` (no seed) and verified `is.null(coefs$seed)` returns TRUE. So the planned `is.null(x$seed)` guard in `print.FETWFE_tes` is correct. `seed = 0` and an integer seed both render correctly.

- Observation: `c(min = 1, max = 2, median = 3, sd = 4)` preserves declaration order, and `expect_named(..., c("min", "max", "median", "sd"))` is order-sensitive.
  Evidence: Plan-review subagent verified empirically that reordering the names breaks the assertion. So the planned construction in `summary.FETWFE_tes` matches what Test 7 asserts.

- Observation: `devtools::check()` on `origin/main` produces 2 pre-existing NOTEs unrelated to issue #16.
  Evidence: One NOTE for hidden top-level directories `.workflow` / `.claude` / `.plans`; one NOTE for top-level files `AGENTS.md` and `paper_arxiv.tex`. None are listed in `.Rbuildignore`. Bundled into this PR's scope (see Purpose / Big Picture, second paragraph).

- Observation: `R = 2` is a real edge case for `stats::sd` of the cohort vector but does not trigger the planned `length(tes) >= 2` guard, since `genCoefs()` already enforces `R >= 2` at line 242–246 of `R/gen_funcs.R`.
  Evidence: Plan-review subagent ran `genCoefs(R = 2, ...)` and confirmed `length(actual_cohort_tes) == 2` and `stats::sd(tes)` returned a real number. The guard is retained as defense-in-depth in case `genCoefs`'s validation ever changes (Greg's directive).

## Decision Log

- Decision: Class name is `FETWFE_tes`.
  Rationale: Matches the existing `FETWFE_coefs` / `FETWFE_simulated` family naming for non-result classes (capital prefix, snake_case suffix). Result classes use lowercase (`fetwfe`, `etwfe`) but the output of `getTes()` is *derived from* a `FETWFE_coefs` object, not an estimator result, so it belongs in the same family as its input.
  Date/Author: 2026-05-10 / Greg + Claude (clarification dialogue).

- Decision: `print.FETWFE_tes` shows medium detail — overall ATT, labeled cohort effects, plus generating parameters (R, T, d, seed) from the source coefs object.
  Rationale: Greg picked the medium option from three offered (minimal / medium / like-print.fetwfe). Self-describing without being overwrought for two pieces of data.
  Date/Author: 2026-05-10 / Greg.

- Decision: Add a `summary.FETWFE_tes` method too, with dispersion stats (min, max, median, sd of cohort effects) as the differentiator from print.
  Rationale: Greg confirmed both methods. Dispersion stats give a quick read on treatment-effect heterogeneity in simulation without forcing the user to compute them by hand.
  Date/Author: 2026-05-10 / Greg.

- Decision: New file `R/getTes_class.R`; do not fold S3 methods into `R/gen_funcs.R`.
  Rationale: Matches the existing `fetwfe_class.R` / `etwfe_class.R` convention; keeps `gen_funcs.R` focused on data-generation utilities.
  Date/Author: 2026-05-10 / Claude.

- Decision: Patch version bump (1.5.0 → 1.5.1), not minor.
  Rationale: No existing API breaks; `att_true` and `actual_cohort_tes` field names unchanged; `expect_type(res, "list")` test continues to pass. The new methods are an addition without behavior change to existing fields. Greg confirmed.
  Date/Author: 2026-05-10 / Greg.

- Decision: Cohort effects in `actual_cohort_tes` remain an unnamed numeric vector at the data layer; cohort labeling (`Cohort 1`, …, `Cohort R`) happens only in the print methods.
  Rationale: Adding names to the vector would be a (small) data-layout change that downstream code accessing by index could in principle notice. Keeping the data layer untouched and labeling only at the presentation layer is the safer, smaller change. No `R/` or `vignettes/` callers were found that depend on either layout.
  Date/Author: 2026-05-10 / Claude.

- Decision: Keep the `length(tes) >= 2` guard in `summary.FETWFE_tes` despite it being unreachable through the current public API.
  Rationale: Greg directed (round 1 feedback response) that the guard stays as defense-in-depth in case `genCoefs()`'s `R >= 2` validation ever changes upstream. A code comment in the new file marks it as defensive. The user-facing string in `print.summary.FETWFE_tes` is changed from `"<NA (R<2)>"` to `"<NA>"` — the `(R<2)` annotation would mislead users into thinking R<2 is reachable.
  Date/Author: 2026-05-10 / Greg.

- Decision: Bundle a small `.Rbuildignore` cleanup with this PR.
  Rationale: Plan-review subagent's `devtools::check()` on `origin/main` surfaced 2 pre-existing CRAN NOTEs (hidden dirs `.workflow` / `.claude` / `.plans`, top-level files `AGENTS.md` and `paper_arxiv.tex`). Adding 5 regex lines to `.Rbuildignore` silences both NOTEs. Splitting into a separate follow-up PR cycle would double workflow overhead for ~5 lines of edit. Greg directed the bundle in chat after seeing the feedback. PR description will call out the bundled scope.
  Date/Author: 2026-05-10 / Greg.

- Decision: Empirically verified `expect_type(res, "list")` continues to pass with class assignment.
  Rationale: Plan-review subagent verified that `typeof(structure(list(), class = "FETWFE_tes"))` returns `"list"`, and testthat's `expect_type` uses `typeof()`. So existing Test 1 in `tests/testthat/test-getTes.R` is not affected by the class change. Recorded here per round-1 feedback action item 4.
  Date/Author: 2026-05-10 / plan-review subagent.

- Decision: Adopt the unconditional per-PR CRAN-readiness gate (`--as-cran` + `spell_check` + `url_check`, plus `build_vignettes` when vignettes change).
  Rationale: Greg directed that every PR cycle leave the package CRAN-ready, not just CRAN-bound ones. Workflow docs updated; this plan's Step 5 and Validation criteria expanded to match. The marginal cost is small (`spell_check` and `url_check` are seconds each) and catches typos / dead URLs that would otherwise accumulate as deferred CRAN-NOTE backlog.
  Date/Author: 2026-05-10 / Greg.

## Outcomes & Retrospective

**What was built.** Issue #16 resolved: `getTes()` output is now an S3 object of class `FETWFE_tes`, with `print()` / `summary()` / `print.summary()` methods. Existing field accessors (`$att_true`, `$actual_cohort_tes`) unchanged; existing tests pass without modification. Bundled scope: 5-line `.Rbuildignore` cleanup (silenced 2 pre-existing CRAN NOTEs), 36-term `inst/WORDLIST` bootstrap (`spell_check()` now passes cleanly), version bump 1.5.0 → 1.5.1 across DESCRIPTION/NEWS.md/inst/CITATION. Post-execution review (subagent) verdict: LGTM, no blockers, no action items.

**Per-PR CRAN gate established.** This is the first PR cycle to run the full per-PR CRAN gate (`devtools::check(args = "--as-cran")` + `devtools::spell_check()` + `urlchecker::url_check()`). The gate is now clean; future PRs inherit a CRAN-ready baseline and need only maintain it.

**`.plans/<branch>/` tracking convention adopted.** Greg established mid-cycle that `.plans/<branch>/` folders are tracked in version control alongside the feature branch, so the design dialogue (plan, subagent feedback, response, post-execution feedback) rides into `main` as part of the PR. Workflow docs were updated to reflect this.

**Lessons learned.**
- The plan-review subagent's empirical "overlay + load_all + document + test + check" methodology caught zero plan-shape bugs but contributed five Surprises & Discoveries entries that materially de-risked implementation. Cost: ~14 minutes of subagent runtime. Worth it.
- The post-execution review's "build the actual tarball and grep its contents" step verified the `.Rbuildignore` cleanup definitively (no `.workflow`, `.claude`, `.plans`, `AGENTS.md`, `paper_arxiv.tex` in the built `.tar.gz`). This is a useful pattern to retain for any PR that touches `.Rbuildignore`.
- The `air format .` doc warning ("preserve tabs by hand because no `air.toml` exists") was empirically wrong — air on this repo is a no-op (auto-detects existing tab indentation somehow). Workflow doc fixed during this cycle.
- Pre-existing `urlchecker::url_check()` redirects (`vignettes/vignette.Rmd:29`, `:31`, both bookdown.org → posit.cloud) and 3 pre-existing `@details`/`@seealso` roxygen warnings in `fetwfe_core.R` (`getPsiRFused`, `getSecondVarTermDataApp`, `getTeResults2`) are flagged here as candidate follow-ups, not in scope for #16.

**Diff size.** 12 files, +310/-38 lines net. Well within the "split if >150 lines" heuristic when measured against source files only (~95 source lines added; the rest is regenerated `man/*.Rd`, `NAMESPACE`, `inst/WORDLIST`, plan updates, version-string updates).

**Match to clarified scope.** All four bullets of the clarified-scope paragraph (the dialog with Greg in `CHOOSING_A_PR.md` § "Issue clarification") realized in code: class set, print method medium-detail, summary with dispersion stats, version bump as patch.

## Context and Orientation

`fetwfe` is an R package implementing fused extended two-way fixed effects (see `paper_arxiv.tex`). The functions `genCoefs()` and `simulateData()` exist for simulation studies: `genCoefs()` generates a coefficient vector with known true treatment effects (returning an object of class `FETWFE_coefs`), `simulateData()` simulates panel data from those coefficients. `getTes()` extracts the true treatment effects from a `FETWFE_coefs` object — the overall true ATT and the per-cohort true ATTs — so a simulation study can compare estimator output against the truth.

Currently `getTes()` returns a plain `list`. Other public objects in the package have S3 classes and dedicated print methods: `fetwfe()` returns class `"fetwfe"` with methods in `R/fetwfe_class.R`; `etwfe()` returns class `"etwfe"` with methods in `R/etwfe_class.R`; `genCoefs()` returns class `"FETWFE_coefs"`; `simulateData()` returns class `"FETWFE_simulated"` (neither of the latter two has a print method, but they carry the class for type-checking via `inherits()`).

Issue #16 asks for the same treatment for `getTes()`: a class plus a nice print method (Greg confirmed in chat that summary is also wanted).

Files involved:

- `R/gen_funcs.R` — defines `getTes()` (lines ~305–379) and the internal helper `getActualCohortTes()` (lines ~1365–1380). Also defines `genCoefs()` which produces the input `FETWFE_coefs` object.
- `R/fetwfe_class.R` — precedent for how S3 classes are structured in this package. Pattern: `coef.<class>`, `print.<class>`, `summary.<class>`, `print.summary.<class>` methods, plus internal helpers like `.truncate_catt` and `.print_catt_tbl`.
- `R/etwfe_class.R` — same pattern for `etwfe()`.
- `tests/testthat/test-getTes.R` — existing four tests for `getTes()`. Two of these check shape (`expect_type(res, "list")`, presence of `att_true` / `actual_cohort_tes` names); the other two check numerical correctness. None check the class — all four continue to pass after the class is added.
- `man/getTes.Rd` — generated by roxygen; will be regenerated when `getTes()`'s `@return` is updated.
- `NAMESPACE` — generated by roxygen; will gain `S3method(print, FETWFE_tes)`, `S3method(summary, FETWFE_tes)`, `S3method(print, summary.FETWFE_tes)` after `devtools::document()`.
- `DESCRIPTION` — `Version: 1.5.0` → `Version: 1.5.1`.
- `NEWS.md` — new section header `## Version 1.5.1 (2026-05-10)` with one bullet.
- `inst/CITATION` — `note = "R package version 1.5.0"` → `note = "R package version 1.5.1"`.

The current `getTes()` body (for context — implementation should not change other than the final `class()` assignment):

    getTes <- function(coefs_obj) {
        if (!inherits(coefs_obj, "FETWFE_coefs")) {
            stop("coefs_obj must be an object of class 'FETWFE_coefs'")
        }
        beta <- coefs_obj$beta
        R <- coefs_obj$R
        T <- coefs_obj$T
        d <- coefs_obj$d
        num_treats <- getNumTreats(R = R, T = T)
        p <- getP(R = R, T = T, d = d, num_treats = num_treats)
        stopifnot(length(beta) == p)
        first_inds <- getFirstInds(R = R, T = T)
        treat_inds <- getTreatInds(R = R, T = T, d = d, num_treats = num_treats)
        actual_cohort_tes <- getActualCohortTes(
            R = R,
            first_inds = first_inds,
            treat_inds = treat_inds,
            coefs = beta,
            num_treats = num_treats
        )
        att_true <- as.numeric(mean(actual_cohort_tes))
        return(list(att_true = att_true, actual_cohort_tes = actual_cohort_tes))
    }

The `coefs_obj` carries `R`, `T`, `d`, and `seed`. The new `getTes()` will additionally carry these forward into the returned object so the print method has self-describing context without needing to be passed `coefs_obj` again.

## Plan of Work

The work breaks into four logical edits, ordered to keep `devtools::test()` passing after each.

**1. Modify `getTes()` in `R/gen_funcs.R` to return a classed object.** Change the final `return(list(...))` to assemble a list that also carries `R`, `T`, `d`, `seed` (copied from `coefs_obj`), and assign `class(out) <- "FETWFE_tes"`. The two original fields (`att_true`, `actual_cohort_tes`) keep their names and types. Update the function's roxygen `@return` block to document the class and the four new fields.

After this edit alone, the existing tests still pass (a list-with-class still satisfies `expect_type(res, "list")`), but `print(getTes(coefs))` prints the raw list because no methods are registered yet.

**2. Add `R/getTes_class.R` with the S3 methods.** New file mirroring the structure of `R/fetwfe_class.R`:

- A roxygen block defining the class doc (`@title Compute True Treatment Effects Output Class`, `@name FETWFE_tes-class`, `NULL`).
- `print.FETWFE_tes(x, ...)` with `#' @export`, producing the medium-detail output shown in the Purpose section.
- `summary.FETWFE_tes(object, ...)` with `#' @export`, returning a `summary.FETWFE_tes`-classed list with fields: `att_true`, `actual_cohort_tes`, `R`, `T`, `d`, `seed`, plus `cohort_te_stats` (a named numeric vector with elements `min`, `max`, `median`, `sd`).
- `print.summary.FETWFE_tes(x, ...)` with `#' @export`, rendering the summary list with the dispersion stats added below the cohort effects block.

All three methods return `invisible(x)` (or `invisible(object)` for `summary`), matching `print.fetwfe`.

**3. Add tests to `tests/testthat/test-getTes.R`.** Append three new `test_that` blocks (Test 5–7):

- Test 5: `getTes()` returns an object of class `FETWFE_tes` and the new `R`, `T`, `d`, `seed` fields are populated correctly.
- Test 6: `print(getTes(coefs))` writes lines containing "Overall true ATT", "Cohort 1", "Cohort 2", "Cohort 3", "Cohorts (R)", "Time periods (T)", "Covariates (d)", "Seed". Use `capture.output()` and grep with `expect_match`.
- Test 7: `summary(getTes(coefs))` returns a `summary.FETWFE_tes` list whose `cohort_te_stats` element is a length-4 numeric vector named `c("min", "max", "median", "sd")`, and whose values agree with `min/max/median/sd` of the cohort vector. Also verify that `print(summary(getTes(coefs)))` includes a dispersion-stats block (grep for "Cohort effect dispersion" or whatever header label is used).

The four existing tests are not modified.

**4. Bump version + NEWS + CITATION.** Edit `DESCRIPTION` `Version:` line. Add a new section to the top of `NEWS.md`:

    ## Version 1.5.1 (2026-05-10)

    - Output of `getTes()` is now an S3 object of class `FETWFE_tes` with `print()` and `summary()` methods. Existing field accessors (`$att_true`, `$actual_cohort_tes`) are unchanged.

Edit `inst/CITATION` to bump the `note = "R package version 1.5.0"` line to `1.5.1`.

## Concrete Steps

All commands run from the repo root `/Users/gregfaletto/Documents/GitHub/fetwfePackage`.

**Step 1 — Edit `R/gen_funcs.R` `getTes()`.** Replace the body's final `return(list(...))` with:

    out <- list(
        att_true = att_true,
        actual_cohort_tes = actual_cohort_tes,
        R = R,
        T = T,
        d = d,
        seed = coefs_obj$seed
    )
    class(out) <- "FETWFE_tes"
    return(out)

Update the function's `@return` roxygen block to read approximately:

    #' @return An object of class \code{"FETWFE_tes"}, which is a list with the
    #' following elements:
    #' \describe{
    #'   \item{att_true}{Numeric scalar; the overall true ATT (equal-weighted
    #'         mean of cohort-specific effects).}
    #'   \item{actual_cohort_tes}{Numeric vector of length \code{R}; the true
    #'         per-cohort treatment effects.}
    #'   \item{R, T, d, seed}{The generating parameters carried over from
    #'         \code{coefs_obj} so the print method can be self-describing.}
    #' }
    #' Use \code{print()} or \code{summary()} on the returned object for a
    #' formatted display.

Update the existing `\dontrun{}` `@examples` block to add a final line demonstrating the new print method:

    #' @examples
    #' \dontrun{
    #' # Generate coefficients
    #' coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)
    #'
    #' # Compute the true treatment effects:
    #' te_results <- getTes(coefs)
    #'
    #' # Overall average treatment effect on the treated:
    #' print(te_results$att_true)
    #'
    #' # Cohort-specific treatment effects:
    #' print(te_results$actual_cohort_tes)
    #'
    #' # Or use the new print method for a self-describing display:
    #' print(te_results)
    #' }

**Step 2 — Create `R/getTes_class.R`.** Full file structure:

    #' @title Compute True Treatment Effects Output Class
    #' @description S3 class for the output of \code{getTes()}.
    #' @name FETWFE_tes-class
    NULL

    #' @export
    print.FETWFE_tes <- function(x, ...) {
        cat(
            "True Treatment Effects (from FETWFE_coefs)\n",
            "==========================================\n",
            sep = ""
        )
        cat(sprintf("Overall true ATT: %.4f\n\n", x$att_true))
        cat("Cohort effects:\n")
        for (r in seq_along(x$actual_cohort_tes)) {
            cat(sprintf("  Cohort %d: %.4f\n", r, x$actual_cohort_tes[r]))
        }
        cat("\n")
        cat("Generated from:\n")
        cat(sprintf("  Cohorts (R)      : %d\n", x$R))
        cat(sprintf("  Time periods (T) : %d\n", x$T))
        cat(sprintf("  Covariates (d)   : %d\n", x$d))
        cat(sprintf("  Seed             : %s\n",
                    if (is.null(x$seed)) "<none>" else as.character(x$seed)))
        invisible(x)
    }

    #' @export
    summary.FETWFE_tes <- function(object, ...) {
        tes <- object$actual_cohort_tes
        # Defensive: genCoefs() currently enforces R >= 2, so length(tes) >= 2
        # is always true in practice. Guard is retained in case that validation
        # ever changes upstream (sd() of a length-1 vector is NA).
        out <- list(
            att_true = object$att_true,
            actual_cohort_tes = tes,
            R = object$R,
            T = object$T,
            d = object$d,
            seed = object$seed,
            cohort_te_stats = c(
                min    = min(tes),
                max    = max(tes),
                median = stats::median(tes),
                sd     = if (length(tes) >= 2) stats::sd(tes) else NA_real_
            )
        )
        structure(out, class = "summary.FETWFE_tes")
    }

    #' @export
    print.summary.FETWFE_tes <- function(x, ...) {
        cat(
            "Summary of True Treatment Effects (from FETWFE_coefs)\n",
            "=====================================================\n",
            sep = ""
        )
        cat(sprintf("Overall true ATT: %.4f\n\n", x$att_true))
        cat("Cohort effects:\n")
        for (r in seq_along(x$actual_cohort_tes)) {
            cat(sprintf("  Cohort %d: %.4f\n", r, x$actual_cohort_tes[r]))
        }
        cat("\n")
        cat("Cohort effect dispersion:\n")
        cat(sprintf("  min    : %.4f\n", x$cohort_te_stats["min"]))
        cat(sprintf("  max    : %.4f\n", x$cohort_te_stats["max"]))
        cat(sprintf("  median : %.4f\n", x$cohort_te_stats["median"]))
        cat(sprintf("  sd     : %s\n",
                    if (is.na(x$cohort_te_stats["sd"])) "<NA>"
                    else sprintf("%.4f", x$cohort_te_stats["sd"])))
        cat("\n")
        cat("Generated from:\n")
        cat(sprintf("  Cohorts (R)      : %d\n", x$R))
        cat(sprintf("  Time periods (T) : %d\n", x$T))
        cat(sprintf("  Covariates (d)   : %d\n", x$d))
        cat(sprintf("  Seed             : %s\n",
                    if (is.null(x$seed)) "<none>" else as.character(x$seed)))
        invisible(x)
    }

Note: the file uses tab indentation per repo convention. The above is shown with spaces for plan readability; convert to tabs before saving. `air format .` will normalize.

Note: `stats::median` and `stats::sd` are not currently in the `@importFrom stats` list at the top of `R/fetwfe.R` — `sd` is already imported (line 38 of NAMESPACE: `importFrom(stats,sd)`), but `median` is not. Either add `@importFrom stats median` to the new file's roxygen block OR call `stats::median(tes)` qualified. Plan uses qualified calls to avoid touching the global imports list.

**Step 3 — Edit `tests/testthat/test-getTes.R`.** Append after line 132:

    # ----------------------------------------------------------------------
    # Test 5: Class + carried-over parameters.
    # ----------------------------------------------------------------------
    test_that("getTes returns an object of class FETWFE_tes with carried params", {
        coefs <- genCoefs(R = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
                          seed = 234)
        res <- getTes(coefs)
        expect_s3_class(res, "FETWFE_tes")
        expect_equal(res$R, 5L)
        expect_equal(res$T, 30L)
        expect_equal(res$d, 12L)
        expect_equal(res$seed, 234)
    })

    # ----------------------------------------------------------------------
    # Test 6: print method writes the expected lines.
    # NOTE: the regexes below intentionally lock the print-output labels
    # ("Cohorts (R)", "Time periods (T)", "Covariates (d)"). Any cosmetic
    # relabeling of print.FETWFE_tes requires updating this test.
    # ----------------------------------------------------------------------
    test_that("print.FETWFE_tes writes the expected sections", {
        coefs <- genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1,
                          seed = 1)
        out <- capture.output(print(getTes(coefs)))
        joined <- paste(out, collapse = "\n")
        expect_match(joined, "Overall true ATT")
        expect_match(joined, "Cohort 1")
        expect_match(joined, "Cohort 2")
        expect_match(joined, "Cohort 3")
        expect_match(joined, "Cohorts \\(R\\)")
        expect_match(joined, "Time periods \\(T\\)")
        expect_match(joined, "Covariates \\(d\\)")
        expect_match(joined, "Seed")
    })

    # ----------------------------------------------------------------------
    # Test 7: summary method returns the documented fields and dispersion.
    # ----------------------------------------------------------------------
    test_that("summary.FETWFE_tes returns expected fields and dispersion stats", {
        coefs <- genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1,
                          seed = 1)
        s <- summary(getTes(coefs))
        expect_s3_class(s, "summary.FETWFE_tes")
        expect_named(s$cohort_te_stats, c("min", "max", "median", "sd"))
        expect_equal(unname(s$cohort_te_stats["min"]),
                     min(s$actual_cohort_tes))
        expect_equal(unname(s$cohort_te_stats["max"]),
                     max(s$actual_cohort_tes))
        expect_equal(unname(s$cohort_te_stats["median"]),
                     stats::median(s$actual_cohort_tes))
        expect_equal(unname(s$cohort_te_stats["sd"]),
                     stats::sd(s$actual_cohort_tes))

        out <- capture.output(print(s))
        joined <- paste(out, collapse = "\n")
        expect_match(joined, "Cohort effect dispersion")
    })

**Step 4 — Bump version files.**

`DESCRIPTION`: change `Version: 1.5.0` → `Version: 1.5.1`.

`NEWS.md`: insert at the top (above the existing `## Version 1.5.0` header):

    # NEWS

    ## Version 1.5.1 (2026-05-10)

    - Output of `getTes()` is now an S3 object of class `FETWFE_tes` with
      `print()` and `summary()` methods. Existing field accessors
      (`$att_true`, `$actual_cohort_tes`) are unchanged.

(The existing `# NEWS` header at the top of `NEWS.md` should be preserved exactly once; insert the new version section between it and the existing 1.5.0 section.)

`inst/CITATION`: change `note = "R package version 1.5.0"` → `note = "R package version 1.5.1"`.

`cran-comments.md` is intentionally **not** bumped. That file tracks CRAN-submission state (last submission notes, reviewer feedback), not the local working version, and is updated separately at submission time.

**Step 4b — Add five lines to `.Rbuildignore`.** Append to the end of `.Rbuildignore`:

    ^\.workflow$
    ^\.claude$
    ^\.plans$
    ^AGENTS\.md$
    ^paper_arxiv\.tex$

These silence the two pre-existing `R CMD check` NOTEs surfaced during plan review (one for hidden top-level dirs, one for the two top-level files). The `^\.plans$` entry is particularly important now that `.plans/<branch>/` is **tracked** in version control (per `.workflow/WORKFLOW_LESSONS.md` §9 — the design dialogue rides along with the feature branch's PR) but should still not enter the CRAN tarball. The other four paths are personal-tooling or methodology-source artifacts that do not belong in a CRAN tarball either. The package's installed surface is unaffected by any of the five entries.

After this step, `devtools::check()` should produce `0 errors / 0 warnings / 0 notes` rather than `0 errors / 0 warnings / 2 notes`.

**Step 5 — Format, regenerate docs, run the full per-PR CRAN gate.** From the repo root:

    air format .

Then in an R session:

    devtools::document()
    devtools::test()
    devtools::check(args = "--as-cran")
    devtools::spell_check()
    urlchecker::url_check()
    # vignettes/*.Rmd is not in this PR's diff, so build_vignettes() is not required.

Expected:

- `devtools::document()` regenerates `man/getTes.Rd` (updated `@return` and `\dontrun{}` example), creates `man/FETWFE_tes-class.Rd`, and updates `NAMESPACE` to add three `S3method(...)` lines. A second run produces no diff (idempotent).
- `devtools::test()` reports all tests passing — the four existing `getTes` tests plus the three new ones, plus everything else in the suite.
- `devtools::check(args = "--as-cran")` exits with `0 errors ✔ | 0 warnings ✔ | 0 notes ✔` (assumes Step 4b's `.Rbuildignore` cleanup landed).
- `devtools::spell_check()` produces empty output — *or*, on first run for this repo, flags domain terms (`FETWFE`, `ETWFE`, `BETWFE`, `ATT`, `CATT`, `Faletto`, `Wooldridge`, `Pesaran`, etc.) which need to be added to `inst/WORDLIST` (one term per line). If a WORDLIST update is needed, it lands in the same commit as the rest of this PR. **Note:** `inst/WORDLIST` may not exist on `origin/main`; if so, this PR creates it. Real typos in any new strings must be fixed (the new strings are: "True Treatment Effects (from FETWFE_coefs)", "Overall true ATT", "Cohort effects", "Generated from", "Cohorts (R)", "Time periods (T)", "Covariates (d)", "Seed", "Summary of True Treatment Effects (from FETWFE_coefs)", "Cohort effect dispersion", "min", "max", "median", "sd").
- `urlchecker::url_check()` reports all URLs reachable. This PR introduces no new URLs, so the result is whatever the baseline is — flag any pre-existing dead URLs in the PR description but don't block merge on them (separate concern).
- `devtools::build_vignettes()` is **not** in the gate for this PR because `vignettes/*.Rmd` is not in the diff. (Per `.workflow/PLANS.md` § "Standard acceptance criteria" #4.)

If `air format .` modifies any file, stage and commit alongside the source edit. If `inst/WORDLIST` is created or updated, add it to the commit's `git add` list (see Step 6).

**Step 6 — Commit.** Stage exactly:

- `R/gen_funcs.R` (modified `getTes()` body and roxygen)
- `R/getTes_class.R` (new file)
- `tests/testthat/test-getTes.R` (appended tests)
- `DESCRIPTION` (version bump)
- `NEWS.md` (new section)
- `inst/CITATION` (version bump)
- `.Rbuildignore` (5 new lines)
- `inst/WORDLIST` (created or updated, only if `spell_check()` flagged legitimate terms — see Step 5)
- `man/getTes.Rd` (regenerated)
- `man/FETWFE_tes-class.Rd` (new, regenerated)
- `NAMESPACE` (regenerated)

Single commit (omit `inst/WORDLIST` from the `git add` list if Step 5's `spell_check()` was empty):

    git add R/gen_funcs.R R/getTes_class.R tests/testthat/test-getTes.R \
            DESCRIPTION NEWS.md inst/CITATION .Rbuildignore \
            inst/WORDLIST \
            man/getTes.Rd man/FETWFE_tes-class.Rd NAMESPACE
    git commit -m "Add S3 class FETWFE_tes for getTes() output (#16)"

Verify with `git status` — only `.claude/` (and any other personal-state directories) should remain untracked. The `.plans/feat-s3-class-getTes-16/` folder is tracked (per `.workflow/WORKFLOW_LESSONS.md` §9) but the plan/feedback/response files were already committed in `8fdba98 update plan` and should not show as modified unless this PR's work caused further plan edits.

## Validation and Acceptance

Per `.workflow/PLANS.md` § "Standard acceptance criteria" — the unconditional per-PR CRAN gate:

1. `devtools::check(args = "--as-cran")` exits with `0 errors ✔ | 0 warnings ✔ | 0 notes ✔`. (Achievable with Step 4b's `.Rbuildignore` cleanup landed.)
2. `devtools::spell_check()` produces empty output. If domain terms were flagged on first run, they have been added to `inst/WORDLIST` in the same commit and a re-run is now empty.
3. `urlchecker::url_check()` reports all URLs reachable. (No new URLs introduced by this PR; pre-existing URL-check state is the baseline.)
4. `devtools::build_vignettes()` is **not required** for this PR: `vignettes/*.Rmd` is not in the diff.
5. `devtools::test()` reports all tests passing. The three new tests (5, 6, 7) fail before the source changes and pass after — verify by running on `origin/main` first, then on this branch.
6. `devtools::document()` is idempotent: a second run produces no diff.
7. `NEWS.md` has a new bullet under a new `## Version 1.5.1` header. `DESCRIPTION` `Version:` is bumped to 1.5.1.
8. `inst/CITATION` `note` line references `1.5.1`.
9. Every named declaration in the "Interfaces and Dependencies" section appears in `git diff origin/main` with matching signature; no extra unplanned exports appear.
10. The `Surprises & Discoveries / Decision Log / Outcomes & Retrospective` sections of this plan are kept current.
11. No math change; no `paper_arxiv.tex` citation needed.

PR-specific acceptance:

- `getTes(genCoefs(...))` returns an object of class `FETWFE_tes`. `inherits(out, "FETWFE_tes")` is `TRUE`.
- `print(getTes(coefs))` produces the medium-detail output shown in the Purpose section.
- `summary(getTes(coefs))` produces a `summary.FETWFE_tes`-classed list with fields `att_true`, `actual_cohort_tes`, `R`, `T`, `d`, `seed`, `cohort_te_stats`. The `cohort_te_stats` element is a length-4 numeric vector named `c("min", "max", "median", "sd")`.
- Existing field access `getTes(coefs)$att_true` and `getTes(coefs)$actual_cohort_tes` continues to work, returning the same numeric values as on `origin/main`.

## Idempotence and Recovery

All steps are idempotent. `devtools::document()`, `devtools::test()`, and `devtools::check()` are safe to re-run; `air format .` is idempotent on already-formatted code. If a step fails:

- If `devtools::document()` produces unexpected diffs (e.g., because the executor's roxygen edit had a syntax error), revert the offending source edit, re-run `document()`, fix, re-run.
- If `devtools::check()` reports a new NOTE, identify the cause from the check log, fix, re-run.
- If a test fails, the implementation is wrong — fix the implementation, not the test.
- If the branch becomes tangled, `git reset --hard origin/main` and start over from the planned edits. The plan is the source of truth; the working tree is recoverable.

## Artifacts and Notes

Expected new file `man/FETWFE_tes-class.Rd` (small, generated):

    \name{FETWFE_tes-class}
    \alias{FETWFE_tes-class}
    \title{Compute True Treatment Effects Output Class}
    \description{S3 class for the output of \code{getTes()}.}

Expected new lines in `NAMESPACE` (after `devtools::document()`):

    S3method(print, FETWFE_tes)
    S3method(print, summary.FETWFE_tes)
    S3method(summary, FETWFE_tes)

Expected `git diff` summary against `origin/main` (counting source only, not regenerated docs):

- `R/gen_funcs.R`: ~15 lines changed (new return-list assembly, new roxygen `@return` block, one-line addition to the `\dontrun{}` example).
- `R/getTes_class.R`: ~78 new lines (with the defensive-guard comment).
- `tests/testthat/test-getTes.R`: ~65 new lines (with the regex-locking note above Test 6).
- `DESCRIPTION`: 1 line.
- `NEWS.md`: 5–6 new lines.
- `inst/CITATION`: 1 line.
- `.Rbuildignore`: 5 new lines.

Total source change: ~160 lines, near but acceptable around the "split if exceeds ~150 lines" heuristic from PLANS.md § "PR scope guidance" (the bundled `.Rbuildignore` cleanup is small and tightly related to CRAN-cleanliness of the same PR).

## Interfaces and Dependencies

In `R/getTes_class.R`, define:

    print.FETWFE_tes         <- function(x, ...) ...
    summary.FETWFE_tes       <- function(object, ...) ...
    print.summary.FETWFE_tes <- function(x, ...) ...

All three carry `#' @export`. The class is documented via:

    #' @title Compute True Treatment Effects Output Class
    #' @description S3 class for the output of \code{getTes()}.
    #' @name FETWFE_tes-class
    NULL

In `R/gen_funcs.R`, modify `getTes()` to:

    getTes <- function(coefs_obj) {
        # ... existing body unchanged through actual_cohort_tes / att_true ...
        out <- list(
            att_true = att_true,
            actual_cohort_tes = actual_cohort_tes,
            R = R,
            T = T,
            d = d,
            seed = coefs_obj$seed
        )
        class(out) <- "FETWFE_tes"
        return(out)
    }

The `@return` roxygen block on `getTes()` is updated as shown in Step 1.

No new package dependencies. `stats::median` and `stats::sd` are called via fully-qualified names (`sd` is already in `@importFrom stats`; `median` is added qualified to avoid a global imports change).

No changes to `Imports:` or `Suggests:` in `DESCRIPTION`.

---

## Revision history

- **2026-05-10 (initial draft)** — Plan drafted from issue-clarification dialogue with Greg. Class name `FETWFE_tes`, medium-detail print, summary with min/max/median/sd dispersion stats, patch version bump 1.5.0 → 1.5.1, all decisions confirmed in chat.

- **2026-05-10 (round 1 plan-review applied)** — Plan-review subagent ran empirical validation (overlay + load_all + document + test + check); no blockers, six minor action items. Greg directed: keep the `length(tes) >= 2` guard for defensive resiliency (with a comment) but change `"<NA (R<2)>"` → `"<NA>"`; bundle a 5-line `.Rbuildignore` cleanup to silence the 2 pre-existing CRAN NOTEs the subagent surfaced; apply remaining items 2–6 (note `cran-comments.md` not bumped, code-comment Test 6 regex sensitivity, Decision Log on `expect_type` continuity, optional `\dontrun{}` example update). Surprises & Discoveries populated with reviewer's empirical findings (`seed = NULL` path, `expect_named` ordering, `R = 2` edge case, `expect_type` typeof behavior, two pre-existing NOTEs). Total source change estimate revised from ~150 to ~160 lines including the bundled `.Rbuildignore` lines.

- **2026-05-10 (CRAN-readiness gate adopted)** — Greg directed that every PR cycle leave the package CRAN-ready, not just CRAN-bound ones. Workflow docs (`.workflow/PLANS.md`, `.workflow/REVIEW_AGENT.md`, `.workflow/R_WORKFLOW.md`, `.workflow/DEV_INSTRUCTIONS.md`) updated to make `devtools::check(args = "--as-cran")` + `devtools::spell_check()` + `urlchecker::url_check()` the unconditional per-PR gate, with `devtools::build_vignettes()` required only when vignettes are touched. This plan's Step 5 and Validation/Acceptance sections updated accordingly: spell_check / url_check added; criterion list re-numbered; commit list now includes `inst/WORDLIST` if first-run spell_check requires bootstrapping it. `build_vignettes()` is skipped for this PR (vignettes untouched).

- **2026-05-10 (implementation + post-execution review)** — Implementation committed as `3688aff` (12 files, +310/-38 lines net). Per-PR CRAN gate verified clean: `--as-cran` 0/0/0, 1065/1065 tests pass, `spell_check` empty (after 36-term WORDLIST bootstrap), `url_check` clean modulo 2 pre-existing redirects in vignettes. Post-execution review subagent verdict: LGTM, no blockers, no action items. Workflow doc cleanup applied during this round: `.workflow/PLANS.md:171` "preserve tabs by hand" warning corrected — `air format .` is empirically a no-op in this repo (auto-detects tabs).
