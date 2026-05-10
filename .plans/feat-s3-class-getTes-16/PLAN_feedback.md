# PLAN.md (feat/s3-class-getTes-16) — review (round 1)

Verdict: **LGTM with two small robustness/streamlining items**. The plan is implementable as written and produces the documented behavior. I overlaid the planned `R/getTes_class.R` and the patched `getTes()` body onto a `devtools::load_all()`ed copy of the package, ran `print()`/`summary()` on `getTes(genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1, seed = 1))`, ran the three planned testthat blocks (all 20 assertions pass), regenerated `man/` and `NAMESPACE` from the planned roxygen in a `/tmp` copy of the package, ran `devtools::check()` on the patched copy, and confirmed the four existing `getTes` tests continue to pass. No blockers. The `seed`-handling guard, the `expect_named` ordering, the `expect_type(res, "list")` continuity, and the regex escapes in Test 6 are all empirically correct.

## 1. Blocking issues

None.

## 2. Streamlining opportunities

1. **`length(tes) >= 2` guard in `summary.FETWFE_tes` is dead code.** `R/gen_funcs.R::genCoefs()` (lines 242–246) already requires `R >= 2` with `stop("R must be a numeric value greater than or equal to 2 ...")`, and `getTes()` (line 349) requires its input to be a `FETWFE_coefs`. So `length(actual_cohort_tes)` is always `>= 2` and the `else NA_real_` branch is unreachable. Either:
   - Drop the guard: `sd = stats::sd(tes)`, and remove the `is.na(...)` / `"<NA (R<2)>"` branch in `print.summary.FETWFE_tes`. This is cleaner and removes a never-tested branch.
   - Or, if you want the defense-in-depth, add a one-line code comment explaining why the guard exists despite `genCoefs` validation, so a future reader doesn't delete it.

   Cosmetic note if you keep the guard: the `print.summary.FETWFE_tes` literal `"<NA (R<2)>"` will never appear in any test output since R<2 is unreachable through the public API. Test 7 doesn't cover it either.

2. **Plan does not mention `cran-comments.md`.** It currently reads `**Version:** 1.5.0`. Per `AGENTS.md`, `cran-comments.md` is updated at CRAN-submission time, not at every patch bump. The plan correctly omits it for this PR. Mention this explicitly in the plan's Concrete Steps so a future contributor doesn't ask why; a one-sentence note in Step 4 ("`cran-comments.md` is intentionally not bumped — it tracks CRAN submission state, not local version") is sufficient.

## 3. Robustness / clarity concerns

1. **Test 6 is brittle to cosmetic relabeling.** Test 6's regexes (`"Cohorts \\(R\\)"`, `"Time periods \\(T\\)"`, `"Covariates \\(d\\)"`) match the exact label format. I confirmed empirically that changing "Cohorts (R)" to "Cohorts: R" in the print method body breaks Test 6 with a regex-mismatch failure. This is fine in the absolute (a behavior test ought to fail when behavior changes), but the prompt flagged it explicitly. Two cheaper options if you want robustness:
   - Replace the parenthesized-label regexes with substring matches: `"Cohorts"`, `"Time periods"`, `"Covariates"`. Keeps Test 6 tied to "the print method names the four parameters" without locking the exact label format.
   - Or split Test 6 into two: one for the must-have content ("Overall true ATT", "Cohort 1", "Cohort 2", "Cohort 3", "Seed") and one (golden-output style) for the exact label format. The current single-test mix is fine; the proposed split is overkill for two pieces of data.

   **My recommendation:** the existing form is acceptable; document in a code comment above Test 6 that the regexes intentionally lock the print-output labels and that any cosmetic change to `print.FETWFE_tes` requires updating the test.

2. **Plan says "the `expect_type(res, "list")` test still passes" but did not state how it was verified.** I verified empirically: `typeof(structure(list(), class = "FETWFE_tes"))` returns `"list"`, so `expect_type(res, "list")` passes (this is testthat's typeof-based assertion, distinct from `expect_s3_class`). Worth a one-line note in the plan's Decision Log: "Verified that `expect_type(res, "list")` still passes after class assignment because `typeof()` is unaffected by S3 class." This is the kind of empirical fact future readers will want to see.

3. **`R = 2` edge case is silently fine, but the plan's `length(tes) >= 2` guard is misleading.** I tested `genCoefs(R = 2, T = 4, d = 0, density = 0.5, eff_size = 5, seed = 1)`: the resulting cohort vector has length 2 and `stats::sd(tes)` returns a real number (~2.36 in my run). The guard's *threshold* is right (`sd` of a length-1 vector is `NA`), but per item 1 above, the guard never triggers for any input that `genCoefs` accepts. If you keep the guard, change the user-facing string in `print.summary.FETWFE_tes` from `"<NA (R<2)>"` to `"<NA>"`; the (R<2) annotation will mislead users into thinking R<2 is somehow possible.

## 4. Cosmetic

1. **`@importFrom stats median` vs. `stats::median(tes)`.** The plan calls `stats::median` qualified to avoid touching `@importFrom`. This is consistent with `R/core_funcs.R` (which uses `stats::predict`, `stats::coef`), `R/etwfe_core.R` (`stats::model.matrix`), and `R/ols_calcs.R`/`R/fetwfe_core.R` (`stats::qnorm`). The package mixes both styles (top-level `@importFrom stats` blocks at `R/fetwfe.R:2` and `R/utility.R:4`, plus per-call `stats::` qualifications). Either is consistent with the repo. Plan's choice is fine.

2. **Decision Log Surprises & Discoveries: empty.** As of round 1 the plan's "Surprises & Discoveries" reads "(None yet.)". After this review, you'll likely add at least one entry (e.g., the `expect_type(res, "list")` continuity verification or the `R=2` guard discussion). That's normal — flag here so it's not forgotten.

3. **Plan example mockup is good but the existing `\dontrun{}` `@examples` block in `getTes` roxygen is not updated.** The current example calls `print(te_results$att_true)` (a numeric scalar). It does not call `print(te_results)` directly. After the change, the example still works but does not demonstrate the new method. A one-line addition like

       # Or use the new print method for a self-describing display:
       print(te_results)

   in the existing `\dontrun{}` block (which already wraps in `\dontrun{}`, so no CRAN example-time concern) would showcase the new behavior without expanding the example's runtime budget. Optional; not a finding if you skip it. (Note: with `\dontrun{}` the example doesn't run during `R CMD check`, so this is purely for the rendered help page.)

4. **Pre-existing `devtools::check()` NOTEs are not the plan's responsibility.** `devtools::check()` on `origin/main` produces 2 NOTEs (one for `.workflow/.claude/.plans` hidden directories, one for `AGENTS.md`/`paper_arxiv.tex` top-level files). The plan claims acceptance criterion "0 errors / 0 warnings / 0 notes", which is a stronger claim than what currently holds. Either:
   - Reword acceptance criterion 1 to "no NEW errors/warnings/notes introduced by this PR."
   - Or accept that the local `devtools::check()` will show 2 pre-existing NOTEs that are unrelated to this work.

   This is a wording fix to the plan, not a behavior change.

## 5. Verifications I performed

All against `db066ac734a1de2c28666570660676a7fb91e164` (origin/main).

- Read the plan, the clarified-scope paragraph, `R/fetwfe_class.R`, `R/etwfe_class.R`, `R/gen_funcs.R::getTes()` and `genCoefs()`, `tests/testthat/test-getTes.R`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, `inst/CITATION`, `cran-comments.md`, `vignettes/simulation_vignette.Rmd` (all `getTes` callsites).
- Confirmed there are no other `getTes` callsites beyond `vignettes/simulation_vignette.Rmd` (which uses `$att_true` / `$actual_cohort_tes` only — not `print(true_tes)` — so the new methods don't change vignette output).
- Wrote planned `R/getTes_class.R` to `/tmp/fetwfe-plan-review/getTes_class_overlay.R` and the planned new `getTes()` body to the same file. Wrote planned tests 5/6/7 to `/tmp/fetwfe-plan-review/test_getTes_planned.R`.
- Ran `devtools::load_all("/Users/gregfaletto/Documents/GitHub/fetwfePackage")`, sourced the overlay into the global env (the loaded fetwfe namespace was locked, so I overlaid in `globalenv()` and called `registerS3method()` for the three methods to mimic how roxygen-generated `S3method(...)` would register them). Then exercised `getTes()`, `print()`, and `summary()` on `genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1, seed = 1)`. Output:

      True Treatment Effects (from FETWFE_coefs)
      ==========================================
      Overall true ATT: 0.0556

      Cohort effects:
        Cohort 1: 0.0000
        Cohort 2: 0.6667
        Cohort 3: -0.5000

      Generated from:
        Cohorts (R)      : 3
        Time periods (T) : 5
        Covariates (d)   : 2
        Seed             : 1

  Matches the mockup in the plan's Purpose section exactly.

- `summary(...)` returns `summary.FETWFE_tes` with `cohort_te_stats` named `c("min", "max", "median", "sd")` in that order; `min`/`max`/`median`/`sd` values agree with `min`/`max`/`stats::median`/`stats::sd` of the cohort vector. `print(summary(...))` produces the dispersion block as described.
- Confirmed `seed = NULL` path: `genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1)` (no seed argument) produces `coefs$seed == NULL`. The planned `is.null(x$seed)` guard correctly renders "<none>". Confirmed `seed = 0` path renders "0" (not "<none>"). Confirmed an "old-style" `FETWFE_coefs` object with no `seed` field at all also renders "<none>" (because `coefs_obj$seed` returns `NULL` for missing names).
- Ran `testthat::test_file("/tmp/fetwfe-plan-review/test_getTes_planned.R")` after the overlay. All 20 assertions pass (3 test_that blocks).
- Copied the package to `/tmp/fetwfe-plan-review/pkg`, applied the planned `getTes()` body change and added the planned `R/getTes_class.R` (and the planned new tests). Ran `devtools::document()` — produced exactly the three expected `S3method(...)` lines and a new `man/FETWFE_tes-class.Rd` file matching the plan's Artifacts section. Ran `devtools::document()` a second time — no diff (idempotent). Ran `devtools::test(filter = "getTes")` and `devtools::test(filter = "getTes_planned")` — both pass cleanly. Ran `devtools::check()` — `0 errors, 0 warnings, 2 notes`, where both notes are pre-existing (hidden directories and `AGENTS.md`/`paper_arxiv.tex`).
- Confirmed planned tests fail on baseline `origin/main` (without the plan's source changes). All three test blocks fail because (a) `expect_s3_class(res, "FETWFE_tes")` fails on a plain list, and (b) `print.FETWFE_tes` and `summary.FETWFE_tes` don't exist so `print(getTes(coefs))` falls through to `print.default` and emits raw `$att_true` / `$actual_cohort_tes` lines that fail the regex assertions.
- Verified `expect_type(res, "list")` still passes after class is added: `typeof()` is unaffected by S3 class assignment, and testthat's `expect_type` calls `typeof()`. So existing Test 1 continues to pass.
- Verified `expect_named(s$cohort_te_stats, c("min", "max", "median", "sd"))` is order-sensitive (different order fails) and that `c(min = ..., max = ..., median = ..., sd = ...)` preserves declaration order. So Test 7 passes against the planned construction.
- Verified regex escapes `"Cohorts \\(R\\)"` etc. work as intended in `expect_match` — they match a literal "(R)" in the print output.
- Probed Test 6 robustness: confirmed the regex `"Cohorts \\(R\\)"` fails after a cosmetic relabeling to "Cohorts: R" (item 3.1 above).
- Verified `R = 2` edge case: `length(actual_cohort_tes) == 2`, and `stats::sd(tes)` returns a real number (the plan's `>= 2` guard does not trigger). Confirmed the guard is correct in principle but unreachable through `genCoefs()` (which requires `R >= 2`).
- Verified `sprintf("%d", coefs$R)` works for the integer-valued doubles that `genCoefs` produces (R, T, d are stored as numeric but always whole-number positive). The pattern errors on non-integer doubles, but the values that flow through `getTes()` cannot be non-integer.
- Verified the `vignettes/simulation_vignette.Rmd` continues to work unchanged: it accesses `true_tes$att_true` and prints `true_tes$actual_cohort_tes` directly, never calls `print(true_tes)` — so the new method does not change vignette output.

## 6. What I did NOT verify

- **`devtools::check(args = "--as-cran")`** — I ran `devtools::check()` but not the `--as-cran` variant, which can surface CRAN-specific NOTEs (chatty examples, non-ASCII, etc.). Given the plan adds only `cat()` calls inside `print()` methods (consistent with existing `print.fetwfe`/`print.etwfe` precedent) and no new examples, I expect `--as-cran` to behave the same as plain `check()`, but I did not confirm.
- **Vignette rebuild.** I did not run `devtools::build_vignettes()`. The vignette code does not call `print(getTes(...))`, so the rebuild output should be byte-identical to baseline, but I did not verify this empirically.
- **Greg's actual approval criteria for the print mockup.** The plan's Purpose-section mockup is what I verified the print method produces. Whether Greg likes the exact label spacing (e.g., "  Cohorts (R)      : 3" with 6 spaces of padding between "(R)" and ":") is a style call; I did not second-guess it.
- **The `inst/CITATION` and `NEWS.md` edits themselves.** The plan describes the version bumps in plain prose (Step 4); I did not apply them in the temp package or run `tools::checkRdaFiles()` on the resulting tarball.
- **The `air format .` step.** I did not run `air format` on the planned files. The plan notes the displayed code is space-indented for readability and tells the implementer to convert to tabs before saving / let `air format` normalize. I did not verify that `air format .` produces no diff against tab-indented planned code.

## 7. Summary of action items

In priority order; none are blockers.

1. (Streamlining 1, optional but recommended.) Drop the `length(tes) >= 2` guard in `summary.FETWFE_tes` and the corresponding `is.na(...)`/`"<NA (R<2)>"` branch in `print.summary.FETWFE_tes`, since `genCoefs()` already requires `R >= 2` and the branch is unreachable. **Or** keep the guard with a one-line code comment explaining why and change the user-facing string from `"<NA (R<2)>"` to `"<NA>"`.

2. (Streamlining 2.) Add a one-sentence note in Step 4 of "Concrete Steps" that `cran-comments.md` is intentionally not bumped because it tracks CRAN-submission state, not local version.

3. (Robustness 1.) Add a code comment above Test 6 in the planned new test code explaining that the `"Cohorts \\(R\\)"` / `"Time periods \\(T\\)"` / `"Covariates \\(d\\)"` regexes intentionally lock the print-output labels, so a cosmetic relabeling of `print.FETWFE_tes` requires updating the test.

4. (Robustness 2.) Add a one-line entry to the Decision Log noting that `expect_type(res, "list")` was empirically verified to still pass after the class assignment (because `typeof()` is independent of S3 class).

5. (Cosmetic 4.) Reword acceptance criterion 1 in "Validation and Acceptance" from "0 errors ✔ | 0 warnings ✔ | 0 notes ✔" to something like "no NEW errors, warnings, or notes introduced by this PR." On `origin/main`, `devtools::check()` already produces 2 NOTEs (hidden dirs `.workflow/.claude/.plans`, top-level files `AGENTS.md`/`paper_arxiv.tex`) that are pre-existing and not this PR's responsibility.

6. (Cosmetic 3, optional.) Update the `\dontrun{}` `@examples` block in `getTes` roxygen to add a single line `print(te_results)` after the existing `print(te_results$att_true)` line, so the rendered help page demonstrates the new print method.

After (1)–(5) the plan is ready for execution. Items (1) and (3)–(5) are local edits to the plan text only; (6) is an additional change to the `R/gen_funcs.R` `getTes` roxygen at execution time. None require another plan-review round; the executor can apply them inline and proceed.
