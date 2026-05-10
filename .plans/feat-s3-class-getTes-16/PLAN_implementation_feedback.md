# PLAN.md (feat/s3-class-getTes-16) — post-execution review (round 1)

**Verdict: LGTM, ready for Greg's review.** No blockers, no robustness findings, two minor cosmetic notes. I ran the full per-PR CRAN gate against the working tree at commit `3688aff`: `devtools::document()` is idempotent, `devtools::check(args = "--as-cran")` reports `0 errors / 0 warnings / 0 notes`, `devtools::test()` reports `1065 / 1065` passing, `devtools::spell_check()` is empty, and `urlchecker::url_check()` reports only the two pre-existing redirects in `vignettes/vignette.Rmd:29` and `:31` that the user flagged as not introduced by this PR. I verified S3 dispatch from both `devtools::load_all()` and `R CMD INSTALL`ed tarball; output matches the plan's Purpose-section mockup exactly. I confirmed the new tests fail against `origin/main` baseline and pass against HEAD. The `^\.plans$` etc. `.Rbuildignore` entries correctly exclude the workflow files from the source tarball.

## 1. Blocking issues

None.

## 2. Streamlining opportunities

None.

## 3. Robustness / clarity concerns

None.

## 4. Cosmetic

1. **`FETWFE_tes-class.Rd` is minimal but acceptable.** The generated `man/FETWFE_tes-class.Rd` has only `\title{}` and `\description{}` — no `\alias{}` beyond the class name itself, no `\seealso{getTes}` linking back to the function. This matches the pattern of `man/fetwfe-class.Rd` and `man/etwfe-class.Rd`, which are also two-line stubs. The renders are sparse on `?FETWFE_tes-class` but consistent with the existing convention. No action needed; mentioned for completeness.

2. **`RoxygenNote: 7.3.2` → `7.3.3` is an unannounced side-effect of `devtools::document()` running under roxygen2 7.3.3.** This is expected when the executor's installed roxygen2 differs from the version that last ran on `origin/main`. Not a finding; both 7.3.2 and 7.3.3 produce the same `.Rd` output for this PR. Mentioned for transparency.

## 5. Verifications I performed

All commands were run against a `cp -R` of the repo at `/tmp/fetwfe-plan-review/pkg`, which was committed to `feat/s3-class-getTes-16` at `3688aff` (matching HEAD). The original repo at `/Users/gregfaletto/Documents/GitHub/fetwfePackage` was never modified.

**Per-PR CRAN gate.**
- `air format .` — NOT run, per `.workflow/PLANS.md:171` which explicitly states "No Styler config; no Air formatter config (despite an `air.toml` line in `.Rbuildignore`, the file does not exist — preserve tabs by hand)." Without an `air.toml`, `air format .` would convert all R-source tabs to spaces across the whole codebase. I instead manually verified that `R/getTes_class.R` uses hard tabs (matching the rest of `R/`) by inspecting with `awk '{gsub(/\t/, "<TAB>"); print}'`. The new file is correctly tab-indented.
- `devtools::document()` — ran twice. First run regenerated `man/FETWFE_tes-class.Rd`, `man/getTes.Rd`, and `NAMESPACE` matching the committed versions byte-for-byte. Second run produced no diff (idempotent). The pre-existing roxygen `@details` / `@seealso` warnings about `getPsiRFused`, `getSecondVarTermDataApp`, `getTeResults2` in `fetwfe_core.R` appeared in output but are not from this PR's diff.
- `devtools::check(args = "--as-cran")` — exits with `Status: OK`, `0 errors ✔ | 0 warnings ✔ | 0 notes ✔` in 30.4s.
- `devtools::test()` — exits with `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1065 ]`.
- `devtools::test(filter = "getTes")` — exits with `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 27 ]`, covering Tests 1–7.
- `devtools::spell_check()` — `No spelling errors found.` Empty data frame.
- `urlchecker::url_check()` — only the two pre-existing redirects at `vignettes/vignette.Rmd:29:218` and `vignettes/vignette.Rmd:31:107` (both bookdown.org → posit.cloud). Vignettes are not in this PR's diff.

**Function traces.**
- `R/gen_funcs.R::getTes()` (lines 357–397): math body unchanged from `origin/main`. Return-list assembly correctly carries `R`, `T`, `d`, `seed` from `coefs_obj`. `class(out) <- "FETWFE_tes"` lands after assembly. The `inherits(coefs_obj, "FETWFE_coefs")` boundary check is preserved; the `stopifnot(length(beta) == p)` check still fires before any class-attached return path.
- `R/getTes_class.R::print.FETWFE_tes` (lines 7–28): `cat()` calls produce labeled output matching the plan's Purpose-section mockup. `is.null(x$seed)` correctly handles both `seed = NULL` and `seed = 0` (verified: NULL → "<none>", 0 → "0"). `seq_along(x$actual_cohort_tes)` iterates correctly for any positive length. Returns `invisible(x)`. `sprintf("%d", x$R)` would error on a non-integer double, but `genCoefs()` enforces integer values for `R`/`T`/`d` so this is safe in practice.
- `R/getTes_class.R::summary.FETWFE_tes` (lines 31–51): assembles the documented field set in the order asserted by Test 7 (`min`, `max`, `median`, `sd`). The `length(tes) >= 2` guard is unreachable through `genCoefs()` (which requires `R >= 2`), retained per Greg's defensive-resiliency directive with the appropriate code comment.
- `R/getTes_class.R::print.summary.FETWFE_tes` (lines 54–88): mirrors `print.FETWFE_tes` with the dispersion block added between Cohort effects and Generated from. The `"<NA>"` literal (changed from `"<NA (R<2)>"` per round-1 feedback) is technically unreachable through the public API but the guard is correct.
- `tests/testthat/test-getTes.R` Test 5 (lines 137–152): asserts class is `FETWFE_tes` and carried params equal generating params. `expect_equal(res$R, 5)` (not `5L`) works because `genCoefs()` stores R/T/d/seed as doubles (verified empirically); both `5` and `5L` would pass with default tolerance.
- Test 6 (lines 160–179): regex-based assertions against `capture.output(print(...))`. Verified against baseline (origin/main without new methods): the test fails as expected because `print.default` produces `$att_true` / `$actual_cohort_tes` output that doesn't match any of the regexes. Against HEAD it passes. The robustness comment at lines 154–158 is in place per plan-review action item 3.
- Test 7 (lines 184–216): asserts `summary.FETWFE_tes` class, the dispersion-stats field set (order-sensitive `expect_named`), and that min/max/median/sd values agree with direct computation on `actual_cohort_tes`. Verifies the print-summary dispersion block via `capture.output`. Passes against HEAD; fails against baseline.

**S3 dispatch.**
- Verified from `devtools::load_all()`: `class(getTes(genCoefs(R = 3, T = 5, d = 2, density = 0.5, eff_size = 1, seed = 1)))` returns `"FETWFE_tes"`; `print(...)` dispatches to `print.FETWFE_tes`; `summary(...)` dispatches to `summary.FETWFE_tes`; `print(summary(...))` dispatches to `print.summary.FETWFE_tes`.
- Verified from `R CMD INSTALL`ed tarball at `/tmp/fetwfe-plan-review/Rlib`: same dispatch behavior. Print output byte-identical to load_all version.
- Verified `NAMESPACE` lines 5, 8, 11 (`S3method(print,FETWFE_tes)`, `S3method(print,summary.FETWFE_tes)`, `S3method(summary,FETWFE_tes)`) match the three `#' @export` tags at `R/getTes_class.R:6`, `:30`, `:53`.

**Backward compatibility.**
- `expect_type(res, "list")` in existing Test 1 still passes (verified directly): `typeof()` is unaffected by S3 class assignment.
- `vignettes/simulation_vignette.Rmd` uses `$att_true` and `$actual_cohort_tes` directly (lines 125–161) and never calls `print(true_tes)`. Vignette output is byte-identical to baseline.

**`.Rbuildignore` cleanup.**
- Ran `R CMD build .` against `/tmp/fetwfe-plan-review/pkg`; ran `tar -tzf fetwfe_1.5.1.tar.gz | grep -E '\.workflow|\.claude|\.plans|AGENTS|paper_arxiv'`. No matching entries in the tarball. The five new `^foo$` regex lines correctly exclude `.workflow/`, `.claude/`, `.plans/`, `AGENTS.md`, and `paper_arxiv.tex` from the source tarball.

**`inst/WORDLIST` correctness.**
- Verified `inst/WORDLIST` is the canonical location for the `spelling` package via `print(spelling:::get_wordfile)`. The 36 entries are loaded by `spelling::get_wordlist(pkg = ".")`.
- Confirmed `FETWFE’s` (with curly apostrophe U+2019) appears in `vignettes/simulation_vignette.Rmd` and is therefore correctly entered in WORDLIST.
- No spell_check flags. None of the new strings introduced in `R/getTes_class.R` ("True Treatment Effects", "Overall true ATT", "Cohort effects", "Cohort effect dispersion", "Generated from", "Cohorts (R)", "Time periods (T)", "Covariates (d)", "Seed") require WORDLIST entries because they are inside `cat()` calls in `R/`, and `spell_check()` only scans DESCRIPTION text, NEWS.md, README, vignettes, and `man/*.Rd` files — not R source.

**Meta-file integrity.**
- `DESCRIPTION:3` `Version: 1.5.1` (was 1.5.0).
- `DESCRIPTION:15` `RoxygenNote: 7.3.3` (was 7.3.2). Side-effect of `devtools::document()` under roxygen2 7.3.3; harmless.
- `NEWS.md:3–7` new `## Version 1.5.1 (2026-05-10)` section.
- `inst/CITATION:6` `note = "R package version 1.5.1"` (was 1.5.0). Three version strings agree.
- `cran-comments.md` correctly NOT bumped (per round-1 feedback decision; tracks CRAN submission state, not local version).

**Commit hygiene.**
- `git show --name-only 3688aff` lists exactly 12 files (matching the plan's Step 6 list, with `.plans/feat-s3-class-getTes-16/PLAN.md` re-staged for the CRAN-gate adoption update). No `.Rhistory`, `.DS_Store`, or other personal-state files. No `.claude/`, no `Meta/`, no `doc/`.

**Anti-pattern scan.**
- `grep -rn "TODO\|FIXME\|XXX" R/getTes_class.R tests/testthat/test-getTes.R` — empty. No new `TODO`/`FIXME`/`XXX` markers.
- `grep 'stop("not implemented"\|.NotYetImplemented\|.NotYetUsed' R/` — empty.
- No `:::` accesses in the new file.
- No unguarded `solve()` (none used in the new code).
- No `Suggests:` package usage in the new code.
- No new `@inheritParams`; no stale ones detected.
- No hand-edits to `man/*.Rd` or `NAMESPACE` (verified via `devtools::document()` idempotence).

**Plan-feedback alignment.**
- Plan-review item 1 (`length(tes) >= 2` guard): kept with code comment at `R/getTes_class.R:33–35`, user-facing string changed from `"<NA (R<2)>"` to `"<NA>"` at line 73. Matches feedback response.
- Plan-review item 2 (`cran-comments.md` note): present in the plan text. Not applicable to source diff.
- Plan-review item 3 (Test 6 robustness comment): present at `tests/testthat/test-getTes.R:154–158`.
- Plan-review item 4 (Decision Log on `expect_type`): present in the plan's Decision Log at lines 108–110 of `PLAN.md`.
- Plan-review item 5 (`.Rbuildignore` cleanup): applied; 5 new lines.
- Plan-review item 6 (`\dontrun{}` example): `print(te_results)` line added at `R/gen_funcs.R:353`, propagated to `man/getTes.Rd:57`.

## 6. What I did NOT verify

- `air format .` — skipped (would mass-convert tabs to spaces across the codebase because no `air.toml` exists in the repo). Verified tab-indentation by direct inspection instead. This is a workflow-doc inconsistency (PLANS.md:171 says "preserve tabs by hand" but REVIEW_AGENT.md § 2 prescribes running `air format .`); not a finding on this PR, but worth flagging separately as a workflow cleanup.
- `devtools::build_vignettes()` standalone — skipped per the plan's acceptance criterion #4 (vignettes/ untouched); the vignettes were nonetheless rebuilt as part of `R CMD check --as-cran` and passed.
- Cross-reference against `paper_arxiv.tex` — not applicable; this PR makes no math changes.
- Manual test of the `cat(sprintf("%d", ...))` failure path for non-integer `R`/`T`/`d` — these cannot flow through `genCoefs()` (which enforces integer inputs), so the failure path is unreachable. Did not construct a synthetic `FETWFE_coefs` object with non-integer R/T/d to probe.
- The two pre-existing URL redirects in `vignettes/vignette.Rmd` — confirmed they appear in `urlchecker::url_check()` output but did not investigate fixing them (out of scope per input).
- The three pre-existing `@details`/`@seealso` roxygen warnings in `fetwfe_core.R` — confirmed they appear in `devtools::document()` output but did not investigate fixing them (out of scope per input).

## 7. Summary of action items

None required. The implementation is ready for Greg's review.

After Greg's checkpoint, push the branch and open the PR. `devtools::check(args = "--as-cran")` is clean, all tests pass, spell_check and url_check are clean modulo the two pre-existing redirects flagged in the PR description.
