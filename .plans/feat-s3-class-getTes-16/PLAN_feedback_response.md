# PLAN.md (feat/s3-class-getTes-16) — feedback response (round 1)

Response to `PLAN_feedback.md` round 1. The reviewer found no blockers and six minor action items. Greg directed: keep item 1's guard for defensive-resiliency reasons; apply items 2–6; additionally clean up `.Rbuildignore` so the pre-existing NOTEs (item 4 / cosmetic 4) actually disappear instead of being merely re-categorized as "pre-existing".

## Item-by-item

**Streamlining 1 — `length(tes) >= 2` guard.** *Partially apply.* Greg directed to keep the guard for resiliency in case `genCoefs()`'s `R >= 2` validation ever changes upstream. Apply the cosmetic half of the recommendation: change the print-method string from `"<NA (R<2)>"` to `"<NA>"` (the `(R<2)` annotation would mislead users since `R<2` is unreachable through the current public API). Also add a one-line code comment in `summary.FETWFE_tes` explaining the guard is defensive. **Decision recorded.**

**Streamlining 2 — note that `cran-comments.md` is not bumped.** *Apply.* Adding a one-sentence note in Step 4 of "Concrete Steps."

**Robustness 1 — code comment on Test 6 regex sensitivity.** *Apply.* Adding the comment above Test 6 in the planned new test code.

**Robustness 2 — Decision Log entry on `expect_type` continuity.** *Apply.* Adding a Decision Log entry plus a Surprises & Discoveries entry recording the empirical fact (`typeof(structure(list(), class = ...)) == "list"`) the reviewer verified.

**Cosmetic 4 — reword "0 notes" acceptance criterion.** *Supersede with a stronger fix.* Greg directed to add the offending paths to `.Rbuildignore` so the NOTEs disappear entirely. New scope addition: bundle a small `.Rbuildignore` cleanup into this PR. Five lines added: `^\.workflow$`, `^\.claude$`, `^\.plans$`, `^AGENTS\.md$`, `^paper_arxiv\.tex$`. With those in place, `devtools::check()` and `devtools::check(args = "--as-cran")` should produce a clean `0 errors / 0 warnings / 0 notes`. Acceptance criterion 1 stays as "0 notes" (not "no new notes"), now achievable.

**Cosmetic 3 (optional) — add `print(te_results)` to `\dontrun{}` example.** *Apply.* Cheap and demonstrates the new method on the rendered help page.

**Other comments from feedback (not action items).** Reviewer's empirical confirmations of the `seed = NULL` path, `expect_named` ordering, regex escapes, and `R = 2` edge case are recorded in the plan's Surprises & Discoveries section so future readers see the verifications rather than re-deriving them.

## Scope expansion: `.Rbuildignore` cleanup

This PR now bundles a small `.Rbuildignore` cleanup alongside the original #16 work. Rationale:

- The two pre-existing CRAN NOTEs the reviewer identified are real and are noise on every `devtools::check()` going forward. Cleaning them up is small (5 lines added to one file), CRAN-friendliness work, and leaving them deferred means the next CRAN submission will need to address them anyway.
- Splitting into two separate PR cycles (#16 + a follow-up `chore/rbuildignore`) doubles the workflow overhead for ~5 lines of edit. Bundling is cheaper and Greg confirmed this approach in chat.
- The PR description will clearly call out the bundled scope ("Resolves #16; also adds five lines to `.Rbuildignore` to silence pre-existing CRAN NOTEs").

The bundled change does not affect the package's installed surface — `.Rbuildignore` only controls what's included in the source tarball.

## Convergence

Per `.workflow/PLANS.md` § "Plan review and iteration":

> If any change applied to the plan is substantive — any blocker, or any non-trivial streamlining or robustness item — invoke the subagent again on the updated plan.

None of the six original items is substantive, and the `.Rbuildignore` addition is mechanically simple (lines into a regex-list file with established conventions). No round 2 of plan review is needed; the executor proceeds to implementation after Greg's checkpoint review of the updated plan.

## Next step

Greg reviews the updated `PLAN.md` before implementation begins.
