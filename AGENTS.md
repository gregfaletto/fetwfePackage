# AGENTS.md

Guidance for AI coding agents working in this repository.

## What this is

`fetwfe` is a CRAN-published R package implementing **fused extended two-way
fixed effects** (FETWFE) and several related estimators for difference-in-
differences with staggered treatment adoption. The methodology is defined in
`paper_arxiv.tex` (Faletto 2025, arXiv:2312.05985); the package is the
reference implementation.

The estimator targets the overall ATT (average treatment effect on the treated)
and CATTs (cohort-specific ATTs) on a balanced panel where treatment is an
absorbing state. It applies a bridge-regression (`L_q`, `0 < q < 1`) fusion
penalty to the ETWFE parameterization to gain efficiency while still admitting
valid standard errors.

If you need methodological context, read `paper_arxiv.tex` (or the arXiv link
in `README.md`) — do not infer the math from the code alone. The code
implements specific equations and lemmas from that paper.

## Working in this repo as a coding agent

For any non-trivial PR, follow the ExecPlan workflow defined in `.workflow/`.
The directory is `.gitignore`d (it is the maintainer's personal workflow
configuration) but it is load-bearing for how the repo is meant to be
developed: it specifies the plan format, the multi-subagent review cycle, the
per-PR CRAN gate, and the recurring failure modes that PRs need to defend
against. The mandate is not optional — `.workflow/PLANS.md` says ExecPlans
must be followed "to the letter" for any plan that lands in this repo.

**"Non-trivial" means anything matching the criteria in
`.workflow/IMPLEMENTER_AGENT.md` § "When to use this subagent":** multi-file
refactors, any change ≥ 50 lines, any change that would generate ≥ 5000
tokens of tool output during implementation, math-touching changes, or new
exported functions. The exemption is small in-thread fixes — ≤ 10-line bug
fixes, trivial CRAN-NOTE cleanups, single-line typo fixes — and even those
must satisfy the per-PR CRAN gate including `air format .`,
`devtools::spell_check()`, and `urlchecker::url_check()`.

The `.workflow/` directory contains the following documents. Read them when
the task at hand falls in their scope:

- `PLANS.md` — ExecPlan format and the per-PR acceptance criteria. The
  authoritative specification for any plan that lands in this repo. The
  living-document sections (`Progress`, `Surprises & Discoveries`,
  `Decision Log`, `Outcomes & Retrospective`) are mandatory and must be kept
  current as work proceeds. The plan filename convention is `PLAN.md`
  (in `.plans/<branch>/`), not `SPEC.md`.
- `CHOOSING_A_PR.md` — how to pick the next target plus the
  **issue-clarification protocol** that is mandatory when working from
  underspecified GitHub issues (which most fetwfe issues are). Read before
  drafting any plan.
- `DEV_INSTRUCTIONS.md` — branch/PR git workflow targeting `main`, PR
  description style, the no-fork convention, and the requirement to draft
  the PR description at `.plans/<branch>/<branch>_pr_description.md`
  before opening the PR.
- `R_WORKFLOW.md` — daily R + RStudio/VS Code edit/test loop, common error
  decoding (rank-condition failures, S3 dispatch errors, `Matrix` /
  `glmnet` / `grpreg` shape mismatches), and the per-PR CRAN gate's exact
  command sequence.
- `IMPLEMENTER_AGENT.md` — implementer subagent role specification. Defines
  what work is implementer-scope vs. in-thread orchestrator-scope, the
  smoke-test cadence, and the escalation criteria.
- `REVIEW_AGENT.md` — post-execution review checklist. Every non-trivial PR
  gets reviewed by a subagent against this specification before opening
  for Greg.
- `SENTINEL_AGENT.md` — drift sentinel that runs twice per cycle (pre- and
  post-implementation). Catches three specific drift classes: copy-paste
  duplication across sibling files, missing `.check_for_<method>(x)`
  preconditions on new methods that read from estimator class objects, and
  anti-pattern recurrence in tests (tautological round-trips, fixtures
  that bypass the code path under test, well-formedness-only assertions).
- `PERIODIC_CODE_REVIEW.md` — quarterly multi-agent review pass for
  cross-PR structural concerns. Not a per-PR gate; run on a separate
  cadence.
- `PRIORITIZATION.md` — tier order for choosing among candidate PRs
  (inference bugs → workflow → other bugs → simpler-code refactors →
  robustness refactors → docs → periodic review → new features) plus
  heuristics for choosing within a tier.
- `WORKFLOW_LESSONS.md` — failure modes that have actually bitten past
  contributions (NAMESPACE drift after roxygen edits, version-string
  drift across `DESCRIPTION` / `NEWS.md` / `inst/CITATION`, copy-paste
  drift between sibling files, etc.). Plans should pre-empt these.

The standard four-subagent cycle for a non-trivial PR is: (1) planning
reviewer reviews the draft `PLAN.md` per `PLANS.md` § "Plan review and
iteration"; (2) drift sentinel runs against the plan; (3) implementer
applies the plan per `IMPLEMENTER_AGENT.md`; (4) drift sentinel and
post-execution reviewer run against the diff in parallel per
`SENTINEL_AGENT.md` and `REVIEW_AGENT.md`. Each review round writes to a
file and iterates to convergence; never overwrite earlier rounds.

## Estimators exported by the package

There are four estimators, each with a `*WithSimulatedData()` wrapper that
takes a `FETWFE_simulated` object from `simulateData()`:

| Function     | What it is                                                | Use it for                           |
|--------------|-----------------------------------------------------------|--------------------------------------|
| `fetwfe()`   | The main FETWFE estimator (bridge-fused ETWFE).           | The recommended estimator.           |
| `etwfe()`    | Plain extended TWFE (Wooldridge 2021-style).              | Baseline / comparison.               |
| `betwfe()`   | Bridge-penalized ETWFE (no fusion transformation).        | Comparison / ablation.               |
| `twfeCovs()` | Standard TWFE with covariates and per-cohort effects.     | Simulation only — biased; do **not** recommend for applied use. The roxygen header says so explicitly. |

Other exported helpers: `genCoefs()` / `genCoefsCore()`, `simulateData()` /
`simulateDataCore()`, `getTes()`, `attgtToFetwfeDf()`, `etwfeToFetwfeDf()`.
The full export list lives in `NAMESPACE` (regenerated by roxygen — don't edit
by hand).

## Repo layout

```
R/                  Package source, grouped by role (one file per unit):
  # Entry points + estimator cores
    fetwfe.R              fetwfe()/etwfe() user API + wrappers; @import block.
    fetwfe_core.R         Input checks + numerical core for fetwfe().
    etwfe_core.R          Input checks + core for etwfe/betwfe/twfeCovs.
    betwfe_core.R         betwfe() API + bridge-penalized ETWFE core.
    twfeCovs.R            twfeCovs() API + core (biased; simulation only).
  # Shared numerical machinery
    input_prep.R          Input validation, design prep, cohort-probability prep.
    gls_machinery.R       Variance/GLS whitening, ridge rows, Gram inverse, Omega^-1/2.
    result_assembly.R     Assemble the final selected-coefficient output object.
    design_matrix.R       Build ETWFE design: covariate/FE/treatment terms.
    fusion_transforms.R   Fusion transforms D_N^-1 (cohort & event-study).
    bridge_selection.R    Bridge back-transform + BIC/CV lambda selection.
    variance_machinery.R  Cohort/overall-ATT variance; cluster-robust SEs.
    cluster_floor.R       PSD floor for the cluster-robust quadratic form.
    process_covs.R        Normalize the covs argument to a character vector.
    utility.R             Misc helpers: idCohorts, input checks, my_scale, etc.
  # S3 classes (methods + validators)
    fetwfe_class.R        print/summary/coef + .validate_fetwfe (class fetwfe).
    etwfe_class.R         print/summary/coef + .validate_etwfe (class etwfe).
    betwfe_class.R        print/summary/coef + .validate_betwfe (betwfe).
    twfeCovs_class.R      print/coef + .validate_twfeCovs (class twfeCovs).
    getTes_class.R        print/summary for FETWFE_tes (true effects).
    sim_classes.R         print for FETWFE_simulated / FETWFE_coefs.
    catt_df_class.R       catt_df column-rename map + helpful errors.
    class_helpers.R       .check_for_* preconditions, validators, printers.
  # Inference + reporting accessors
    event_study.R         eventStudy(): pooled event-time estimates.
    cohort_study.R        cohortStudy(): per-cohort ATTs.
    simultaneous_cis.R    simultaneousCIs() generic + methods (max-t bands).
    plot.R                plot() methods for CATT / event-study estimates.
    broom_methods.R       tidy/glance/augment for all estimator classes.
  # Simulation, data generation + conversion
    gen_coefs.R           genCoefs/genCoefsCore/getTes: coefs + true effects.
    gen_data.R            simulateData/simulateDataCore: panel generation.
    sim_helpers.R         Internal simulation draws (covariates, FE, etc.).
    cohort_assignment.R   Covariate-dependent cohort-assignment DGP.
    convert_dfs.R         attgtToFetwfeDf / etwfeToFetwfeDf converters.

man/                roxygen2-generated .Rd files. NEVER edit directly.
NAMESPACE           roxygen2-generated. NEVER edit directly.
DESCRIPTION         Package metadata; bump Version here for releases.
NEWS.md             Per-version changelog (user-visible).
tests/testthat/     testthat suite, one test-*.R per exported function.
vignettes/          knitr/rmarkdown vignettes (vignette.Rmd is the main one).
inst/CITATION       Citation metadata; keep version in sync with DESCRIPTION.
cran-comments.md    Submission notes for the current CRAN release.
CRAN-SUBMISSION     CRAN-managed; don't hand-edit.
paper_arxiv.tex     Source for the methodology paper (R-build-ignored).
```

## Code conventions

- **Indentation: hard tabs.** Every R file in `R/` uses tab indentation. Do
  not introduce spaces.
- **Documentation: roxygen2 (markdown), on every function.**
  `Roxygen: list(markdown = TRUE)` is set in DESCRIPTION.
  **Every function — public or internal — gets a roxygen block** with a
  one-line title, `@param` for each argument, and `@return`. Public
  functions additionally get `@export` and `@examples`. Internal helpers
  additionally get `@keywords internal` and `@noRd` (so they're documented
  for code readers but don't produce a user-facing man page). See
  `R/utility.R::idCohorts`, `R/gen_coefs.R::getActualCohortTes`, and
  `R/variance_machinery.R::getCohortATTsFinal` for the internal-helper pattern
  in practice.
  **Exception: S3 methods registered for a generic in another package**
  (e.g., `print.fetwfe`, `summary.etwfe`, `coef.fetwfe`) carry only
  `#' @export` — the generic in the parent package owns the `@param` /
  `@return` contract, and a per-method block would only restate it. See
  `R/fetwfe_class.R` and `R/etwfe_class.R` for the existing pattern;
  match it.
- **Imports are declared at the top of `R/fetwfe.R`** via
  `@import` / `@importFrom`. If you add a new dependency, also update
  `Imports:` in DESCRIPTION.
- **S3 classes**: outputs of `fetwfe()` and `etwfe()` carry classes
  `"fetwfe"` / `"etwfe"`. Methods (`print`, `summary`, `coef`) live in
  `fetwfe_class.R` / `etwfe_class.R` and are registered via
  `S3method(...)` in `NAMESPACE` (roxygen handles this).
- **Verbosity gating**: use `message()` (not `print()` / `cat()`) for progress
  output, and gate it on a `verbose` argument. CRAN requires no chatter when
  `verbose = FALSE`.
- **No compiled code.** The package is pure R. Don't add C/C++ unless the user
  explicitly asks — it would change the CRAN submission profile.

## Domain terminology cheat sheet

- **N, T, R, d, p**: units, time periods, treated cohorts, covariates,
  total parameters. These names are used consistently across the codebase.
- **Cohort**: set of units that adopt treatment at the same time. The
  never-treated group is the reference cohort and must be non-empty in the
  final period.
- **Absorbing-state treatment**: once treated, always treated. Units treated
  in the first period are dropped automatically with a warning (see
  `idCohorts()` in `utility.R`).
- **`q`**: bridge penalty exponent. `q < 1` gives sparsity and valid standard
  errors; `q = 1` is lasso; `q = 2` is ridge. Default is `0.5`.
- **`indep_counts`**: optional independent-sample cohort counts. When supplied,
  the overall-ATT SE is asymptotically exact rather than conservative.
- **`sig_eps_sq`, `sig_eps_c_sq`**: observation-level and unit-level variance
  components. Estimated via Pesaran (2015, §26.5.1) ridge regression when
  not supplied.

## Build / test / check

This is a standard R package, so the usual commands apply. Run from the repo
root:

```r
# In R, with devtools / usethis available:
devtools::document()      # regenerate man/*.Rd and NAMESPACE from roxygen
devtools::test()          # run testthat suite
devtools::check()         # full R CMD check
devtools::check(args = "--as-cran")   # CRAN-strict check
devtools::build_vignettes()
```

Or from the shell:

```sh
R CMD build .
R CMD check --as-cran fetwfe_<version>.tar.gz
```

After **any** change to roxygen comments, function signatures, or exports,
run `devtools::document()` so `man/` and `NAMESPACE` stay in sync. Stale
generated files will fail `R CMD check`.

## Things that will break a CRAN submission

- Editing `man/*.Rd` or `NAMESPACE` by hand — they're regenerated.
- Examples that take longer than ~5 seconds (use `\donttest{}` for slow ones).
- Writing to anywhere other than `tempdir()` from examples / tests / vignettes.
- Console output from non-interactive code paths (gate behind `verbose`).
- Adding a dependency without updating DESCRIPTION's `Imports:` / `Suggests:`.
- Bumping `Version:` in DESCRIPTION without a matching entry in `NEWS.md` and
  `inst/CITATION`.
- Uncommitted CRAN-only artifacts (`fetwfe.Rcheck/`, `..Rcheck/`, `Meta/`,
  `doc/`, `paper_arxiv.tex`) leaking into the build — they're already in
  `.Rbuildignore`; keep it that way.

## Tests

`tests/testthat/` has one `test-<function>.R` per exported function. The
larger files (`test-fetwfe.R`, `test-etwfe.R`, `test-betwfe.R`,
`test-twfeCovs.R`) build small synthetic panels via the shared helpers
`generate_panel_data()`, `generate_bad_panel_data()`, and
`generate_minimal_panel_data()` defined in
`tests/testthat/helper-panel-fixture.R` (sourced by testthat before any
`test-*.R` runs). When you add or change a public function, add or
update its `test-*.R` file in the same PR.

## When in doubt

- Methodology questions → `paper_arxiv.tex` is authoritative; the code is
  derived from it.
- API behavior → the roxygen header on the user-facing function is
  authoritative; `man/*.Rd` is regenerated from it.
- Recent changes → `NEWS.md` and `git log`.

## Class constructor validators

Each S3-classed estimator object (`fetwfe`, `etwfe`, `betwfe`, `twfeCovs`) is
checked against a runtime `.validate_<class>(x)` helper called at the bottom of
the entry-point function before class assignment. The validators encode the
documented cross-slot contracts of each class — slot inventory, SE consistency
(if `calc_ses = FALSE` then `att_se` and `catt_ses` must all be NA), selection
consistency, p-value NA-derivation, catt_df shape, cohort-probability
structural sanity, dimensions, lambda monotonicity, type sanity. See issue #85
for the design rationale.

When adding a new public function that constructs an instance of one of these
classes, **call the appropriate `.validate_<class>(out)` helper immediately
before `class(out) <- "<class>"`**. The validator catches malformed output at
construction time rather than at downstream-method-confusion time.

When adding a new method that READS from a class object (`eventStudy`,
`augment`, `tidy`, `glance`, `plot`, `coef`, future `predict`, etc.), call the
matching `.check_for_<method>(x)` helper at the top of the method
(infrastructure added by issue #86 — see `R/class_helpers.R`). The helper
re-validates the object via `.assert_estimator_object(x)` and dispatches to
the right `.validate_<class>` internally based on `inherits()`. For methods
that need a derived contract (e.g., `.check_for_event_study(x)` returns
`list(has_valid_ses = ...)` to gate SE computation), the helper returns a
small named list; otherwise it returns `invisible(x)`. The pattern is what
fixed #73 (`eventStudy()` reporting finite SEs when the fit's `calc_ses`
was FALSE).

When adding a new slot to a class's `@return`:

1. Add the slot to `.EXPECTED_SLOTS_<CLASS>` in `R/<class>_class.R` (or
   `R/twfeCovs.R` for twfeCovs).
2. Add the corresponding `\item{name}{description}` to the `@return` block.
3. Add a contract assertion in the validator if there's a non-trivial
   invariant.
4. The `tests/testthat/test-doc-slot-parity.R` (#70) will catch docs-vs-code
   drift; the `tests/testthat/test-class-validators.R` cross-class consistency
   block will catch validator-vs-code drift.
