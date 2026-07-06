# Tidy a `cohortStudy` object

Returns a `broom`-style tidy data frame for the output of
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md).
Renames the snake_case columns of `catt_df` to broom conventions (`se`
-\> `std.error`, `p_value` -\> `p.value`, `ci_low` / `ci_high` -\>
`conf.low` / `conf.high`) and adds a `term` column
(`"cohort_<cohort label>"`) plus a `statistic` column
(`estimate / std.error`) so the schema matches
[`tidy.eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.eventStudy.md)
for downstream `bind_rows()` consumers. When the input carries a
`selected` column (`fetwfe` / `betwfe`), it is passed through as the
final column.

## Usage

``` r
# S3 method for class 'cohortStudy'
tidy(x, ...)
```

## Arguments

- x:

  A `cohortStudy` object returned by
  [`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md).

- ...:

  Unused; present for S3 compatibility.

## Value

A data frame with one row per treated cohort and columns `term`,
`estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
`conf.high`, and (if present in the input) `selected`.

## Details

Confidence intervals come from the cohort fit's stored bounds (which
encode the alpha passed at fit time); unlike
[`tidy.eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.eventStudy.md),
this method does not recompute the CIs at a custom `conf.level` because
the standard errors in `catt_df` are already paired with the fit-time
bounds (`ci_low` / `ci_high`), so re-emitting those is the
minimum-surprise behavior.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- fetwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(cohortStudy(res))
} # }
```
