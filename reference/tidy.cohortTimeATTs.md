# Tidy a `cohortTimeATTs` object

Returns a `broom`-style tidy data frame for the output of
[`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md).
Renames the snake_case columns to broom conventions (`se` -\>
`std.error`, `p_value` -\> `p.value`, `ci_low` / `ci_high` -\>
`conf.low` / `conf.high`), keeps the `time` column, and adds a `term`
column (`"cohort_<cohort label>_time_<time>"`) plus a `statistic` column
(`estimate / std.error`) so the schema parallels
[`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md)
and
[`tidy.eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.eventStudy.md)
for downstream `bind_rows()` consumers. When the input carries a
`selected` column (`fetwfe` / `betwfe`), it is passed through as the
final column.

## Usage

``` r
# S3 method for class 'cohortTimeATTs'
tidy(x, ...)
```

## Arguments

- x:

  A `cohortTimeATTs` object returned by
  [`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md).

- ...:

  Unused; present for S3 compatibility.

## Value

A data frame with one row per `(cohort, time)` cell and columns `term`,
`time`, `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
`conf.high`, and (if present in the input) `selected`.

## Details

Confidence intervals are the pointwise `1 - alpha` Wald bounds
[`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md)
computed (encoding the `alpha` passed there); like
[`tidy.cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.cohortStudy.md)
this method passes them through rather than recomputing at a custom
`conf.level`.

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- fetwfeWithSimulatedData(
    simulateData(genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2),
                 N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  )
  broom::tidy(cohortTimeATTs(res))
} # }
```
