# Convert data prepared for `etwfe::etwfe()` to the format required by `fetwfe()` and `fetwfe::etwfe()`

`etwfeToFetwfeDf()` reshapes and renames a panel dataset that is already
formatted for `etwfe::etwfe()` (McDermott 2024) so that it can be passed
directly to
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
or
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
from the `fetwfe` package. In particular, it

- creates an *absorbing-state* treatment dummy that equals 1 from the
  first treated period onward\* and 0 otherwise,

- (optionally) drops units that are already treated in the very first
  period of the sample (because
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  removes them internally), and

- returns a tidy dataframe whose column names match the arguments that
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)/[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  expect.

## Usage

``` r
etwfeToFetwfeDf(
  data,
  yvar,
  tvar,
  idvar,
  gvar,
  covars = character(0),
  drop_first_period_treated = TRUE,
  out_names = list(time = "time_var", unit = "unit_var", treatment = "treatment",
    response = "response"),
  verbose = FALSE
)
```

## Arguments

- data:

  A long-format data.frame that you could already feed to
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md).

- yvar:

  Character. Column name of the outcome (left-hand side in your `fml`).

- tvar:

  Character. Column name of the time variable that you pass to
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  as `tvar`.

- idvar:

  Character. Column name of the unit identifier (the variable you would
  cluster on, or pass to `etwfe(..., ivar = idvar)` if you were using
  unit FEs).

- gvar:

  Character. Column name of the "first treated" cohort variable passed
  to
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  as `gvar`. Must be `0` (or `Inf`, which is mapped to `0`) for
  never-treated units, or the (strictly positive) first treated period.

- covars:

  Character vector of *additional* covariate columns to keep (default
  `character(0)`).

- drop_first_period_treated:

  Logical. Should units already treated in the very first sample period
  be removed?
  ([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  will drop them internally anyway, but doing it here keeps the returned
  dataframe clean.) Default `TRUE`.

- out_names:

  Named list giving the column names that the returned dataframe should
  have. The default (`time_var`, `unit_var`, `treatment`, `response`)
  matches the arguments usually supplied to
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).
  **Do not change the *names* of this list** – only the *values* – and
  keep all four.

- verbose:

  Logical. If `TRUE`, a
  [`message()`](https://rdrr.io/r/base/message.html) reports the count
  of first-period-treated unit-period rows dropped when
  `drop_first_period_treated = TRUE`. Default `FALSE` (silent).

## Value

A tidy `data.frame` with (in this order)

- `time_var` integer,

- `unit_var` character,

- `treatment` integer 0/1 absorbing-state dummy,

- `response` numeric outcome,

- any covariates requested in `covars`. Ready to pass straight to
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  or
  [`fetwfe::etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md).

## References

McDermott G (2024). *etwfe: Extended Two-Way Fixed Effects*.
doi:10.32614/CRAN.package.etwfe
[doi:10.32614/CRAN.package.etwfe](https://doi.org/10.32614/CRAN.package.etwfe)
, R package version 0.5.0, <https://CRAN.R-project.org/package=etwfe>.

## Examples

``` r
## toy example ---------------------------------------------------------------
if (FALSE) { # \dontrun{
library(did)  # provides the mpdta example dataframe
data(mpdta)

head(mpdta)

tidy_df <- etwfeToFetwfeDf(
  data  = mpdta,
  yvar = "lemp",
  tvar = "year",
  idvar = "countyreal",
  gvar = "first.treat",
  covars = c("lpop"))

head(tidy_df)

} # }
## Now you can call fetwfe()  ------------------------------------------------
# res <- fetwfe(
#   pdata      = tidy_df,
#   time_var   = "time_var",
#   unit_var   = "unit_var",
#   treatment  = "treatment",
#   response   = "response",
#   covs       = c("lpop"))
```
