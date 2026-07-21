# Convert data formatted for `att_gt()` to a dataframe suitable for `fetwfe()` / `etwfe()`

`attgtToFetwfeDf()` reshapes and renames a panel dataset that is already
formatted for
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
(Callaway and Sant'Anna 2021) so that it can be passed directly to
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
attgtToFetwfeDf(
  data,
  yname,
  tname,
  idname,
  gname,
  covars = character(0),
  drop_first_period_treated = TRUE,
  out_names = list(time = "time_var", unit = "unit_var", treatment = "treatment",
    response = "response"),
  verbose = FALSE
)
```

## Arguments

- data:

  A `data.frame` in **long** format containing at least the four columns
  used by
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html):
  outcome `yname`, time `tname`, unit id `idname`, and the
  first-treatment period `gname` (which is 0 for the never-treated
  group).

- yname:

  Character scalar. Name of the outcome column.

- tname:

  Character scalar. Name of the time variable (numeric or integer). This
  becomes `time_var` in the returned dataframe.

- idname:

  Character scalar. Name of the unit identifier. Converted to character
  and returned as `unit_var`.

- gname:

  Character scalar. Name of the *group* variable holding the first
  period of treatment. Values must be 0 (or `Inf`, a common
  never-treated encoding, which is mapped to 0) for never-treated, or a
  positive integer representing the first treated period.

- covars:

  Character vector of additional covariate column names to carry through
  (default `character(0)`). These columns are left untouched and appear
  *after* the required columns in the returned dataframe.

- drop_first_period_treated:

  Logical. If `TRUE` (default), units that are already treated in the
  first sample period are removed *before* creating the treatment dummy.
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  would do this internally, but dropping them here keeps the returned
  dataframe cleaner.

- out_names:

  A named list giving the column names to use in the resulting
  dataframe. Defaults are
  `list(time = "time_var", unit = "unit_var", treatment = "treatment", response = "response")`.
  Override if you prefer different names (for instance, to keep the
  original `yname`). The vector *must* contain exactly these four names.

- verbose:

  Logical. If `TRUE`, a
  [`message()`](https://rdrr.io/r/base/message.html) reports the count
  of first-period-treated unit-period rows dropped when
  `drop_first_period_treated = TRUE`. Default `FALSE` (silent).

## Value

A `data.frame` with columns `time_var`, `unit_var`, `treatment`,
`response`, and any covariates requested in `covars`, ready to be fed to
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)/[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md).
All required columns are of the correct type: `time_var` is integer,
`unit_var` is character, `treatment` is integer 0/1, and `response` is
numeric.

## References

Callaway, Brantly and Pedro H.C. Sant'Anna. "Difference-in- Differences
with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2,
pp. 200-230, 2021.
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)
, <https://arxiv.org/abs/1803.09015>.

## Examples

``` r
## toy example ---------------------------------------------------------------
if (FALSE) { # \dontrun{
library(did)  # provides the mpdta example dataframe
data(mpdta)

head(mpdta)

tidy_df <- attgtToFetwfeDf(
  data  = mpdta,
  yname = "lemp",
  tname = "year",
  idname = "countyreal",
  gname = "first.treat",
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
