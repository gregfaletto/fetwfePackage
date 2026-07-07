# Coming from did or etwfe? Bringing your data into fetwfe

The two most widely used R packages for difference-in-differences with
staggered adoptions are [**`did`**](https://bcallaway11.github.io/did/)
(Callaway and Sant’Anna 2021) and
[**`etwfe`**](https://grantmcdermott.com/etwfe/) (Wooldridge 2021;
implemented by Grant McDermott). If your panel is already set up for
either, `fetwfe` ships one-line converters, so you can try the fused
estimator without re-shaping your data:

- [`attgtToFetwfeDf()`](https://gregfaletto.github.io/fetwfePackage/reference/attgtToFetwfeDf.md)
  — for a panel formatted for
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html),
- [`etwfeToFetwfeDf()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfeToFetwfeDf.md)
  — for a panel formatted for `etwfe::etwfe()`.

Both take a long panel plus the column names you already use; each
returns a data frame with the `time_var` / `unit_var` / `treatment` /
`response` columns that
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
expects (building the absorbing-state treatment indicator, dropping any
units already treated in the first period, and renaming the columns for
you).

*(The runnable examples below use the `did` package’s `mpdta` dataset;
install `did` to reproduce them.)*

## From `did`

[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
identifies a panel by outcome (`yname`), time (`tname`), unit
(`idname`), and the **first-treated cohort** (`gname`, `0` for
never-treated units). The canonical example is `mpdta` — county-level
teen employment and the minimum wage:

``` r

library(did)
data(mpdta)
head(mpdta[, c("countyreal", "year", "first.treat", "lemp")])
#>     countyreal year first.treat     lemp
#> 866       8001 2003        2007 8.461469
#> 841       8001 2004        2007 8.336870
#> 842       8001 2005        2007 8.340217
#> 819       8001 2006        2007 8.378161
#> 827       8001 2007        2007 8.487352
#> 937       8019 2003        2007 4.997212
```

One line converts it, and then you fit
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
as usual:

``` r

library(fetwfe)

fdf <- attgtToFetwfeDf(
    mpdta,
    yname = "lemp",
    tname = "year",
    idname = "countyreal",
    gname = "first.treat"
)

res <- fetwfe(
    pdata = as.data.frame(fdf),
    time_var = "time_var",
    unit_var = "unit_var",
    treatment = "treatment",
    response = "response"
)

round(c(ATT = res$att_hat, SE = res$att_se), 4)
#>     ATT      SE 
#> -0.0387  0.0123
```

## From `etwfe`

`etwfe::etwfe()` uses the same long-panel shape, identified by `yvar` /
`tvar` / `idvar` / `gvar` — where `gvar` is again the first-treated
period (`0` for never-treated).
[`etwfeToFetwfeDf()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfeToFetwfeDf.md)
is the parallel converter. Since `mpdta` is already in that shape
(`first.treat` is the cohort variable), it converts identically:

``` r

edf <- etwfeToFetwfeDf(
    mpdta,
    yvar = "lemp",
    tvar = "year",
    idvar = "countyreal",
    gvar = "first.treat"
)

# Byte-identical to the fetwfe-ready panel from the did converter above:
identical(edf, fdf)
#> [1] TRUE
```

## What `fetwfe` adds

`did` and `etwfe` both return the full set of cohort-by-time
(group-time) treatment effects. FETWFE starts from that *same* saturated
extended-two-way-fixed-effects model, and then:

- **fuses** effects that are statistically similar across neighboring
  cohorts and time periods — a data-driven bias–variance trade-off that
  spends fewer effective parameters when the data support pooling, and
  retains heterogeneity when they don’t; and
- returns an **overall ATT with asymptotically valid standard errors and
  confidence intervals *after* that selection**, along with per-cohort
  and event-study breakdowns
  ([`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md),
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md))
  and family-wise simultaneous confidence bands
  ([`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)).

So the converters let you keep your existing data pipeline and add the
fused estimate alongside your `did` / `etwfe` results. For the full
workflow see the introductory vignette
([`vignette("fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/fetwfe.md)),
and for the choice of fusion geometry see *“Choosing a fusion structure:
cohort vs. event-study penalties”*
([`vignette("fusion_structure_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/fusion_structure_vignette.md)).

## References

- Callaway, B. and Sant’Anna, P. H. C. (2021). Difference-in-Differences
  with multiple time periods. *Journal of Econometrics*, 225(2),
  200–230.
- Wooldridge, J. M. (2021). Two-way fixed effects, the two-way Mundlak
  regression, and difference-in-differences estimators. SSRN Working
  Paper No. 3906345.
- Faletto, G. (2025). *Fused Extended Two-Way Fixed Effects for
  Difference-in-Differences with Staggered Adoptions*.
  [arXiv:2312.05985](https://arxiv.org/abs/2312.05985).
