# Double-bracket access on a `catt_df` object

Intercepts access by the pre-1.11.0 Title-Case column names (`Cohort`,
`Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`) and stops
with a migration message pointing to the new snake_case name. All other
access falls through to the `data.frame` method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x[[i, ...]]
```

## Arguments

- x:

  A `catt_df` object (data frame with class
  `c("catt_df", "data.frame")`).

- i:

  Index; passed through to `[[.data.frame`. Character indices matching
  an old column name are intercepted; other indices fall through.

- ...:

  Further arguments passed through.

## Value

The column or value, as for `[[.data.frame`.
