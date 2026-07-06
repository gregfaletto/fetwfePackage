# Single-bracket access on a `catt_df` object

Intercepts column selection by the pre-1.11.0 Title-Case names and stops
with a migration message. The check fires when an old name appears in
the column-selector position: `j` for the `df[i, j]` two-index form, or
`i` for the `df[j]` one-index column-selection form. Row-only access
(`df[i, ]`) and access by new column names fall through to the
`data.frame` method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x[i, j, ...]
```

## Arguments

- x:

  A `catt_df` object.

- i, j:

  Row / column selectors; see `[.data.frame`.

- ...:

  Further arguments passed through.

## Value

The selected subset, as for `[.data.frame`.
