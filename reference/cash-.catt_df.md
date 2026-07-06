# Dollar-sign access on a `catt_df` object

Intercepts access by the pre-1.11.0 Title-Case column names and stops
with a migration message. All other access falls through to the
`data.frame` method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x$name
```

## Arguments

- x:

  A `catt_df` object.

- name:

  Character; the column name being accessed via `$`.

## Value

The column, as for `$.data.frame`.
