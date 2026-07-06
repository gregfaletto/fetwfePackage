# Dollar-sign assignment on a `catt_df` object

Intercepts assignment by the pre-1.11.0 Title-Case column names (e.g.,
`df$Cohort <- v`) and stops with a migration message. All other
assignment falls through to the `data.frame` method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x$name <- value
```

## Arguments

- x:

  A `catt_df` object.

- name:

  Character; the column name being assigned via `$`.

- value:

  The value to assign.

## Value

The modified `catt_df` object, as for `$<-.data.frame`.
