# Double-bracket assignment on a `catt_df` object

Intercepts assignment by the pre-1.11.0 Title-Case column names (e.g.,
`df[["Estimated TE"]] <- v`) and stops with a migration message pointing
to the new snake_case name. This closes the gap where a partial
migration (RHS updated to new name, LHS still old) would silently append
a new column rather than overwriting the renamed one. All other
assignment falls through to the `data.frame` method via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x[[i, ...]] <- value
```

## Arguments

- x:

  A `catt_df` object.

- i:

  Index; passed through to `[[<-.data.frame`. Character indices matching
  an old column name are intercepted.

- ...:

  Further arguments passed through.

- value:

  The value to assign.

## Value

The modified `catt_df` object, as for `[[<-.data.frame`.
