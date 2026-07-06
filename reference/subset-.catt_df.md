# Single-bracket assignment on a `catt_df` object

Intercepts column assignment by the pre-1.11.0 Title-Case names and
stops with a migration message. The check fires when an old name appears
in the column-selector position. Row-only assignment (`df[i, ] <- v`)
and assignment by new column names fall through to the `data.frame`
method via [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Usage

``` r
# S3 method for class 'catt_df'
x[i, j, ...] <- value
```

## Arguments

- x:

  A `catt_df` object.

- i, j:

  Row / column selectors; see `[<-.data.frame`.

- ...:

  Further arguments passed through.

- value:

  The value to assign.

## Value

The modified `catt_df` object, as for `[<-.data.frame`.

## Details

The [`nargs()`](https://rdrr.io/r/base/nargs.html) distinction here
mirrors the read-side `[.catt_df`, but the threshold shifts by one
because `[<-` carries an extra positional `value` argument: `df[i] <- v`
has `nargs() == 3` and i is the column selector; `df[i, j] <- v` and
`df[i, ] <- v` have `nargs() == 4` and j (if non-missing) is the column
selector.
