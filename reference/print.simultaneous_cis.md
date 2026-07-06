# Print a simultaneous-confidence-interval object

Compact display of a
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
result: the family, `K`, and `alpha`; the simultaneous / pointwise /
Bonferroni critical values; the high-dimensional desparsified-band
diagnostics (only in the `p >= NT` regime); and the per-effect interval
table.

## Usage

``` r
# S3 method for class 'simultaneous_cis'
print(x, ...)
```

## Arguments

- x:

  A `"simultaneous_cis"` object from
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md).

- ...:

  Ignored.

## Value

`x`, invisibly.
