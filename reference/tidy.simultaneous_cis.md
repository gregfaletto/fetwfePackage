# Tidy a simultaneous-confidence-interval object

One row per effect for a
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
result, in `broom` form: the simultaneous band as `conf.low` /
`conf.high`, plus the pointwise band as `pointwise.conf.low` /
`pointwise.conf.high`.

## Usage

``` r
# S3 method for class 'simultaneous_cis'
tidy(x, ...)
```

## Arguments

- x:

  A `"simultaneous_cis"` object from
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md).

- ...:

  Ignored.

## Value

A data frame with one row per effect and columns `term`, `estimate`,
`conf.low`, `conf.high`, `pointwise.conf.low`, and
`pointwise.conf.high`.
