# Tidy a debiased overall-ATT estimate

One-row
[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) data
frame for a
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
result, following the package's `broom` column convention (`term`,
`estimate`, `std.error`, `conf.low`, `conf.high`) – the same columns
[`tidy.fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.fetwfe.md)
/
[`tidy.eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.eventStudy.md)
return.

## Usage

``` r
# S3 method for class 'debiased_att'
tidy(x, ...)
```

## Arguments

- x:

  A `"debiased_att"` object from
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md).

- ...:

  Ignored.

## Value

A one-row data frame with columns `term`, `estimate`, `std.error`,
`conf.low`, `conf.high`.
