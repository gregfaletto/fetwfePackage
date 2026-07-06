# Augment is not defined for a twfeCovs fit (documented omission)

`augment()` is intentionally not provided for
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
objects (#58): the fit's coefficient vector lives in a reduced
cohort-level basis that the shared fitted-value path (`X %*% beta_hat`)
does not match, so a meaningful `.fitted` / `.resid` cannot be
reconstructed. Use
[`tidy.twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.twfeCovs.md),
[`glance.twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/glance.twfeCovs.md),
or [`summary()`](https://rdrr.io/r/base/summary.html) instead. Calling
this method always raises an error.

## Usage

``` r
# S3 method for class 'twfeCovs'
augment(x, data, ...)
```

## Arguments

- x:

  An object of class `"twfeCovs"`.

- data:

  Ignored.

- ...:

  Ignored.

## Value

(none; raises an error).
