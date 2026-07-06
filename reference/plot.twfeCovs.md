# Plot is not defined for a twfeCovs fit (documented omission)

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is
intentionally not provided for
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
objects (#58):
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
estimates one pooled effect per cohort, so there is no per-(cohort,
time) / event-study structure to plot. Use
[`summary()`](https://rdrr.io/r/base/summary.html) or
[`tidy.twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/tidy.twfeCovs.md)
for the cohort effects. Calling this method always raises an error.

## Usage

``` r
# S3 method for class 'twfeCovs'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"twfeCovs"`.

- ...:

  Ignored.

## Value

(none; raises an error).
