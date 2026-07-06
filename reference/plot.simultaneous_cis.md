# Plot a simultaneous-confidence-interval object

Draws the per-effect estimates with both bands overlaid – the wider
simultaneous band (error bars) on top of the narrower pointwise band
(thick line range) – for a
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
result, so the family-wise vs per-effect coverage trade-off is visible.

## Usage

``` r
# S3 method for class 'simultaneous_cis'
plot(x, ...)
```

## Arguments

- x:

  A `"simultaneous_cis"` object from
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md).

- ...:

  Ignored.

## Value

A `ggplot` object.

## Note

Requires the `ggplot2` package (in Suggests). Install via
`install.packages("ggplot2")` if it is not already installed. Matches
`R/plot.R`'s `.plot_estimator` precedent for ggplot2-dependent plot
methods.
