# Simulated Panel-Data Class

S3 class for objects returned by
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md).
Compact `print` method summarizes the panel's dimensions and cohort
structure instead of dumping the full `N*T x p` design matrix (which the
default `print.list` would do).
