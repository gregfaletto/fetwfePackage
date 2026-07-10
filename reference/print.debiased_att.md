# Print a debiased overall-ATT estimate

Compact display of
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)'s
estimate, standard error, and Wald confidence interval. In the
high-dimensional (`p >= NT`) regime it also surfaces the
nodewise-direction diagnostics (`lambda_c` / `lambda_node`, KKT
feasibility, convergence) the documentation tells users to inspect; the
fixed-`p` print shows only the estimate block.

## Usage

``` r
# S3 method for class 'debiased_att'
print(x, ...)
```

## Arguments

- x:

  A `"debiased_att"` object from
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md).

- ...:

  Ignored.

## Value

`x`, invisibly.
