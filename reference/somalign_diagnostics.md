# Extract somalign diagnostics

Extract somalign diagnostics

## Usage

``` r
somalign_diagnostics(fit)
```

## Arguments

- fit:

  A `somalign_fit` object.

## Value

A named list of solver, OT, node, and projection diagnostics.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
fit <- somalign_fit(qry, ref)
somalign_diagnostics(fit)
} # }
```
