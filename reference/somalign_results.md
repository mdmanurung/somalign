# Return per-sample somalign results

Direct reference projection columns are canonical. Corrected projection
columns are auxiliary and should be used for visualisation, annotation,
and triage rather than feature-level differential testing.

## Usage

``` r
somalign_results(fit, data = NULL)
```

## Arguments

- fit:

  A `somalign_fit` object.

- data:

  Optional data frame to append after result columns.

## Value

A data frame with direct and corrected projection columns.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
fit <- somalign_fit(qry, ref)
somalign_results(fit)
} # }
```
