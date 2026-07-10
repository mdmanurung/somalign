# Run an OT sensitivity grid

Run an OT sensitivity grid

## Usage

``` r
somalign_sensitivity_grid(
  query,
  reference,
  epsilon,
  rho_query,
  rho_ref,
  solver = c("internal", "auto"),
  parallel = FALSE,
  ...
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Numeric vector of entropic regularisation values.

- rho_query:

  Numeric vector of query-side mass relaxation values.

- rho_ref:

  Numeric vector of reference-side mass relaxation values.

- solver:

  Solver passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  `"auto"` is accepted as a compatibility alias for the internal pure-R
  solver.

- parallel:

  Logical. When `TRUE`, grid rows are evaluated in parallel using
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  with `mc.cores = getOption("mc.cores", 1L)`. On Windows `mclapply`
  falls back to a single core automatically. When `FALSE` (default) a
  sequential for-loop is used, which is fully reproducible across
  platforms.

- ...:

  Additional arguments passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

## Value

A data frame with one row per parameter combination.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
somalign_sensitivity_grid(qry, ref,
                          epsilon = c(0.05, 0.1),
                          rho_query = c(0.5, 1),
                          rho_ref = 1)
} # }
```
