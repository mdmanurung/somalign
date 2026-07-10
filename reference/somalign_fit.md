# Align a query SOM to a reference SOM

Align a query SOM to a reference SOM

## Usage

``` r
somalign_fit(
  query,
  reference,
  epsilon = 0.05,
  rho_query = 1,
  rho_ref = 1,
  solver = c("auto", "pot", "internal"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Entropic regularisation strength.

- rho_query:

  Query-side unbalanced mass relaxation.

- rho_ref:

  Reference-side unbalanced mass relaxation.

- solver:

  `"auto"`, `"pot"`, or `"internal"`.

- min_match_fraction:

  Minimum transported fraction required before a query node label
  transfer is accepted.

- confidence_threshold:

  Minimum top-label probability required before a query node label
  transfer is accepted.

- correction_min_mass:

  Minimum transported node mass required before a correction shift is
  applied. Corrections also require the node match fraction to pass
  `min_match_fraction`.

- max_iter:

  Maximum internal Sinkhorn iterations.

- tol:

  Internal Sinkhorn convergence tolerance.

- chunk_size:

  Integer. Number of samples to project per chunk when computing nearest
  reference node. Use `Inf` or `NULL` for no chunking (allocates a full
  n_samples x n_nodes matrix). Default `10000L`.

## Value

A `somalign_fit` object.

## Details

The transport plan row sums will not equal `query$node_masses` exactly —
this is by design. Unbalanced optimal transport allows mass destruction,
so some query mass may be absorbed rather than transported. Deviation
grows with lower `rho_query` / `rho_ref` values and higher `epsilon`. At
the defaults (`rho_query = 1`, `rho_ref = 1`, `epsilon = 0.05`), row-sum
deviation can reach approximately 13%. Use
`diagnostics$ot$max_row_mass_error` to quantify the deviation in a given
fit; for near-balanced data, increase `rho_query` (e.g.
`rho_query = 10`) to enforce tighter marginal constraints.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
fit <- somalign_fit(qry, ref)
} # }
```
